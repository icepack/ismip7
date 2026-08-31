"""thermo_model.py - depth-averaged enthalpy (Stefan) model for a physical
fluidity prior mean.

Ported from mismip_time-dependent-da/scripts/thermo_model.py (Recinos/fenics_ice
line of work) into the ISMIP7 tree, with ONE adaptation for real Antarctica:
the surface-energy boundary condition ``E_srf`` may be a spatially-varying
Function (a mean-annual surface-temperature field) rather than the single
scalar the idealized MISMIP+ setup used. Everything else is the original
physics, verbatim.

Model per solve (steady, linear given lagged inputs; E in DG0 on the plan-view
mesh):
    upwind DG0 transport of h*E  -  D*E  (discrete-divergence correction so the
        operator annihilates constants exactly; D = per-cell facet divergence
        of the CURRENT (h, u))
  + kappa*alpha*(E - E_srf)/h            (HeatTransport3D-style surface exchange)
  + a_plus*(E - E_srf)                   (accumulation deposits ice at surface
                                          energy; basal melt cancels)
  = h*q_strain + q_fric + q_geo
with q_strain = 2 A^(-1/n) eps_e^(1/n+1), q_fric = grounded * C|u|^(1/m+1), and
the Stefan relations T = min(E/rho c, Tm), w = (E - rho c Tm)+ / rho L.

Shear-layer closure: solve once with and once without strain heating (the
operator is linear), amplify only the strain anomaly:
    E_eff = E_base + shear_amp*(E_full - E_base)
    A = rate_factor(T(E_eff)) * (1 + duval*min(w(E_eff), w_max)).
All nonlinear couplings (A <-> u <-> heating) are Picard-lagged by the caller.
"""
import numpy as np
import firedrake as fd
from firedrake import (Constant, inner, grad, dx, dS, ds, sqrt, sym,
                       max_value, min_value, avg, jump)
from icepack.constants import (ice_density as rho_I, heat_capacity as c_heat,
                               thermal_diffusivity as alpha_th, latent_heat as L_heat,
                               melting_temperature as Tm, glen_flow_law as n_glen)
from icepack.models.viscosity import rate_factor

DEFAULTS = dict(kappa=4.0, T_srf=263.15, q_geo=50.0, shear_amp=30.0,
                duval=181.25, w_max=0.01, friction_exp=3.0, melt_delta=50.0,
                u_scale=1.0)


# 1 mW/m^2 = 3.15576e4 Pa m/yr = 3.15576e-2 MPa m/yr (icepack MPa-m-yr units)
def q_geo_icepack(q_geo_mW):
    return q_geo_mW * 3.15576e4 / 1e6


def energy_surface(p, T_srf_field=None):
    r"""Surface enthalpy E_srf = rho_I c T_srf. Pass ``T_srf_field`` (a Function
    or UFL expr, K) for a spatially-varying surface BC; otherwise the scalar
    ``p['T_srf']`` is used (the idealized-MISMIP default)."""
    T = p["T_srf"] if T_srf_field is None else T_srf_field
    return float(rho_I * c_heat) * T


def temperature(E):
    return min_value(E / (rho_I * c_heat), Tm)


def meltwater(E):
    return max_value(E - rho_I * c_heat * Tm, 0.0) / (rho_I * L_heat)


def fluidity_from(E_eff, p):
    A_cold = rate_factor(temperature(E_eff))
    w = min_value(meltwater(E_eff), Constant(p["w_max"]))
    return A_cold * (1 + Constant(p["duval"]) * w)


def grounded_frac(h, s, bed, p):
    cavity = (s - h) - bed
    return 1.0 - 0.5 * (1.0 + fd.tanh(cavity / Constant(p["melt_delta"])))


def solve_energy(u_adv, h, s, C, A_k, bed, acc, p, with_strain=True,
                 T_srf_field=None):
    """One steady linear enthalpy solve at the CURRENT fields. A_k, acc: DG0.

    ``T_srf_field`` (Function/UFL, K) sets a spatially-varying surface BC; None
    falls back to the scalar ``p['T_srf']``."""
    mesh = h.function_space().mesh()
    QD = A_k.function_space()
    u = Constant(p["u_scale"]) * u_adv
    E = fd.TrialFunction(QD)
    psi = fd.TestFunction(QD)
    nu = fd.FacetNormal(mesh)
    un = 0.5 * (inner(u, nu) + abs(inner(u, nu)))
    E_srf = energy_surface(p, T_srf_field)
    if T_srf_field is not None:
        E_srf = fd.Function(QD).project(E_srf)   # materialize the field once
    else:
        E_srf = Constant(E_srf)

    a_adv = (jump(psi) * (un('+') * avg(h) * E('+') - un('-') * avg(h) * E('-'))) * dS \
        + psi * un * h * E * ds
    L_in = -psi * (inner(u, nu) - un) * h * E_srf * ds
    # discrete-divergence correction (operator must annihilate constants exactly)
    F_div = (jump(psi) * avg(h) * (un('+') - un('-'))) * dS + psi * h * inner(u, nu) * ds
    D_fn = fd.Function(QD)
    cellvol = fd.assemble(fd.TestFunction(QD) * dx).dat.data_ro
    D_fn.dat.data[:] = fd.assemble(F_div).dat.data_ro / cellvol
    a_adv = a_adv - D_fn * E * psi * dx

    a_plus = max_value(acc, Constant(0.0))
    a_exch = (Constant(p["kappa"]) * alpha_th / h + a_plus) * E * psi * dx
    L_exch = (Constant(p["kappa"]) * alpha_th / h + a_plus) * E_srf * psi * dx

    eps = sym(grad(u_adv))
    eps_e = sqrt(0.5 * inner(eps, eps) + Constant(1e-16))
    q_strain = 2 * A_k ** (-1.0 / n_glen) * eps_e ** (1.0 / n_glen + 1)
    speed = sqrt(inner(u_adv, u_adv) + Constant(1e-12))
    q_fric = grounded_frac(h, s, bed, p) * C * speed ** (1.0 / p["friction_exp"] + 1)
    strain_on = Constant(1.0 if with_strain else 0.0)
    L_src = (strain_on * h * q_strain + q_fric
             + Constant(q_geo_icepack(p["q_geo"]))) * psi * dx

    E_out = fd.Function(QD, name="energy")
    fd.solve(a_adv + a_exch == L_in + L_exch + L_src, E_out,
             solver_parameters={"ksp_type": "preonly", "pc_type": "lu",
                                "pc_factor_mat_solver_type": "mumps"})
    return E_out


def thermo_update(u, h, s, C, A_k_dg, bed, acc_dg, p, T_srf_field=None):
    """Full thermo step: two linear solves + shear-layer superposition.
    Returns (E_eff, E_full) as DG0 Functions; caller maps fluidity_from(E_eff, p)."""
    E_full = solve_energy(u, h, s, C, A_k_dg, bed, acc_dg, p, with_strain=True,
                          T_srf_field=T_srf_field)
    E_base = solve_energy(u, h, s, C, A_k_dg, bed, acc_dg, p, with_strain=False,
                          T_srf_field=T_srf_field)
    E_eff = fd.Function(E_full.function_space(), name="energy_eff")
    E_eff.dat.data[:] = (E_base.dat.data_ro
                         + p["shear_amp"] * (E_full.dat.data_ro - E_base.dat.data_ro))
    return E_eff, E_full


def compute_fluidity_prior(u, h, s, bed, C, acc, T_srf, p=None, max_picard=60,
                           rtol=1e-3, A_init=20.0, A_clip=(1.0, 1e4),
                           h_min=10.0, relax=0.4, verbose=True, comm=None):
    r"""Physical fluidity prior mean ``A_prior(x)`` on ``h``'s CG space, from a
    FIXED-velocity thermomechanical Picard loop (A <-> heating at the observed
    u/geometry). This is the real-Antarctica analogue of MISMIP+
    ``thermo_fluidity.py``'s ``--fixed_u`` calibration; there is no A->u
    feedback, so it is the reference fluidity at the observations.

    Parameters
    ----------
    u, h, s, bed, C, acc : Function
        Velocity (obs), thickness, surface, bed, friction coefficient
        (balance/Weertman anchor), accumulation. ``h`` is floored at ``h_min``
        for the thermal operator so the ice-free buffer (h->0, where 1/h in the
        surface-exchange term would blow up) stays finite; buffer fluidity is
        irrelevant (no ice) and pinned to the surface value there.
    T_srf : Function or float
        Mean-annual surface temperature [K] (the englacial upper BC); a field
        for real Antarctica, a scalar for an idealized test.
    p : dict or None
        Overrides for ``DEFAULTS`` (kappa, q_geo, shear_amp, duval, ...).

    Returns
    -------
    A_prior : Function
        Fluidity on ``h.function_space()`` (icepack MPa^-n yr^-1 units), the
        prior mean for the inversion control ``phi = log(A / A_prior)``.
    """
    from firedrake.petsc import PETSc
    _print = (PETSc.Sys.Print if verbose else (lambda *a, **k: None))
    mesh = h.function_space().mesh()
    Q = h.function_space()
    QD = fd.FunctionSpace(mesh, "DG", 0)
    P = dict(DEFAULTS)
    if p:
        P.update(p)

    h_th = fd.Function(Q).interpolate(max_value(h, Constant(h_min)))
    acc_dg = fd.Function(QD).project(acc)
    if isinstance(T_srf, (int, float)):
        T_field = Constant(float(T_srf))
    else:
        T_field = fd.Function(QD).project(T_srf)

    A_k = fd.Function(Q, name="fluidity").interpolate(Constant(A_init))  # cold uniform start
    E = None
    _print(f"  Thermo fluidity prior (fixed u; kappa={P['kappa']}, q_geo={P['q_geo']} "
           f"mW/m2, shear_amp={P['shear_amp']}, duval={P['duval']}):")
    for it in range(1, max_picard + 1):
        A_dg = fd.Function(QD).project(A_k)
        E, _ = thermo_update(u, h_th, s, C, A_dg, bed, acc_dg, P, T_srf_field=T_field)
        A_new = fd.Function(Q, name="fluidity").project(fluidity_from(E, P))
        A_new.dat.data[:] = np.clip(A_new.dat.data_ro, *A_clip)
        # Under-relaxation: A <-> strain-heating is a stiff fixed point at
        # fixed velocity (higher A -> less shear heat -> colder -> lower A),
        # so a plain Picard iterate flip-flops. Damping with relax in (0,1]
        # converges it. A_(k+1) = A_k + relax*(A_new - A_k).
        A_relaxed = A_k.dat.data_ro + relax * (A_new.dat.data_ro - A_k.dat.data_ro)
        dloc = (np.abs(A_relaxed - A_k.dat.data_ro).max()
                / max(np.abs(A_k.dat.data_ro).max(), 1.0)) if A_k.dat.data_ro.size else 0.0
        dA = mesh.comm.allreduce(float(dloc), op=__import__("mpi4py").MPI.MAX)
        A_k.dat.data[:] = A_relaxed
        lo = mesh.comm.allreduce(float(A_k.dat.data_ro.min() if A_k.dat.data_ro.size else 1e30),
                                 op=__import__("mpi4py").MPI.MIN)
        hi = mesh.comm.allreduce(float(A_k.dat.data_ro.max() if A_k.dat.data_ro.size else -1e30),
                                 op=__import__("mpi4py").MPI.MAX)
        _print(f"    it {it:2d}: A [{lo:6.1f}, {hi:8.1f}]  dA={dA:.2e}")
        if dA < rtol:
            _print("    thermo prior converged.")
            break
    else:
        PETSc.Sys.Print(
            f"    WARNING: thermo fluidity prior did NOT converge after "
            f"{max_picard} Picard iterations (last dA={dA:.2e} >= rtol={rtol:.1e}); "
            f"using the partially-converged A_prior."
        )
    return A_k

