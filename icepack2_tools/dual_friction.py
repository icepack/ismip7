r"""Regularized-Coulomb (Joughin/Schoof) basal friction in the icepack2 DUAL,
written in VARIATIONAL RESIDUAL form.

Why a residual instead of the usual ``model.minimization.friction_power``
action: regularized Coulomb caps the basal stress at ``tau_c = c0 * N`` (N =
effective pressure), which is *not* a power law and so cannot be expressed by
``friction_power``.  Closing the basal stress directly as ``tau = tau_RC(u)``
makes the basal-stress block of the dual system the IDENTITY -- non-singular at
``tau = 0`` -- where the power-law dual is rank-deficient.  Two consequences the
power-law (Budd/Weertman) dual cannot deliver:

  * EXACT-zero basal drag on floating ice.  The Coulomb cap ``tau_c = c0 * N``
    goes to zero at flotation (N -> 0), so ``tau_b -> 0`` *exactly* -- no
    ``phi_eff`` floor, no ~1.7% residual shelf drag.  ``N`` is built from the
    *model* surface, so it vanishes at precisely the hydrostatic flotation
    criterion: genuinely frictionless shelves.

  * GROUNDED-ONLY friction inference.  The Weertman branch coefficient carries
    ``exp(theta * He)`` with ``He`` the smooth grounded indicator, so the
    adjoint gradient ``dJ/dtheta`` is identically zero on floating ice -- the
    optimizer physically cannot place friction on a shelf.  Friction is inferred
    where grounded; the shelves are masked by ``N -> 0``.

Together these make the diagnostic solve robust at a frictionless, migrating
grounding line -- the property that makes prognostic experiments tractable.

The viscous (membrane-stress) block is the SAME Goldsby-Kohlstedt composite that
``antarctica/scripts/inversion_icepack2.py`` uses: dislocation creep
(exponent ``n_flow``) plus an ``alpha``-weighted linear (n=1) regularizer
evaluated at a CONSTANT reference thickness ``H_ref`` so the block stays
positive-definite as ``h -> 0`` at the calving front.  Momentum balance and the
calving-front back-pressure are the canonical ``icepack2.model.variational``
forms (verified byte-equivalent to the hand-rolled versions they replace).

Provenance: ports ``gia-icepack/scripts/ase_model.py:build_F_rc`` (itself
derived from this repo's composite recipe) back into ISMIP7, swapping in the
n_flow=4 / H_ref composite viscous block.  Companion to
:mod:`icepack2_tools.grounding`.
"""

from firedrake import (
    Constant,
    Function,
    TestFunction,
    FacetNormal,
    split,
    sym,
    grad,
    tr,
    inner,
    sqrt,
    max_value,
    min_value,
    conditional,
    gt,
    eq,
    exp,
    avg,
    jump,
    dx,
    dS,
)
from icepack2 import model
from icepack2.constants import (
    ice_density as rho_I,
    water_density as rho_W,
    gravity as g,
)

from .geometry import surface_slope
from .grounding import height_above_flotation, smooth_heaviside

# Grounding-zone Heaviside width [m height-above-flotation].  Wider => more
# sub-element smoothing of the grounded/floating transition that gates theta.
GL_WIDTH = 10.0


def effective_pressure(H, s):
    r"""Ocean-connected effective pressure ``N = max(p_I - p_W, 0)`` [MPa].

    ``p_I = rho_I g H``, ``p_W = rho_W g max(0, H - s)``.  Uses the *model*
    surface ``s`` (hydrostatic on the shelf), so ``N`` vanishes EXACTLY at
    flotation -- the property that gives exact-zero Coulomb drag on shelves.
    A 1 m floor on H keeps it finite where the BedMachine thickness is 0.
    """
    Hs = max_value(H, Constant(1.0))
    p_I = rho_I * g * Hs
    p_W = rho_W * g * max_value(Constant(0.0), H - s)
    return max_value(p_I - p_W, Constant(0.0))


def grounded_mask(H, b, gl_width=GL_WIDTH):
    r"""Smooth grounded indicator ``He in [0, 1]`` (1 grounded, 0 floating) from
    height-above-flotation, smoothed over ``gl_width`` m HAF.  Gates the friction
    control ``theta`` to grounded ice (``dJ/dtheta -> 0`` as ``He -> 0``)."""
    haf = height_above_flotation(H, b)
    return smooth_heaviside(haf, kH=1.0 / gl_width)


def weertman_anchor(H, s, u_obs, m_slide, Q):
    r"""Balance Weertman coefficient ``C_w0 = tau_d / |u_obs|^(1/m)`` as a
    ``Function`` on ``Q`` -- a fixed anchor so the inverted ``theta`` is an O(1)
    log-adjustment rather than carrying the full friction magnitude.

    ``tau_d = rho_I g H |grad s|`` is the driving stress; ``|u_obs|`` is floored
    at 1 m/yr.  At ``u = u_obs`` and ``theta = 0`` the Weertman branch returns
    ``tau_W = tau_d`` -- i.e. the balance is driving-stress-consistent.

    ``s`` may be DG0 (see :func:`surface_slope`).
    """
    grad_s = surface_slope(s)
    tau_d = rho_I * g * H * sqrt(inner(grad_s, grad_s) + Constant(1e-12))
    u_speed = max_value(sqrt(inner(u_obs, u_obs)), Constant(1.0))
    return Function(Q, name="C_w0").interpolate(tau_d / u_speed ** (1.0 / m_slide))


def build_rc_residual(
    z,
    theta,
    phi,
    *,
    H,
    s,
    b,
    C_w0,
    A4_base,
    n_flow,
    n_flow_val,
    m_slide,
    tau_c,
    alpha,
    H_ref,
    fric_law="regularized_coulomb",
    N_ref=None,
    nhat_floor=0.0,
    nhat_cap=3.0,
    alpha_gl=0.0,
    c0=0.5,
    u_min=1.0,
    eps_tauc=0.0,
    c_w0_floor=0.0,
    h_visc_floor=0.0,
    ocean_drag=0.0,
    h_ocean=10.0,
    u_lim=0.0,
    k_lim=1e-3,
    gl_width=GL_WIDTH,
    calving_ids=None,
):
    r"""Assemble the icepack2 dual regularized-Coulomb residual ``F`` for the
    mixed state ``z = (u, M, tau)`` on ``Z = V x Sigma x T``.

    Controls
        theta : log Weertman/Coulomb-coefficient adjustment, gated to grounded
                ice via ``exp(theta * He)`` -- inferred grounded-only.
        phi   : log-fluidity adjustment scaling the composite flow law.

    Viscous block (composite, matches inversion_icepack2.py)
        main    : ``A4_base exp(phi)`` dislocation creep, exponent ``n_flow``
                  (a mutable continuation ``Constant``, ramped 1 -> n_flow_val).
        regular : ``alpha`` x linear (n=1) at constant ``H_ref`` -> PD at h=0.

    Friction block (regularized Coulomb, residual closure)
        tau_W = C_w0 exp(theta He) |u|^(1/m_slide)      (Weertman, grounded-gated)
        tau_c = max(c0 N, eps_tauc)                     (Coulomb cap -> 0 afloat)
        tau_b = tau_W tau_c / (tau_W + tau_c)           (Weertman low-u, cap high-u)
        closes ``tau = -tau_b u / |u|_reg`` LINEARLY (identity tau-block).

    Parameters
    ----------
    z, theta, phi : Function
        Mixed state and the two log-controls (live references; the returned form
        tracks in-place updates to ``theta``/``phi`` and to the geometry).
    H, s, b, C_w0 : Function
        Thickness, surface, bed, and the precomputed Weertman anchor.
    A4_base : Constant
        Composite creep prefactor (``A0 * a4_factor``).
    n_flow : Constant
        Mutable viscous exponent for n-continuation (1 -> n_flow_val).
    n_flow_val : float
        Full viscous exponent; fixes the linearization power ``tau_c**(n-1)``.
    m_slide : float
        Sliding exponent -- FIXED (the residual needs no m-continuation).
    tau_c, alpha, H_ref : Constant
        Linearization stress, viscous-regularizer weight, reference thickness.
    c0 : float
        Coulomb-cap coefficient (``tau_c = c0 * N``); ~0.5 (Tsai/Schoof).
    u_min : float
        Velocity regularization [m/yr] -> finite stress at u=0.
    eps_tauc : float
        Yield floor [MPa]; 0 gives bit-exact zero drag where N=0.
    c_w0_floor : float
        Numerical coercivity floor on the Weertman coefficient (same units as
        ``C_w0``).  Needed only when the geometry has ice-free buffer nodes
        (``H = 0``, ``h_clamp = 0``): there ``C_w0 -> 0`` AND the viscous
        coupling ``H*M`` vanishes, so the velocity loses all coercivity and the
        Jacobian goes singular.  Flooring ``C_w0`` restores a vanishingly small
        Weertman drag against the (nonzero) buffer Coulomb cap.  It does NOT leak
        onto real shelves -- those keep ``N = 0 -> tau_cap = 0 -> tau_b = 0``
        exactly -- and it is absorbed by ``theta`` on grounded ice (a
        multiplicative log-control), so the inferred grounded friction is
        unchanged.  0 disables it (correct when ``H`` is clamped > 0).
    h_visc_floor : float
        Thickness floor [m] applied ONLY to the membrane/viscous coupling
        (the ``-h M:eps`` terms), NOT to the driving stress, effective pressure,
        or calving back-pressure.  This is the robust cure for ice-free buffer
        nodes (``H=0``): the membrane term supplies the velocity coercivity the
        vanished ``H*M`` coupling and floored friction cannot, while the TRUE
        ``H=0`` driving stress keeps the buffer velocity ~0 (clamping ``H`` in
        the driving term instead fabricates ``rho g H_floor grad s`` and blows
        the buffer velocity up).  Real ice (``H >> floor``) is unaffected.  0
        disables it (correct when ``H`` is clamped > 0 upstream).
    ocean_drag : float
        Linear drag [MPa yr/m] applied ONLY where ``H < h_ocean`` (ice-free
        buffer / freshly-calved floor cells), ramping linearly to EXACTLY zero
        at ``H = h_ocean``.  Frictionless floor cells (``C_w0 ~ 0``, ``N ~ 0``)
        have no velocity coercivity, and the thick front against them pumps
        momentum through the ``dS`` jump term into the degenerate side - the
        per-solve velocity-runaway equilibrium behind the Jul 2026 forward
        blow-up (RC and Budd alike; dt-independent, collar-immune).  Ports
        gia-icepack ``ase_model.momentum_F`` ocean_drag, which cured the same
        runaway on ASE.  Real ice (incl. 10 m shelves) is untouched; 0 off.
    u_lim, k_lim : float
        Soft speed limiter (gia ``momentum_F`` u_lim/k_lim): supplemental drag
        ``k_lim * max(0, |u| - u_lim)`` [MPa], EXACTLY zero below ``u_lim``
        [m/yr].  ``u_lim ~ 2e4`` is ~5x the fastest observed Antarctic flow, so
        physical flow never feels it; it bounds the retreat-cliff momentum
        pumping wherever the floor-cell pathology sits.  ``u_lim = 0`` off.
    calving_ids : tuple or None
        Outflow boundary ids for the calving-front back-pressure (None to skip).

    Returns
    -------
    F : ufl.Form
        The dual residual; solve ``F == 0`` for ``z``.
    """
    Z = z.function_space()
    mesh = Z.mesh()
    u, M, tau = split(z)
    v, Mt, sig = split(TestFunction(Z))
    d = 2  # depth-averaged (SSA): horizontal dimension is 2

    He = grounded_mask(H, b, gl_width=gl_width)
    N = effective_pressure(H, s)
    # Membrane/viscous coupling thickness: floored for buffer coercivity ONLY;
    # driving stress, N, and calving keep the true H (see h_visc_floor docstring).
    H_visc = max_value(H, Constant(h_visc_floor))

    # ---- composite viscous (membrane-stress) block ----
    # M-stationarity of [viscous_power(glen, H_visc) + alpha*viscous_power(lin,
    # H_ref) + momentum_balance].  Linear regularizer keeps the constant H_ref so
    # the block is positive-definite at h -> 0; both creep terms use H_visc.
    eps = sym(grad(u))
    M2 = (inner(M, M) - tr(M) ** 2 / (d + 1)) / 2
    Mn = conditional(eq(n_flow, 1), Constant(1.0), M2 ** ((n_flow - 1) / 2))
    A_eff = A4_base * exp(phi)
    dev_MMt = inner(M, Mt) - tr(M) * tr(Mt) / (d + 1)
    lin_reg = H_ref * A_eff * tau_c ** (n_flow_val - 1) * dev_MMt   # linear (n=1) @ H_ref
    F = (
        H_visc * A_eff * Mn * dev_MMt                              # dislocation creep (n_flow)
        + alpha * lin_reg                                         # small global linear regularizer (h->0 PD)
        - H_visc * inner(eps, Mt)                                  # strain-rate coupling
    ) * dx
    # GL-gated coercivity collar (Andrew/gia composite_viscosity_gated): a STRONG
    # linear regularizer gated to the near-flotation / floating band via (1 - He),
    # ~0 on the well-grounded trunk (so it keeps its true Glen rheology).  This is
    # what buys back the velocity damping that EXACT-zero shelf drag removes -- the
    # frictionless shelf otherwise has nothing to pin velocity noise (the RC-forward
    # blow-up).  Purely viscous (membrane), so basal drag stays exactly zero afloat.
    if alpha_gl > 0.0:
        F += alpha_gl * (Constant(1.0) - He) * lin_reg * dx

    # ---- regularized-Coulomb friction (basal-stress) block ----
    # C_w0 floored only for solver coercivity in ice-free buffer nodes (see the
    # c_w0_floor docstring): harmless on real shelves (N=0 -> tau_b=0 exactly)
    # and absorbed by theta on grounded ice.
    u_reg = sqrt(inner(u, u) + Constant(u_min) ** 2)
    C_w0_eff = max_value(C_w0, Constant(c_w0_floor))
    tau_W = C_w0_eff * exp(theta * He) * u_reg ** (1.0 / m_slide)   # Weertman, theta grounded-gated

    if fric_law == "budd":
        # Budd: tau_b = tau_W * N_hat, with N_hat the NORMALIZED effective
        # pressure N_eff/N_ref (=1 at the reference/inversion geometry, so the
        # inverted friction is preserved at t=0 and the effective-pressure
        # feedback is a RELATIVE change as the geometry evolves).  Exact-zero
        # shelf via conditional(N>0, ., 0): N_eff=0 afloat -> tau_b=0 exactly.
        # PISM-delta floor (Bueler & van Pelt 2015, till_effective_fraction_
        # overburden, delta ~ 0.02): on GROUNDED ice N_hat >= delta*P_o/N_ref
        # (P_o = local overburden, so the floor evolves and decays as ice
        # thins), removing the near-flotation frictionless-GL degeneracy.  The
        # Joughin cap (reduceNearGLBeta) bounds N_hat above.
        # Floor the denominator ALWAYS (not just for exact-zero guarding): UFL
        # evaluates both conditional branches, so an unguarded N/Nr = 0/0 on the
        # shelf poisons the Jacobian with NaN even though the conditional selects
        # 0.  1e-6 is tiny vs grounded N (~rho_I g H), so N/Nr = 1 stands on all
        # grounded ice when N_ref=None (inversion); only the (conditional-zeroed)
        # shelf sees the floor.
        Nr = max_value(N if N_ref is None else N_ref, Constant(1e-6))
        p_I = rho_I * g * max_value(H, Constant(1.0))
        N_hat = min_value(N / Nr, Constant(nhat_cap))
        if nhat_floor > 0.0:
            novb = Constant(nhat_floor) * p_I / Nr                  # delta * local overburden (normalized)
            N_hat = min_value(max_value(N_hat, novb), Constant(nhat_cap))
        N_hat = conditional(gt(N, Constant(0.0)), N_hat, Constant(0.0))  # EXACT-zero shelf
        tau_b = tau_W * N_hat
    else:  # regularized_coulomb
        tau_cap = max_value(Constant(c0) * N, Constant(eps_tauc))   # Coulomb cap = c0 N -> 0 afloat
        tau_b = tau_W * tau_cap / max_value(tau_W + tau_cap, Constant(1e-15))  # Weertman low-u, cap high-u; ->0 afloat (NaN-safe)

    # Floor-cell coercivity drags (see docstring): both are EXACTLY zero on
    # real ice / physical speeds, so genuine shelves keep tau_b = 0.
    if ocean_drag > 0.0:
        tau_b = tau_b + (Constant(ocean_drag)
                         * max_value(Constant(0.0), Constant(1.0) - H / Constant(h_ocean))
                         * u_reg)
    if u_lim > 0.0:
        # k_lim may be a live Constant: with k_lim = 0 the term (and its
        # Jacobian) vanish identically - no kink, no effect - so a caller
        # can keep the term structurally present and raise k_lim only for
        # rescue solves at runaway-front geometries (the 1873/2024 walls:
        # frozen-iterate speeds of 1e10 m/yr at floating h ~ 0-60 m front
        # nodes; physical flow never reaches u_lim ~ 2e4).
        k_lim_c = k_lim if isinstance(k_lim, Constant) else Constant(k_lim)
        tau_b = tau_b + k_lim_c * max_value(Constant(0.0), u_reg - Constant(u_lim))

    F += inner(tau + tau_b * u / u_reg, sig) * dx                   # tau = -tau_b u/|u| (identity tau-block)

    # ---- momentum balance (hand-written: H_visc on the membrane term, TRUE H on
    # the driving stress and facet term) + calving-front back-pressure ----
    #
    # Geometry may be CG1 or DG0, and the pair of terms below is correct for
    # both: they represent the DISTRIBUTIONAL gradient of a possibly
    # discontinuous surface as broken cell gradient + facet jump.  With CG1 s
    # the jump vanishes and the cell term carries everything; with DG0 s the
    # cell gradient vanishes and the facet term carries everything.  Do not
    # "fix" the seemingly dead grad(s) under DG0 -- deleting it would break the
    # CG1 path, and the facet term is not an add-on but the whole driving
    # stress.  No exterior-boundary counterpart is needed: the distributional
    # gradient is taken over the open domain, and the terminus traction below
    # supplies the boundary condition.
    nu = FacetNormal(mesh)
    F += (-H_visc * inner(M, sym(grad(v))) + inner(tau - rho_I * g * H * grad(s), v)) * dx
    F += rho_I * g * avg(H) * inner(jump(s, nu), avg(v)) * dS
    if calving_ids:
        F += model.variational.calving_terminus(
            velocity=u, thickness=H, surface=s, outflow_ids=tuple(calving_ids)
        )
    return F
