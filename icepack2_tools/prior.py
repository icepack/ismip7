"""prior.py - fenics_ice-style Whittle-Matern prior for the ISMIP7 inversion.

Ported from mismip_time-dependent-da/prior.py (Recinos et al. 2023 / fenics_ice
convention). For a control field ``theta`` (= log(param / param_prior), a
log-deviation from a PHYSICAL prior mean) with strength ``gamma`` and
correlation length ``L_reg``, the prior energy added to the inversion cost is

    R(theta; gamma) = 0.5/area * gamma * \\int [ theta^2 + L_reg^2 |grad theta|^2 ] dx

Its Hessian is the prior precision A = delta_eff*M + gamma_eff*K with
delta_eff = gamma/area (mass) and gamma_eff = gamma*L_reg^2/area (stiffness),
so the correlation length sqrt(gamma_eff/delta_eff) = L_reg is fixed and gamma
is the single per-control strength knob.

WHY THIS FIXES THE ISMIP7 n=3 BLOW-UP: the old ISMIP7 regularization was
PURE SMOOTHNESS (gamma*L^2*|grad theta|^2 only, no mass term), so a large
smooth control field was unpenalized -- and with the fluidity control on a
CONSTANT baseline A0, phi had to carry all the spatial fluidity structure and
blew up to +-36 at n=3. The fenics_ice fix is two-fold and does NOT penalize
parameter amplitude: (1) put the control on a PHYSICAL prior mean (thermo
fluidity A_prior, balance friction C_w0), so the amplitude lives in the mean
and theta is a small deviation; (2) use this proper prior whose mass term just
removes the null space, at a physical variance (gamma ~ 1e4 -> prior std
~0.2-0.3, the physical log-deviation scale). The regularization then constrains
the DEVIATION, not the amplitude.
"""

L_REG = 7500.0  # correlation length (m); the one definition used everywhere


def regularization_form(theta, gamma, area, L_reg=L_REG):
    """Whittle-Matern prior energy 0.5/area * gamma * (theta^2 + L_reg^2
    |grad theta|^2) dx, as a UFL form. Add to the inversion cost per control."""
    from firedrake import dx, inner, grad
    coef = 0.5 * float(gamma) / float(area)
    return coef * (theta ** 2
                   + float(L_reg) ** 2 * inner(grad(theta), grad(theta))) * dx


def regularization_gradient_form(theta, test, gamma, area, L_reg=L_REG):
    """Gateaux derivative of regularization_form wrt theta (the assembled
    cost-gradient contribution): gamma/area * (theta*test + L^2 grad.grad)."""
    from firedrake import dx, inner, grad
    coef = float(gamma) / float(area)
    return coef * (theta * test
                   + float(L_reg) ** 2 * inner(grad(theta), grad(test))) * dx
