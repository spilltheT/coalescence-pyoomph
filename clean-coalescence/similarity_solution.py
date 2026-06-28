#!/usr/bin/env python3
r"""Self-similar master curve for clean viscous sessile-drop coalescence.

In the lubrication (thin-film) limit, two clean sessile drops merging on a
substrate grow a bridge whose neck height obeys the viscous-coalescence law
of Hernandez-Sanchez et al. (PRL 2012),

    h_0(t) = c theta^4 t,   c = 0.272 ,

and whose whole profile collapses onto a single self-similar master shape when
written in the measured-neck-height variables

    H(xi) = h / h_0(t),     xi = theta x / h_0(t) .

Substituting  h = h_0(t) F(eta),  eta = x / w(t),  with h_0 = c theta^4 t and
w = h_0/theta = c theta^3 t  into the clean lubrication equation

    h_t = -(h^3 h_xxx / 3)_x          (sigma = 1, no slip, no Marangoni)

makes the explicit time dependence cancel and yields the autonomous similarity
ODE (identical to the clean limit of Talukdar et al.)

    (F^3 F''')' = 3 c (eta F' - F) ,

with boundary conditions set by symmetry at the neck and matching to the outer
wedge of slope theta:

    F(0) = 1 ,  F'(0) = 0 ,  F'''(0) = 0      (symmetric minimum)
    F'(inf) = 1 ,  F''(inf) = 0               (unit-slope wedge wings)

Because the simulations independently fix the prefactor c = 0.272, we treat c
as known and solve a one-parameter shooting problem for the neck curvature
F''(0): integrate from eta = 0 and adjust F''(0) until F'(eta_max) = 1. The
solution then satisfies F''(eta_max) ~ 0 automatically (|F''| < 4e-3), which
confirms that 0.272 is indeed the eigenvalue of the boundary-value problem.

Running this file writes ``similarity_master.csv`` (columns xi, H over the full
symmetric range) for use by the collapse plot (plot_collapse.py).

Reference: J. F. Hernandez-Sanchez, L. A. Lubbers, A. Eddi & J. H. Snoeijer,
"Symmetric and asymmetric coalescence of drops on a substrate",
Phys. Rev. Lett. 109, 184502 (2012). arXiv:1207.2635.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

C_STAR = 0.272      # viscous-coalescence prefactor h0 = c theta^4 t
ETA_MAX = 10.0      # integration domain for the shoot (asymptotic wedge)


def _rhs(eta: float, y, c: float):
    F, F1, F2, F3 = y
    F4 = (3.0 * c * (eta * F1 - F) - 3.0 * F**2 * F1 * F3) / F**3
    return [F1, F2, F3, F4]


def _integrate(s: float, c: float, eta_max: float):
    """Integrate the similarity ODE from the neck with F''(0) = s.

    Returns the raw ``solve_ivp`` result. Mis-bracketed shoots legitimately
    pinch the film (F -> 0) and stop early with ``success == False`` while still
    yielding a usable truncated end-slope; the shooting residual relies on that,
    so the success check is enforced on the *accepted* solution in
    :func:`solve_master` rather than here.
    """
    return solve_ivp(_rhs, [0.0, eta_max], [1.0, 0.0, s, 0.0], args=(c,),
                     rtol=1e-11, atol=1e-13, dense_output=True, max_step=5e-3)


def solve_master(c: float = C_STAR, eta_max: float = ETA_MAX):
    """Return (F2_0, sol, residual) for the matched self-similar solution.

    F2_0 is the neck curvature F''(0); ``sol`` is the dense ODE solution;
    ``residual`` is F''(eta_max), which should be ~0 if c is the eigenvalue.
    """
    def slope_defect(s):
        return _integrate(s, c, eta_max).y[1, -1] - 1.0   # F'(eta_max) - 1

    s_lo, s_hi = 0.2, 2.0
    d_lo, d_hi = slope_defect(s_lo), slope_defect(s_hi)
    if d_lo * d_hi > 0.0:
        raise RuntimeError(
            f"shooting bracket [{s_lo}, {s_hi}] for F''(0) does not straddle a "
            f"root of F'(eta_max)-1 (defects {d_lo:+.3g}, {d_hi:+.3g}) for "
            f"c={c:.6g}, eta_max={eta_max:.6g}; widen the bracket or adjust c.")
    s = brentq(slope_defect, s_lo, s_hi, xtol=1e-12)
    sol = _integrate(s, c, eta_max)
    if not sol.success:
        raise RuntimeError(
            f"matched shoot F''(0)={s:.6g} failed to integrate to eta_max="
            f"{eta_max:.6g} (c={c:.6g}): {sol.message}; the dense solution would "
            f"extrapolate past the integrated range.")
    return s, sol, float(sol.y[2, -1])


def master_curve(xi_max: float = 5.0, n: int = 1201, c: float = C_STAR):
    """Self-similar master curve H(xi) on the symmetric range [-xi_max, xi_max].

    H is even in xi, so we solve on [0, xi_max] and mirror.
    """
    _, sol, _ = solve_master(c=c, eta_max=max(ETA_MAX, xi_max + 1.0))
    xi_pos = np.linspace(0.0, xi_max, n // 2 + 1)
    H_pos = sol.sol(xi_pos)[0]
    xi = np.concatenate([-xi_pos[::-1][:-1], xi_pos])
    H = np.concatenate([H_pos[::-1][:-1], H_pos])
    return xi, H


def main() -> None:
    here = Path(__file__).resolve().parent
    s, sol, resid = solve_master()
    print(f"c = {C_STAR}  ->  F''(0) = {s:.6f},  F'(eta_max) = {sol.y[1,-1]:.6f},"
          f"  F''(eta_max) = {resid:.2e}  (should be ~0)")
    xi, H = master_curve(xi_max=6.0, n=1601)
    out = here / "similarity_master.csv"
    np.savetxt(out, np.column_stack([xi, H]), delimiter=",",
               header="xi,H  self-similar master curve (clean coalescence, c=0.272)",
               comments="# ")
    print(f"wrote {out}  ({len(xi)} points, xi in [{xi.min():.1f}, {xi.max():.1f}])")


if __name__ == "__main__":
    main()
