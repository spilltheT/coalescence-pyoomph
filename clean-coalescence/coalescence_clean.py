#!/usr/bin/env python3
"""Clean (surfactant-free) sessile-drop coalescence in the lubrication limit.

This is the beta = Pe = 0 limit of the sessile-drop surfactant model of
Talukdar et al.: two viscous spherical-cap drops sit on a substrate, meet at a
microscopic bridge, and the capillary-driven neck grows self-similarly. Solved
with pyoomph (finite elements + adaptive mesh refinement).

Nondimensionalisation: lengths by the contact radius L, time by L mu / sigma_0.
A precursor film h_inf/L = 1e-4 regularises the contact line.

Run one contact angle, e.g.

    python coalescence_clean.py --theta-deg 10 --t-end 100 --outstep 0.1

Representative contact angles, with the (N, refine, t_end, outstep, maxstep)
used to span the self-similar regime are listed in
README.md. Output (domain/domain_*.txt with columns x, h, p) is written to
``coalescence_clean_theta<deg>/`` and reduced by ``reduce_clean.py``.
"""
import argparse
import os
import sys

from pyoomph import *
from pyoomph.expressions import *
from lubrication_clean import LubricationEquations


class DropletCoalescence(Problem):
    def __init__(self, theta_deg=20.0, N=5000, max_refine=6, hp=1e-4, Lx=6.0):
        super().__init__()
        self.quiet()
        self.L = 1                                  # contact radius (length scale)
        self.theta = theta_deg * pi / 180           # contact angle
        self.R = self.L / sin(self.theta)           # sphere radius
        self.H = self.R - self.L / tan(self.theta)  # apex height
        self.hp = hp                                # precursor film thickness
        self.Lx = Lx                                # domain size [-Lx/2, Lx/2]
        self.N = N                                  # base elements
        self.max_refinement_level = max_refine
        self.sigma = 1                              # surface tension (dimensionless)
        self._step_count = 0
        self._progress_interval = 100

    def actions_after_newton_solve(self):
        self._step_count += 1
        if self._step_count % self._progress_interval == 0:
            sys.stderr.write(f"\rt = {float(self.get_current_time()):.2f}")
            sys.stderr.flush()

    def define_problem(self):
        self.add_mesh(LineMesh(minimum=-self.Lx / 2, size=self.Lx, N=self.N))
        eqs = LubricationEquations(sigma=self.sigma)
        eqs += MeshFileOutput()
        eqs += TextFileOutput()
        # two spherical caps centred at x = +/- sqrt(2 R H - H^2)
        x_center = (2 * self.R * self.H - self.H**2) ** 0.5
        arg1 = maximum(0, self.R**2 - (var("coordinate_x") + x_center) ** 2)
        arg2 = maximum(0, self.R**2 - (var("coordinate_x") - x_center) ** 2)
        h1 = -self.R + self.H + arg1**0.5
        h2 = -self.R + self.H + arg2**0.5
        h_init = maximum(maximum(h1, h2), self.hp)
        eqs += InitialCondition(h=h_init)
        eqs += SpatialErrorEstimator(h=1)
        self.add_equations(eqs @ "domain")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--theta-deg", type=float, default=20.0)
    ap.add_argument("--N", type=int, default=5000)
    ap.add_argument("--max-refine", type=int, default=6)
    ap.add_argument("--t-end", type=float, default=100.0)
    ap.add_argument("--outstep", type=float, default=0.1)
    ap.add_argument("--maxstep", type=float, default=50.0)
    ap.add_argument("--outdir", default=None,
                    help="output directory (default coalescence_clean_theta<deg>)")
    a = ap.parse_args()

    outdir = a.outdir or f"coalescence_clean_theta{int(round(a.theta_deg))}"
    with DropletCoalescence(theta_deg=a.theta_deg, N=a.N,
                            max_refine=a.max_refine) as problem:
        problem.set_output_directory(outdir)
        with open(os.devnull, "w") as devnull:        # silence C-level mesh chatter
            old = os.dup(1)
            os.dup2(devnull.fileno(), 1)
            try:
                problem.run(a.t_end, outstep=a.outstep, startstep=0.001,
                            maxstep=a.maxstep, temporal_error=1, spatial_adapt=1)
            finally:
                os.dup2(old, 1)
                os.close(old)
        print(f"\nSimulation complete -> {outdir}/")


if __name__ == "__main__":
    main()
