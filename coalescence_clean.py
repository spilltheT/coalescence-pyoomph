import os
import sys
from pyoomph import *
from pyoomph.expressions import *
from lubrication_clean import LubricationEquations

class DropletCoalescence(Problem):
	def __init__(self):
		super(DropletCoalescence, self).__init__()
		self.quiet()  # suppress mesh refinement messages
		# Geometry (see paper §2.1)
		self.L = 1                      # contact line radius (length scale)
		self.theta = 20 * pi / 180      # contact angle (20° as in paper)
		self.R = self.L / sin(self.theta)  # sphere radius
		self.H = self.R - self.L / tan(self.theta)  # apex height
		self.hp = 1e-4                  # precursor film thickness h_∞/L = 10⁻⁴
		self.Lx = 6                     # domain size [-3, 3]
		self.N = 1000                   # number of elements
		self.max_refinement_level = 6

		# Physical parameters
		self.sigma = 1  # surface tension (dimensionless)
		self._step_count = 0
		self._progress_interval = 100  # print progress every N timesteps

	def actions_after_newton_solve(self):
		self._step_count += 1
		if self._step_count % self._progress_interval == 0:
			t = float(self.get_current_time())
			sys.stderr.write(f"\rt = {t:.2f}")
			sys.stderr.flush()

	def define_problem(self):
		# 1D mesh from x = -3 to x = 3
		self.add_mesh(LineMesh(minimum=-3, size=self.Lx, N=self.N))

		h = var("h")

		# Clean lubrication equations (no surfactant, no disjoining pressure)
		eqs = LubricationEquations(sigma=self.sigma)
		eqs += MeshFileOutput()
		eqs += TextFileOutput()

		# Spherical cap height profile for two droplets
		# Droplet centers at x = ±sqrt(2RH - H^2)
		x_center = (2 * self.R * self.H - self.H**2)**(0.5)
		# Clamp sqrt arguments to avoid complex values outside droplet footprint
		arg1 = maximum(0, self.R**2 - (var("coordinate_x") + x_center)**2)
		arg2 = maximum(0, self.R**2 - (var("coordinate_x") - x_center)**2)
		h1 = -self.R + self.H + arg1**(0.5)
		h2 = -self.R + self.H + arg2**(0.5)
		h_init = maximum(maximum(h1, h2), self.hp)

		eqs += InitialCondition(h=h_init)
		eqs += SpatialErrorEstimator(h=1)

		self.add_equations(eqs @ "domain")

if __name__ == "__main__":
	with DropletCoalescence() as problem:
		# Redirect C-level stdout to suppress mesh refinement messages
		# Progress updates go to stderr which remains visible
		with open(os.devnull, 'w') as devnull:
			old_stdout_fd = os.dup(1)
			os.dup2(devnull.fileno(), 1)
			try:
				problem.run(100, outstep=0.1, startstep=0.001, maxstep=50,
				            temporal_error=1, spatial_adapt=1)
			finally:
				os.dup2(old_stdout_fd, 1)
				os.close(old_stdout_fd)
		print("\nSimulation complete.")
