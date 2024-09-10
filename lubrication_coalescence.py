from pyoomph.generic.problem import Problem
from pyoomph.meshes.meshdatacache import MeshDataEigenModes
from pyoomph.typings import List, Optional, Union
from lubrication_spreading import * # Import the previous example problem
from pyoomph.output.plotting import *

class PlotterTry(MatplotlibPlotter):
	def __init__(self, problem: Problem, filetrunk: str = "plot_{:05d}", fileext: str | List[str] = "png", eigenvector: int | None = None, eigenmode: MeshDataEigenModes = "abs", add_eigen_to_mesh_positions: bool = True, position_eigen_scale: float = 1):
		super().__init__(problem, filetrunk, fileext, eigenvector, eigenmode, add_eigen_to_mesh_positions, position_eigen_scale)

	def define_plot(self):
		p = self.get_problem()
		xmin=-p.Lx
		xmax = p.Lx
		colorbar_1 = self.add_colorbar("h", cmap ='seismic', position = 'top center')

		self.set_view(-3, 0, 0 , 2.05)

		self.add_plot("domain/h", colorbar=colorbar_1)

class DropletCoalescence(DropletSpreading):	
	def __init__(self):
		super(DropletCoalescence,self).__init__()
		self.distance=2.5 # droplet distance
		self.Lx=7.5
		self.max_refinement_level=6
		self.plotter = PlotterTry(self)
			
					
	def define_problem(self):
		self.add_mesh(RectangularQuadMesh(N=[10,5],size=[self.Lx,self.Lx/2],lower_left=[-self.Lx*0.5,0])) 
		
		h=var("h") # Building disjoining pressure
		disjoining_pressure=5*self.sigma*self.hp**2*self.theta_eq**2*(h**3 - self.hp**3)/(3*h**6)
		
		eqs=LubricationEquations(sigma=self.sigma,disjoining_pressure=disjoining_pressure) # equations
		eqs+=MeshFileOutput() # output	
		x=var("coordinate")
		dist1=x-vector(-self.distance/2,0) # distance to the centers of the droplets
		dist2=x-vector(self.distance/2,0)
		h1=self.h_center*(1-dot(dist1,dist1)/self.R**2) # height functions of the droplets
		h2=self.h_center*(1-dot(dist2,dist2)/self.R**2)		
		h_init=maximum(maximum(h1,h2),self.hp) # Initial height: maximum of h1, h2 and precursor
		eqs+=InitialCondition(h=h_init) 
		
		eqs+=SpatialErrorEstimator(h=1) # refine based on the height field
		
		self.add_equations(eqs@"domain") # adding the equation

		
if __name__=="__main__":
	with DropletCoalescence() as problem:
		problem.run(100,outstep=True,startstep=0.01,maxstep=10,temporal_error=1,spatial_adapt=1)
