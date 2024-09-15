from pyoomph.generic.problem import Problem
from pyoomph.meshes.meshdatacache import MeshDataEigenModes
from pyoomph.typings import List, Optional, Union
from lubrication_spreading import * 
from pyoomph.output.plotting import *
from pyoomph.output.meshio import TextFileOutputAlongLine

class PlotterTry(MatplotlibPlotter):
	def __init__(self, problem: Problem, filetrunk: str = "plot_{:05d}", fileext: str | List[str] = "png", eigenvector: int | None = None, eigenmode: MeshDataEigenModes = "abs", add_eigen_to_mesh_positions: bool = True, position_eigen_scale: float = 1):
		super().__init__(problem, filetrunk, fileext, eigenvector, eigenmode, add_eigen_to_mesh_positions, position_eigen_scale)

	def define_plot(self):
		p = self.get_problem()
		colorbar_1 = self.add_colorbar("h", cmap ='seismic', position = 'top center')

		self.set_view(-3, 0, 0 , 2.05)

		self.add_plot("domain/h", colorbar=colorbar_1)

class DropletCoalescence(DropletSpreading):	
	def __init__(self):
		super(DropletCoalescence,self).__init__()
		self.distance=2.0
		self.hp=0.001
		self.h_center=0.2
		self.Lx=7.5
		self.max_refinement_level=6
		self.plotter = PlotterTry(self)
			
					
	def define_problem(self):
		self.add_mesh(LineMesh(N=500,size=self.Lx)) 
		
		h=var("h") 

		disjoining_pressure=0

		eqs=LubricationEquations(sigma=self.sigma,disjoining_pressure=disjoining_pressure) # equations
		eqs+=MeshFileOutput() # output	
		x=var("coordinate")
		h1=self.h_center*(1-((var("coordinate_x")+1)/self.R)**2) # height functions of the droplets
		h2=self.h_center*(1-((var("coordinate_x")-1)/self.R)**2)		
		h_init=h_init=maximum(maximum(h1,h2),self.hp) # Initial height: maximum of h1, h2 and precursor
		eqs+=InitialCondition(h=h_init) 
		
		eqs+=SpatialErrorEstimator(h=1) # refine based on the height field

		# self+=TextFileOutputAlongLine(filename='profile', start=(-3,0), end=(3,0), N =200) @ "domain"
		
		self.add_equations(eqs@"domain") # adding the equation

		
if __name__=="__main__":
	with DropletCoalescence() as problem:
		problem.run(100,outstep=0.1,startstep=0.01,maxstep=10,temporal_error=1,spatial_adapt=1)
