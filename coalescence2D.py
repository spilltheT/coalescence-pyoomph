import sys
from pyoomph.materials import * 
import pyoomph.materials.default_materials

from pyoomph.generic.problem import Problem
from pyoomph.meshes.meshdatacache import MeshDataEigenModes
from pyoomph.typings import List, Optional, Union
from lubrication_spreading import * 
from pyoomph.output.plotting import *
from pyoomph.output.meshio import TextFileOutputAlongLine


from pyoomph.expressions.units import *  # units
from pyoomph.expressions.phys_consts import gas_constant  # and the gas constant

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
		self.L=1
		self.theta=10*pi/180
		self.R=self.L/(sin(self.theta))
		self.H=self.R-self.L/(tan(self.theta))
		self.hp=0.0001
		self.Lx=6
		self.max_refinement_level=6

		self.Ca=1e0
		self.Ma=float(sys.argv[1])
		# self.plotter = PlotterTry(self)
			
					
	def define_problem(self):
		self.add_mesh(LineMesh(minimum=-3,size=self.Lx, N=1000)) 
		
		h=var("h") 
		Gamma=var("Gamma")

		# self.sigma=self.sigma-self.ma*Gamma

		eqs=LubricationEquations(Ca=self.Ca,Ma=self.Ma) # equations
		eqs+=MeshFileOutput() # output	
		x=var("coordinate")

		h1=-self.R + self.H + (self.R**2 - (var("coordinate_x") + (2 * self.R * self.H - self.H**2)**(0.5))**2)**(0.5)
		h2=-self.R + self.H + (self.R**2 - (var("coordinate_x") - (2 * self.R * self.H - self.H**2)**(0.5))**2)**(0.5)
		h_init=h_init=maximum(maximum(h1,h2),self.hp) 

		Gamma_init = Gamma_init = -0.5*tanh(1/self.hp * var("coordinate_x")) + 0.5
		
		eqs+=InitialCondition(Gamma=Gamma_init)
		eqs+=InitialCondition(h=h_init) 
		
		eqs+=SpatialErrorEstimator(h=1) # refine based on the height field
		eqs += TextFileOutput()

		# self+=TextFileOutputAlongLine(filename='profile', start=(-3,0), end=(3,0), N =200) @ "domain"
		
		self.add_equations(eqs@"domain") # adding the equation
		
if __name__=="__main__":
	with DropletCoalescence() as problem:
		problem.run(100,outstep=0.1,startstep=0.01,maxstep=50,temporal_error=1,spatial_adapt=1)
