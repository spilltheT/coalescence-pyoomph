from pyoomph.generic.problem import Problem
from lubrication import * # Import the lubrication equations
# packages for plotting function 
from pyoomph.meshes.meshdatacache import MeshDataEigenModes
from pyoomph.typings import List, Optional, Union
from pyoomph.output.plotting import *

class DropletSpreading(Problem):	
	def __init__(self):
		super(DropletSpreading,self).__init__()
		self.hp=0.0075 # precursor height
		self.sigma=1 # surface tension
		self.R,self.h_center=1,0.5 # initial radius and height of the droplet
		self.theta_eq=pi/8  # equilibrium contact angle
			
					
	def define_problem(self):
		self.set_coordinate_system("axisymmetric")	
		self.add_mesh(LineMesh(N=500,size=5)) # simple line mesh		
		
		h=var("h") # Building disjoining pressure
		disjoining_pressure=5*self.sigma*self.hp**2*self.theta_eq**2*(h**3 - self.hp**3)/(3*h**6)
		
		eqs=LubricationEquations(sigma=self.sigma,disjoining_pressure=disjoining_pressure) # equations
		eqs+=TextFileOutput() # output	
		h_init=maximum(self.h_center*(1-(var("coordinate_x")/self.R)**2),self.hp) # Initial height
		eqs+=InitialCondition(h=h_init) 
		
		self.add_equations(eqs@"domain") # adding the equation

class DropletCoalescence(DropletSpreading):	
	def __init__(self):
		super(DropletCoalescence,self).__init__()
		self.hp=0.0075 # precursor height
		self.sigma=1 # surface tension
		self.R,self.h_center=1,0.5 # initial radius and height of the droplet
		self.theta_eq=pi/8  # equilibrium contact angle
		self.distance=2.5 # droplet distance
		self.Lx=7.5
		self.max_refinement_level=6
					
	def define_problem(self):
		
		# Set coordinate system
		self.set_coordinate_system("axisymmetric")
		
		# Add line mesh to problem
		mesh = LineMesh(N=500,size=5)
		self.add_mesh(mesh)

		h=var("h") # Building disjoining pressure
		disjoining_pressure=5*self.sigma*self.hp**2*self.theta_eq**2*(h**3 - self.hp**3)/(3*h**6)
		
		# Add equations to the problem
		eqs=LubricationEquations(sigma=self.sigma,disjoining_pressure=disjoining_pressure) # equations
		eqs+=TextFileOutput(filetrunk="result") # output
		
		# Height functions of the droplets
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
