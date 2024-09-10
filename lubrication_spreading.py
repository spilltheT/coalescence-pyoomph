from lubrication import *
		
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

		
if __name__=="__main__":
	with DropletSpreading() as problem:
		problem.run(1000,outstep=True,startstep=0.01,maxstep=10,temporal_error=1)
