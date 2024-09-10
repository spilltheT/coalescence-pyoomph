from pyoomph import *
from pyoomph.expressions import *

class LubricationEquations(Equations):
	def __init__(self,sigma=1,mu=1,disjoining_pressure=0):
		super(LubricationEquations, self).__init__()
		self.sigma=sigma
		self.mu=mu
		self.disjoining_pressure=disjoining_pressure
		
	def define_fields(self):
		self.define_scalar_field("h","C2")
		self.define_scalar_field("p","C2")		
		
	def define_residuals(self):
		h,eta=var_and_test("h")
		p,q=var_and_test("p")		
		self.add_residual(weak(partial_t(h),eta)+weak(1/self.mu*(h**3/3*grad(p)-h**2/2*grad(self.sigma)),grad(eta)))
		self.add_residual(weak(p-self.disjoining_pressure,q)-weak(self.sigma*grad(h),grad(q)))
		
		
class LubricationProblem(Problem):	
	def define_problem(self):
		self.add_mesh(LineMesh(N=100)) # simple line mesh		
		eqs=LubricationEquations() # equations
		eqs+=TextFileOutput() # output	
		eqs+=InitialCondition(h=0.05*(1+0.25*cos(2*pi*var("coordinate_x"))))  # small height with a modulation
		self.add_equations(eqs@"domain") # adding the equation

		
if __name__=="__main__":
	with LubricationProblem() as problem:
		problem.run(50,outstep=True,startstep=0.25)
