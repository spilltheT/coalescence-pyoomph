from pyoomph import *
from pyoomph.expressions import *

class LubricationEquations(Equations):
	def __init__(self,Ca=1,Ma=1):
		super(LubricationEquations, self).__init__()
		self.Ca=Ca
		self.Ma=Ma
		
	def define_fields(self):
		self.define_scalar_field("h","C2")
		self.define_scalar_field("p","C2")		
		self.define_scalar_field("Gamma","C2")
		
	def define_residuals(self):
		h,eta=var_and_test("h")
		p,q=var_and_test("p")
		Gamma,phi=var_and_test("Gamma")
		self.add_residual(weak(partial_t(h),eta) + weak((h**3/3*grad(p) - h**2/2*grad(Gamma)),grad(eta)))
		self.add_residual(weak(self.Ca*p,q)-weak(grad(h) + self.Ca*Gamma*grad(h),grad(q)))
		self.add_residual(weak(partial_t(Gamma),phi) - weak((h*Gamma*grad(Gamma) - h**2/2*Gamma*grad(p)) - 1/self.Ma*grad(Gamma),grad(phi)))
		
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
