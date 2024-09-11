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
	