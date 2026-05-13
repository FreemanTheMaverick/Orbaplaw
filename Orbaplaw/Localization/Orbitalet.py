import numpy as np
import Maniverse as mv
from . import FosterBoysObj
from . import FockObj

class OrbitaletObj(mv.Objective):
	def __init__(self, Wrefs, W2refSum, Eref, gamma_e):
		super().__init__()
		self.FosterBoysObj = FosterBoysObj(Wrefs, W2refSum)
		self.FockObj = FockObj(Eref)
		self.Wfb = 1. - gamma_e
		self.Wf = gamma_e * 1000

	def Calculate(self, Us, _):
		self.FosterBoysObj.Calculate(Us, _)
		self.FockObj.Calculate(Us, _)
		self.Value = self.Wfb * self.FosterBoysObj.Value + self.Wf * self.FockObj.Value
		self.Gradient = [ self.Wfb * self.FosterBoysObj.Gradient[0] + self.Wf * self.FockObj.Gradient[0] ]

	def Hessian(self, Vs):
		return [ self.Wfb * self.FosterBoysObj.Hessian(Vs)[0] + self.Wf * self.FockObj.Hessian(Vs)[0] ]

def Orbitalet(Wrefs, W2refSum, Eref, gamma_e):
	obj = OrbitaletObj(Wrefs, W2refSum, Eref, gamma_e)
	M = mv.Iterate(obj, [mv.Orthogonal(np.eye(len(Eref)))], True)
	tr = mv.TrustRegion()
	cg = mv.ConjugateGradient(M, 0, 1, (1e-3, 1e-3), M.getDimension(), 1)
	tol0 = 1e-8 * M.getDimension()
	tol1 = 1e-5 * M.getDimension()
	tol2 = 10
	mv.Newton(M, tr, cg, (tol0, tol1, tol2), 1000, 1)
	return M.Point
