import numpy as np
import Maniverse as mv

class FockObj(mv.Objective):
	def __init__(self, Eref):
		super().__init__()
		self.Fref = np.diag(Eref)
		self.F2ref = np.diag(Eref ** 2)
		self.DiagF = None
		self.twoF2ref = None
		self.fourFref = None
		self.eightFrefU = None
		self.UTFref = None

	def Calculate(self, U, _):
		F = U[0].T @ self.Fref @ U[0]
		F2 = U[0].T @ self.F2ref @ U[0]
		self.DiagF = np.diag(np.diag(F))
		self.Value = np.trace(F2) - np.linalg.norm(self.DiagF) ** 2
		self.Gradient = [ 2 * self.F2ref @ U[0] - 4 * self.Fref @ U[0] @ self.DiagF ]
		self.twoF2ref = 2 * self.F2ref
		self.fourFref = 4 * self.Fref
		self.eightFrefU = 8 * self.Fref @ U[0]
		self.UTFref = U[0].T @ self.Fref

	def Hessian(self, V):
		return [[ self.twoF2ref @ V[0] - self.fourFref @ V[0] @ self.DiagF - self.eightFrefU @ np.diag(np.diag(self.UTFref @ V[0])) ]]

def Fock(Eref):
	obj = FockObj(Eref)
	M = mv.Iterate(obj, [mv.Orthogonal(np.eye(len(Eref)))], True)
	tr_setting = mv.TrustRegion()
	tol0 = 1e-8 * M.getDimension()
	tol1 = 1e-6 * M.getDimension()
	tol2 = 10
	mv.TruncatedNewton(
			M, tr_setting, (tol0, tol1, tol2),
			0.001, 1000, 1
	)
	return M.Point
