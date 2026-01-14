import numpy as np
import Maniverse as mv

class FosterBoysObj(mv.Objective):
	def __init__(self, Wrefs, W2refSum):
		super().__init__()
		self.Wrefs = Wrefs.copy()
		self.W2refSum = W2refSum.copy()
		self.twoW2refSum = 2 * self.W2refSum
		self.fourWrefs = [ 4 * Wref for Wref in Wrefs ]
		self.DiagWs = None
		self.eightWrefUs = None
		self.UTWrefs = None

	def Calculate(self, Us, _):
		U = Us[0]
		Ws = [ U.T @ Wref @ U for Wref in self.Wrefs]
		self.DiagWs = [ np.diag(np.diag(W)) for W in Ws ]
		W2Sum = U.T @ self.W2refSum @ U
		self.Value = np.trace(W2Sum)
		for DiagW in self.DiagWs:
			self.Value -= np.linalg.norm(DiagW) ** 2
		G = 2 * self.W2refSum @ U
		for Wref, DiagW in zip(self.Wrefs, self.DiagWs):
			G -= 4 * Wref @ U @ DiagW
		self.Gradient = [ G ]
		self.eightWrefUs = [ 8 * Wref @ U for Wref in self.Wrefs ]
		self.UTWrefs = [ U.T @ Wref for Wref in self.Wrefs ]

	def Hessian(self, Vs):
		V = Vs[0]
		HV = self.twoW2refSum @ V
		for fourWref, DiagW, eightWrefU, UTWref in zip(self.fourWrefs, self.DiagWs, self.eightWrefUs, self.UTWrefs):
			HV -= fourWref @ V @ DiagW + eightWrefU @ np.diag(np.diag(UTWref @ V ))
		return [[ HV ]]

def FosterBoys(Wrefs, W2refSum):
	obj = FosterBoysObj(Wrefs, W2refSum)
	M = mv.Iterate(obj, [mv.Orthogonal(np.eye(Wrefs[0].shape[0]))], True)
	tr_setting = mv.TrustRegion()
	tol0 = 1e-8 * M.getDimension()
	tol1 = 1e-6 * M.getDimension()
	tol2 = 10
	mv.TruncatedNewton(
			M, tr_setting, (tol0, tol1, tol2),
			0.001, 1000, 1
	)
	return M.Point
