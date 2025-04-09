import numpy as np
import Maniverse as mv

def FB_func(U, order, Wrefs, W2refSum):
	Ws = [ U.T @ Wref @ U for Wref in Wrefs]
	DiagWs = [ np.diag(np.diag(W)) for W in Ws ]
	W2Sum = U.T @ W2refSum @ U
	L = np.trace(W2Sum)
	for DiagW in DiagWs:
		L -= np.linalg.norm(DiagW) ** 2
	Ge = 2 * W2refSum @ U
	for Wref, DiagW in zip(Wrefs, DiagWs):
		Ge -= 4 * Wref @ U @ DiagW
	twoW2refSum = 2 * W2refSum
	fourWrefs = [ 4 * Wref for Wref in Wrefs ]
	eightWrefUs = [ 8 * Wref @ U for Wref in Wrefs ]
	UTWrefs = [ U.T @ Wref for Wref in Wrefs ]
	def He(v):
		return v
	if order == 2:
		def He(v):
			Hv = twoW2refSum @ v
			for fourWref, DiagW, eightWrefU, UTWref in zip(fourWrefs, DiagWs, eightWrefUs, UTWrefs):
				Hv -= fourWref @ v @ DiagW + eightWrefU @ np.diag(np.diag(UTWref @ v))
			return Hv
	return L, Ge, He

def FosterBoys(Wrefs, W2refSum):
	M = mv.Orthogonal(np.eye(Wrefs[0].shape[0]), True)
	def func(U, order):
		L, Ge, He = FB_func(U, order, Wrefs, W2refSum)
		return L, Ge, He
	L = 0
	tr_setting = mv.TrustRegionSetting()
	tol0 = 1e-8 * M.getDimension()
	tol1 = 1e-6 * M.getDimension()
	tol2 = 10
	mv.TrustRegion(
			func, tr_setting, (tol0, tol1, tol2),
			0.001, 1, 1000, L, M, 1
	)
	return M.P
