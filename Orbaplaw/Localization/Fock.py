import numpy as np
import Maniverse as mv

def Fock_func(U, order, Eref):
	Fref = np.diag(Eref)
	F2ref = np.diag(Eref**2)
	F = U.T @ Fref @ U
	DiagF = np.diag(np.diag(F))
	F2 = U.T @ F2ref @ U
	L = np.trace(F2) - np.linalg.norm(DiagF) ** 2
	Ge = 2 * F2ref @ U - 4 * Fref @ U @ DiagF
	twoF2ref = 2 * F2ref
	fourFref = 4 * Fref
	eightFrefU = 8 * Fref @ U
	UTFref = U.T @ Fref
	def He(v):
		return v
	if order == 2:
		def He(v):
			Hv = twoF2ref @ v - fourFref @ v @ DiagF - eightFrefU @ np.diag(np.diag(UTFref @ v))
			return Hv
	return L, Ge, He

def Fock(Eref):
	M = mv.Orthogonal(np.eye(len(Eref)), True)
	def func(U, order):
		L, Ge, He = Fock_func(U, order, Eref)
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
