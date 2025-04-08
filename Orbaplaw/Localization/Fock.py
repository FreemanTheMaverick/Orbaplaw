import numpy as np
import scipy.linalg as sl
import Maniverse as mv
from Orbaplaw import Optimization as opt

def Fock_func(U,C0,F,F2):
	C=C0@U
	return np.sum(np.diag(C.T@F2@C)-np.diag(C.T@F@C)**2)

def Fock_jac(U,C0,F,F2):
	C=C0@U
	Gamma=C0.T@(2*(F2@C)-4*(F@C)*np.diag(C.T@F@C)) # Euclidean derivative
	return Gamma@U.T-U@Gamma.T # Riemannian derivative

def Fock_conv(x,f,g,xlast,flast,glast):
	return np.max(np.abs(g))<5e-2 and abs(f-flast)<1e-4

def Fock(C0,S,F,conv):
	epsilon,Cmo=sl.eigh(F,b=S)
	a=S@(Cmo*epsilon)
	F2=a@a.T
	func=lambda x:-Fock_func(x,C0,F,F2)
	jac=lambda x:-Fock_jac(x,C0,F,F2)
	conv_=Fock_conv if conv is None else conv
	return opt.Lehtola(C0,0,func,jac,4,conv_)

'''
def Fock(C, E, S):
	Fref = C.T @ S @ C @ np.diag(E) @ C.T @ S @ C
	F2ref = C.T @ S @ C @ np.diag(E)**2 @ C.T @ S @ C
	M = mv.Orthogonal(np.eye(Fref.shape[0]), True)
	def func(U, order):
		F = U.T @ Fref @ U
		F2 = U.T @ F2ref @ U
		L = np.trace(F2) - np.linalg.norm(np.diag(F)) ** 2;
		Fdiag = np.diag(np.diag(F))
		Ge = 2 @ F2ref @ U - 4 * Fref @ U @ Fdiag
		twoF2ref = 2 * F2ref
		fourFref = 4 * Fref
		eightFrefU = 8 * Fref @ U
		UTFref = U.T @ Fref
		def He(v):
			return v
		if order == 2:
			def He(v):
				return twoF2ref @ v
					- fourFref @ v @ Fdiag
					- eightFrefU @ np.diag(np.diag(UTFref @ v))
		return L, Ge, He
	L = 0
	tr_setting = mv.TrustRegionSetting()
	mv.TrustRegion(
			func, tr_setting, (1.e-6, 1.e-4, 1.e-7),
			1, 100, L, M, True
	)
	return M.P
	'''
