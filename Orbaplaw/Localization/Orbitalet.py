import numpy as np
import Maniverse as mv
from . import FB_func
from . import Fock_func

def Orbitalet(Wrefs, W2refSum, Eref, gamma_e):
	M = mv.Orthogonal(np.eye(len(Eref)), True)
	Scale1 = 1. - gamma_e
	Scale2 = gamma_e * 1000
	def func(U, order):
		L1, Ge1, He1 = FB_func(U, order, Wrefs, W2refSum)
		L2, Ge2, He2 = Fock_func(U, order, Eref)
		L = Scale1 * L1 + Scale2 * L2
		Ge = Scale1 * Ge1 + Scale2 * Ge2
		def He(v):
			return Scale1 * He1(v) + Scale2 * He2(v)
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
