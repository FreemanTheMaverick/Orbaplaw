import numpy as np
from Orbaplaw import Optimization as opt

def FB_func(U,C0,Ws,W2s):
	C=C0@U
	L=0
	for W,W2 in zip(Ws,W2s):
		L+=np.sum(np.diag(C.T@W2@C)-np.diag(C.T@W@C)**2)
	return L

def FB_jac(U,C0,Ws,W2s):
	C=C0@U
	Gamma=0
	for W,W2 in zip(Ws,W2s):
		Gamma+=C0.T@(2*(W2@C)-4*(W@C)*np.diag(C.T@W@C)) # Euclidean derivative
	return Gamma@U.T-U@Gamma.T # Riemannian derivative

def FB_conv(x,f,g,xlast,flast,glast):
	return np.max(np.abs(g))<5e-2 and abs(f-flast)<1e-4

def oldFosterBoys(C0,Ws,W2s,conv):
	func=lambda x:-FB_func(x,C0,Ws,W2s)
	jac=lambda x:-FB_jac(x,C0,Ws,W2s)
	conv_=FB_conv if conv is None else conv
	return opt.Lehtola(C0,0.1,func,jac,4,conv_)

import Maniverse as mv

def FosterBoys(C, Wao, W2aoSum):
	Wrefs = [ C.T @ W @ C for W in Ws ]
	W2refSum = C.T @ W2Sum @ C
	M = mv.Orthogonal(np.eye(Wxref.shape[0]), True)
	def func(U, order):
		Ws = [ U.T @ Wref @ U for Wref in Wrefs]
		DiagWs = [ np.diag(W) for W in Ws ]
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
				for fourWref, DiagW, eightWrefU, UTWref in zip(fourWrefs, DiagWs, eigeWrefUs, UTWrefs):
					Hv -= fourWref @ v @ DiagW + eightWrefU @ np.diag(np.diag(UTWref @ v))
				return Hv
		return L, Ge, He
	L = 0
	tr_setting = mv.TrustRegionSetting()
	mv.TrustRegion(
			func, tr_setting, (1.e-6, 1.e-4, 1.e-7),
			1, 100, L, M, True
	)
	return M.P
