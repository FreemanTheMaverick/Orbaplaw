import numpy as np
import scipy.linalg as sl
from Orbaplaw import Optimization as opt
from Orbaplaw import Population as pop


def PM_func(U,C0,S,basis_indices_by_center,charge_type):
    C=C0
    if U is not None:
        C=C0@U
    Qs=None
    if charge_type=="Lowdin":
        Qs=pop.Lowdin_func(C,S,basis_indices_by_center)
    L=0
    for Qa in Qs:
        L+=np.sum(np.diag(Qa)**2)
    return L

def PM_jac(U,C0,S,basis_indices_by_center,charge_type):
    C=C0@U
    Qs=None
    QrrUkrs=None
    if charge_type=="Lowdin":
        Qs=pop.Lowdin_func(C,S,basis_indices_by_center)
        QrrUkrs=pop.Lowdin_jac(C,C0,S,basis_indices_by_center)
    Gamma=np.zeros_like(U)
    for Qa,QrrUkra in zip(Qs,QrrUkrs):
        Gamma+=QrrUkra*np.diag(Qa)
    Gamma*=2 # Euclidean derivative
    return Gamma@U.T-U@Gamma.T # Riemannian derivative

def PM_conv(x,f,g,xlast,flast,glast):
    return np.max(np.abs(g))<1e-3 or abs(f-flast)<5e-6

def oldPipekMezey(C0,S,basis_indices_by_center,charge_type,conv):
    func=lambda x:PM_func(x,C0,S,basis_indices_by_center,charge_type)
    jac=lambda x:PM_jac(x,C0,S,basis_indices_by_center,charge_type)
    conv_=PM_conv if conv is None else conv
    return opt.Lehtola(C0,0.001,func,jac,4,conv_)

import Maniverse as mv

class PipekMezeyObj(mv.Objective):
	def __init__(self, Qrefs):
		super().__init__()
		self.Qrefs = Qrefs
		self.Qdiags = None
		self.QrefUs = None

	def Calculate(self, Us, _):
		U = Us[0]
		self.Qdiags = [ np.diag(np.diag(U.T @ Qref @ U)) for Qref in self.Qrefs ]
		self.Value = 0
		for Qdiag in self.Qdiags:
			self.Value -= np.linalg.norm(np.diag(Qdiag)) ** 2
		G = np.zeros_like(U)
		for Qref, Qdiag in zip(self.Qrefs, self.Qdiags):
			G -= Qref @ U @ Qdiag
		self.Gradient = [ 4 * G ]
		self.QrefUs = [ Qref @ U for Qref in self.Qrefs ]

	def Hessian(self, Vs):
		V = Vs[0]
		hess1 = np.zeros_like(V)
		hess2 = np.zeros_like(V)
		for Qref, Qdiag, QrefU in zip(self.Qrefs, self.Qdiags, self.QrefUs):
			hess1 += Qref @ V @ Qdiag
			hess2 += QrefU @ np.diag(np.diag(QrefU.T @ V))
		return [[ - 4 * hess1 - 8 * hess2 ]]

def PipekMezey(Qrefs):
	obj = PipekMezeyObj(Qrefs)
	M = mv.Iterate(obj, [mv.Orthogonal(np.eye(Qrefs[0].shape[1]))], True)
	tr_setting = mv.TrustRegion()
	tol0 = 1e-8 * M.getDimension()
	tol1 = 1e-6 * M.getDimension()
	tol2 = 10
	mv.TruncatedNewton(
			M, tr_setting, (tol0, tol1, tol2),
			0.001, 1000, 1
	)
	return M.Point
