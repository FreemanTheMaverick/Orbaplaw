import numpy as np
import scipy.linalg as sl
from Orbaplaw import Population as pop
from Orbaplaw import Optimization as opt


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

def PipekMezey(C0,S,basis_indices_by_center,charge_type,conv):
    func=lambda x:PM_func(x,C0,S,basis_indices_by_center,charge_type)
    jac=lambda x:PM_jac(x,C0,S,basis_indices_by_center,charge_type)
    conv_=PM_conv if conv is None else conv
    return opt.Lehtola(C0,0.001,func,jac,4,conv_)
