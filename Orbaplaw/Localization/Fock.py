import numpy as np
import scipy.linalg as sl
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
    return opt.Lehtola(C0,func,jac,4,conv_)
