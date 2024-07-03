import numpy as np
import scipy.linalg as sl
from . import FB_func, FB_jac
from . import Fock_func, Fock_jac
from Orbaplaw import Optimization as opt


def Orbitalet_conv(x,f,g,xlast,flast,glast):
    return np.max(np.abs(g))<5e-2 and abs(f-flast)<1e-4

def Orbitalet(C0,Ws,W2s,S,F,gamma,CC,conv):
    epsilon,Cmo=sl.eigh(F,b=S)
    a=S@(Cmo*epsilon)
    F2=a@a.T
    func=lambda x:-(1.-gamma)*FB_func(x,C0,Ws,W2s)-gamma*CC*Fock_func(x,C0,F,F2)
    jac=lambda x:-(1.-gamma)*FB_jac(x,C0,Ws,W2s)-gamma*CC*Fock_jac(x,C0,F,F2)
    conv_=Orbitalet_conv if conv is None else conv
    return opt.Lehtola(C0,0,func,jac,4,conv_)
