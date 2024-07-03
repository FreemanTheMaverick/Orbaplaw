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

def FosterBoys(C0,Ws,W2s,conv):
    func=lambda x:-FB_func(x,C0,Ws,W2s)
    jac=lambda x:-FB_jac(x,C0,Ws,W2s)
    conv_=FB_conv if conv is None else conv
    return opt.Lehtola(C0,0.1,func,jac,4,conv_)
