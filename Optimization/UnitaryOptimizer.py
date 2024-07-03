import numpy as np
import scipy.linalg as sl
from . import DirectionSearch as dc
from . import LineSearch as ls


def Lehtola(C0,guess,func,jac,q,conv,func_para=(),jac_para=(),dir_optn={},line_optn={}):
    maxiter=1000
    np.random.seed(0)
    X=np.random.rand(C0.shape[1],C0.shape[1])
    X-=X.T
    X*=guess
    Ulast=np.zeros_like(X)+114514
    U=sl.expm(X)
    Llast=-1919810
    L=-1919810
    Glast=None
    G=None
    Hlast=None
    H=None
    for iiter in range(maxiter):
        if iiter>0:
            Ulast=U.copy()
            Llast=L
            Glast=G.copy()
            Hlast=H.copy()
            if iiter%C0.shape[1]==0:
                Glast=None
                Hlast=None
        L=func(U)
        G=jac(U) # Riemannian derivative
        H=dc.ConjugateGradient(Hlast,G,Glast,dir_optn.get("cgtype","PR"))
        wmax=np.max(np.abs(np.linalg.eigvalsh(1.j*H)))
        x1=np.zeros_like(H)
        x2=H*2*np.pi/wmax/q
        v1=L
        v2=func(sl.expm(x2)@U)
        g1=G
        g2=jac(sl.expm(x2)@U)
        func_for_line_search=lambda x:func(sl.expm(x)@U)
        best_x=ls.MeshCubic(x1,x2,v1,v2,g1,g2,func_for_line_search,(),line_optn.get("num",3))
        U=sl.expm(best_x)@U
        print("Iteration %d:  L = %f; Gmax = %f" % (iiter,L,np.max(np.abs(G))))
        if conv(U,L,G,Ulast,Llast,Glast):
            return U
    print("Convergence not fully achieved!")
    return U
