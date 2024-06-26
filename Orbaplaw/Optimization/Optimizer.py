import numpy as np
import scipy.linalg as sl
from . import DirectionSearch as dc
from . import LineSearch as ls


def UnitaryOptimizer(func,C0,jac,func_para,jac_para,q,grad_optn={},line_optn={}):
    maxiter=500
    Llast=-114514
    L=None
    Glast=None
    G=None
    Hlast=None
    H=None
    np.random.seed(0)
    X=np.random.rand(C0.shape[1],C0.shape[1])
    X-=X.T
    X*=0.01
    U=sl.expm(X)
    for iiter in range(maxiter):
        if iiter>0:
            Llast=L
            Glast=G.copy()
            Hlast=H.copy()
            if iiter%C0.shape[1]==0:
                Glast=None
                Hlast=None
        L=func(U,*func_para)
        G=jac(U,*jac_para) # Riemannian derivative
        H=dc.ConjugateGradient(Hlast,G,Glast,grad_optn.get("cgtype","PR"))
        wmax=np.max(np.abs(np.linalg.eigvalsh(1.j*H)))
        x1=np.zeros_like(H)
        x2=H*2*np.pi/wmax/q
        v1=L
        v2=func(sl.expm(x2)@U,*func_para)
        g1=G
        g2=jac(sl.expm(x2)@U,*jac_para)
        func_for_line_search=lambda x:func(sl.expm(x)@U,*func_para)
        best_x=ls.MeshCubic(x1,x2,v1,v2,g1,g2,func_for_line_search,(),line_optn.get("num",3))
        U=sl.expm(best_x)@U
        print("Iteration %d:  L = %f; Gmax = %f" % (iiter,L,np.max(np.abs(G))))
        if np.max(np.abs(G))<1e-3 or abs(L-Llast)<5e-6:
            return U
    raise RuntimeError("Convergence failed!")
