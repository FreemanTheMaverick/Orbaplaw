import numpy as np
import scipy.optimize as so
import scipy.linalg as sl
import random as rd
import itertools as it
import copy as cp

def LowdinCharge(C,S,basis_indices_by_center):
    nbasis=C.shape[0]
    norbitals=C.shape[1]
    assert (nbasis,nbasis)==S.shape,"Dimensions of coefficient and overlap matrices do not match!"
    S12C=sl.sqrtm(S)@C
    Qs=[]
    for basis_indices in basis_indices_by_center:
        Qa=np.zeros([norbitals,norbitals])
        for i in range(norbitals):
            for j in range(norbitals):
                for mu in basis_indices:
                    Qa[i,j]+=S12C[mu,i]*S12C[mu,j]
        Qs.append(Qa)
    return Qs

def Lower2Skew(x):
    nrows=int((1+np.sqrt(1+8*len(x)))/2)
    X=np.zeros([nrows,nrows])
    i=0
    for k in range(nrows):
        for l in range(k):
            X[k,l]=x[i]
            X[l,k]=-x[i]
            i+=1
    return X

def Skew2Lower(X):
    nrows=X.shape[0]
    x=[]
    for k in range(nrows):
        for l in range(k):
            x.append(X[k,l])
    return x

def PM_function(x,*args):
    X=Lower2Skew(x)
    C,S,basis_indices_by_center,charge_type,reverse=args
    Qs=charge_type(C@sl.expm(-X),S,basis_indices_by_center)
    function=0
    for Q in Qs:
        for i in range(Q.shape[0]):
            function+=Q[i,i]**2
    return function if reverse else -function

def PM_gradient(x,*args):
    X=Lower2Skew(x)
    C,S,basis_indices_by_center,charge_type,reverse=args
    Qs=charge_type(C@sl.expm(-X),S,basis_indices_by_center)
    gradient=np.array([])
    for k in range(Qs[0].shape[0]):
        for l in range(k):
            element=0
            for Q in Qs:
                element+=4*(Q[k,k]-Q[l,l])*Q[k,l]
            gradient=np.append(gradient,element)
    return gradient if reverse else -gradient

def PM_hessian(x,y,*args):
    X=Lower2Skew(x)
    Y=Lower2Skew(y)
    C,S,basis_indices_by_center,charge_type,reverse=args
    Qs=charge_type(C@sl.expm(-X),S,basis_indices_by_center)
    hessian=np.array([])
    Qks=[]
    kD=Y@Lower2Skew(PM_gradient(x,C,S,basis_indices_by_center,charge_type,reverse))
    for Q in Qs:
        Qks.append(Q@Y)
    for k in range(Qs[0].shape[0]):
        for l in range(k):
            element=0
            for Q,Qk in zip(Qs,Qks):
                element+=4*(Q[l,l]-Q[k,k])*(Qk[k,l]+Qk[l,k])
                element-=8*(Qk[k,k]*Q[k,l]-Qk[l,l]*Q[l,k])
                element+=0.25*(kD[l,k]-kD[k,l])
            hessian=np.append(hessian,element)
    return hessian if reverse else -hessian

def generatePipekMezey(C,S,basis_indices_by_center,charge_type,reverse):
    norbitals=C.shape[1]
    U=np.eye(norbitals)
    max_gamma_error=114514
    while max_gamma_error>1: # Jacobi sweep
        max_gamma_error=0
        for k in range(norbitals):
            for l in range(k):
                Ast=0
                Bst=0
                S12C=sl.sqrtm(S)@C@U
                for basis_indices in basis_indices_by_center:
                    Qkk=0
                    Qll=0
                    Qkl=0
                    for mu in basis_indices:
                        Qll+=S12C[mu,l]*S12C[mu,l]
                        Qkk+=S12C[mu,k]*S12C[mu,k]
                        Qkl+=S12C[mu,k]*S12C[mu,l]
                    Ast+=Qkl**2-0.25*(Qkk-Qll)**2
                    Bst+=Qkl*(Qkk-Qll)
                if np.sqrt(Ast**2+Bst**2)<1e-10:
                    continue
                sin4a=Bst/np.sqrt(Ast**2+Bst**2)
                cos4a=-Ast/np.sqrt(Ast**2+Bst**2)
                four_a=np.sign(sin4a)*np.arccos(cos4a)
                a=four_a/4.
                gamma=a+ (0.25 if reverse else 0) *np.pi
                gamma_error=(gamma-0.25*np.pi)**2 if reverse else abs((gamma-0.5*np.pi)*gamma)
                max_gamma_error=max_gamma_error if max_gamma_error > gamma_error else gamma_error
                U[:,[l,k]]=U[:,[l,k]]@np.array([[np.cos(gamma),np.sin(gamma)],[-np.sin(gamma),np.cos(gamma)]])
    x0=[0 for i in range(int(norbitals*(norbitals-1)/2))]
    PM=so.minimize(PM_function,x0,args=(C@U,S,basis_indices_by_center,charge_type,reverse),method="trust-ncg",jac=PM_gradient,hess="2-point")
    X=Lower2Skew(PM.x)
    return U@sl.expm(-X)

def PipekMezey(mo_mwfn,occ=True,vir=True,charge_type=LowdinCharge):
    S=mo_mwfn.Overlap_matrix
    basis_indices_by_center=mo_mwfn.getBasisIndexByCenter()
    pm_mwfn=cp.deepcopy(mo_mwfn)
    if mo_mwfn.Wfntype==0:
        C=mo_mwfn.getCoefficientMatrix()
        nocc=mo_mwfn.Naelec
        if occ: # Localizing occupied orbitals
            Uocc=generatePipekMezey(C[:,:nocc],S,basis_indices_by_center,charge_type,False)
            C[:,:nocc]=C[:,:nocc]@Uocc
        if vir: # Localizing virtual orbitals
            Uvir=generatePipekMezey(C[:,nocc:],S,basis_indices_by_center,charge_type,False)
            C[:,nocc:]=C[:,nocc:]@Uvir
        pm_mwfn.setCoefficientMatrix(C)
        pm_mwfn.setEnergy([0 for i in range(mo_mwfn.getNumIndBasis())])
        return pm_mwfn
