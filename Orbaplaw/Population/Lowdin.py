import numpy as np
import scipy.linalg as sl


def Lowdin(C, S, basis_indices_by_center):
	SsqrtC = sl.sqrtm(S) @ C
	Qs = []
	for basis_indices in basis_indices_by_center:
		tmp = SsqrtC[basis_indices, :]
		Qs.append(tmp.T @ tmp)
	return Qs

def Lowdin_func(C,S,basis_indices_by_center):
    nbasis=C.shape[0]
    norbitals=C.shape[1]
    assert (nbasis,nbasis)==S.shape,"Dimensions of coefficient and overlap matrices do not match!"
    S12C=sl.sqrtm(S)@C
    Qs=[]
    for basis_indices in basis_indices_by_center:
        S12Ca=S12C[basis_indices,:]
        Qa=S12Ca.T@S12Ca
        Qs.append(Qa)
    return Qs

def Lowdin_jac(C,C0,S,basis_indices_by_center):
    nbasis=C.shape[0]
    norbitals=C.shape[1]
    assert (nbasis,nbasis)==S.shape,"Dimensions of coefficient and overlap matrices do not match!"
    S12=sl.sqrtm(S)
    S12C0=S12@C0
    S12C=S12@C
    QrrUkrs=[]
    for basis_indices in basis_indices_by_center:
        S12C0a=S12C0[basis_indices,:]
        S12Ca=S12C[basis_indices,:]
        QrrUkrs.append(2*S12C0a.T@S12Ca)
    return QrrUkrs

