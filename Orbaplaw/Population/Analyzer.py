import numpy as np
from . import Lowdin_func


def Analyzer(mo_mwfn,method="Lowdin",method_optn={}):
    natoms=mo_mwfn.getNumCenters()
    S=mo_mwfn.Overlap_matrix
    basis_indices_by_center=mo_mwfn.getBasisIndexByCenter()
    method_string=None
    method_function=None
    if method.upper()=="LOWDIN":
        method_string="Lowdin"
        method_function=Lowdin_func
    print(method_string+" population:")
    if mo_mwfn.Wfntype==0 or mo_mwfn.Wfntype==1:
        for spin in ([0] if mo_mwfn.Wfntype==0 else [1,2]):
            print("Spin "+str(spin))
            nocc=int(mo_mwfn.getNumElec(spin)/(2 if spin==0 else 1))
            C=mo_mwfn.getCoefficientMatrix(spin)[:,:nocc]
            Qs=Lowdin_func(C,S,basis_indices_by_center)
            for iatom in range(natoms):
                Q=Qs[iatom]
                pop=np.sum(np.diag(Q))
                print(iatom,pop*(2 if spin==0 else 1))

