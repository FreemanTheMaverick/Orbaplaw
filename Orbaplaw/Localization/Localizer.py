import numpy as np
import copy as cp
from . import PipekMezey
from . import Fock


def Localizer(mo_mwfn,space="occ",method="PipekMezey",method_optn={}):
    loc_mwfn=cp.deepcopy(mo_mwfn)
    method_string=None
    method_function=None
    if method.upper()=="PM" or "PIPEK" in method.upper() or "MEZEY" in method.upper():
        if space=="mix":
            raise RuntimeError("Fractionally occupied orbitals and mixing occupied and virtual orbitals in Pipek-Mezey localization is not supported!")
        method_string="Pipek-Mezey"
        method_function=PipekMezey
        S=mo_mwfn.Overlap_matrix
        basis_indices_by_center=mo_mwfn.getBasisIndexByCenter()
        charge_type=method_optn.get("charge_type","Lowdin")
        conv=method_optn.get("conv",None)
        args=[(S,basis_indices_by_center,charge_type,conv) for i in range(3)]
    elif method.upper()=="FOCK":
        method_string="Fock"
        method_function=Fock
        S=mo_mwfn.Overlap_matrix
        F=[None for i in range(3)]
        for spin in ([0] if mo_mwfn.Wfntype==0 else [1,2]):
            match spin:
                case 0:
                    F[0]=mo_mwfn.Hamiltonian_matrix
                case 1:
                    F[1]=mo_mwfn.Alpha_Hamiltonian_matrix
                case 2:
                    F[2]=mo_mwfn.Beta_Hamiltonian_matrix
        conv=method_optn.get("conv",None)
        args=[(S,F[i],conv) for i in range(3)]
    print(method_string+" localization:")
    if mo_mwfn.Wfntype==0 or mo_mwfn.Wfntype==1:
        for spin in ([0] if mo_mwfn.Wfntype==0 else [1,2]):
            print("Spin "+str(spin))
            C=mo_mwfn.getCoefficientMatrix(spin)
            arg=args[spin]
            nocc=mo_mwfn.Naelec if spin==1 else mo_mwfn.Nbelec
            if space=="mix":
                print("Localizing occupied orbitals")
                C=C@method_function(C,*arg)
            else:
                if space=="occ" or space=="both":
                    print("Localizing occupied orbitals")
                    C[:,:nocc]=C[:,:nocc]@method_function(C[:,:nocc],*arg)
                if space=="vir" or space=="both":
                    print("Localizing virtual orbitals")
                    C[:,nocc:]=C[:,nocc:]@method_function(C[:,nocc:],*arg)
            loc_mwfn.setCoefficientMatrix(spin,C)
            loc_mwfn.setEnergy(spin,[0 for i in range(mo_mwfn.getNumIndBasis())])
    return loc_mwfn
