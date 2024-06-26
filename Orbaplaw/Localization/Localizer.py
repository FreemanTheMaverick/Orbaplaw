import numpy as np
import copy as cp
from . import PipekMezey


def Localizer(mo_mwfn,space="occ",method="PipekMezey",method_optn={}):
    S=mo_mwfn.Overlap_matrix
    basis_indices_by_center=mo_mwfn.getBasisIndexByCenter()
    loc_mwfn=cp.deepcopy(mo_mwfn)
    method_string=None
    method_function=None
    charge_type=None
    if method.upper()=="PM" or "PIPEK" in method.upper() or "MEZEY" in method.upper():
        if space=="mix":
            raise RuntimeError("Fractionally occupied orbitals and mixing occupied and virtual orbitals in Pipek-Mezey localization is not supported!")
        method_string="Pipek-Mezey"
        method_function=PipekMezey
        charge_type=method_optn.get("charge_type","Lowdin")
    print(method_string+" localization:")
    if mo_mwfn.Wfntype==0 or mo_mwfn.Wfntype==1:
        for spin in ([0] if mo_mwfn.Wfntype==0 else [1,2]):
            print("Spin "+str(spin))
            C=mo_mwfn.getCoefficientMatrix(spin)
            nocc=mo_mwfn.Naelec if spin==1 else mo_mwfn.Nbelec
            if space=="mix":
                pass
            else:
                if space=="occ" or space=="both":
                    print("Localizing occupied orbitals")
                    C[:,:nocc]=C[:,:nocc]@method_function(C[:,:nocc],S,basis_indices_by_center,charge_type)
                if space=="vir" or space=="both":
                    print("Localizing virtual orbitals")
                    C[:,nocc:]=C[:,nocc:]@method_function(C[:,nocc:],S,basis_indices_by_center,charge_type)
            loc_mwfn.setCoefficientMatrix(spin,C)
            loc_mwfn.setEnergy(spin,[0 for i in range(mo_mwfn.getNumIndBasis())])
    return loc_mwfn
