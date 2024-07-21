import numpy as np
import copy as cp
from . import PipekMezey
from . import FosterBoys
from . import Fock
from . import Orbitalet

calcS='''
if mo_mwfn.Overlap_matrix is None:
    mo_mwfn.calcOverlap()
'''

def Localizer(mo_mwfn,space="occ",method="PipekMezey",method_optn={}):
    loc_mwfn=cp.deepcopy(mo_mwfn)
    method_string=None
    method_function=None
    if method.upper()=="PM" or "PIPEK" in method.upper() or "MEZEY" in method.upper():
        if space=="mix":
            raise RuntimeError("Fractionally occupied orbitals and mixing occupied and virtual orbitals in Pipek-Mezey localization is not supported!")
        method_string="Pipek-Mezey"
        method_function=PipekMezey
        exec(calcS)
        S=mo_mwfn.Overlap_matrix
        basis_indices_by_center=mo_mwfn.getBasisIndexByCenter()
        charge_type=method_optn.get("charge_type","Lowdin")
        conv=method_optn.get("conv",None)
        args=[(S,basis_indices_by_center,charge_type,conv) for i in range(3)]
    elif method.upper()=="FB" or "FOSTER" in method.upper() or "BOYS" in method.upper():
        method_string="Foster-Boys"
        method_function=FosterBoys
        Ws=[-mo_mwfn.X_electric_dipole_moment_matrix,-mo_mwfn.Y_electric_dipole_moment_matrix,-mo_mwfn.Z_electric_dipole_moment_matrix]
        W2s=[-mo_mwfn.XX_electric_quadrupole_moment_matrix,-mo_mwfn.YY_electric_quadrupole_moment_matrix,-mo_mwfn.ZZ_electric_quadrupole_moment_matrix]
        conv=method_optn.get("conv",None)
        args=[(Ws,W2s,conv) for i in range(3)]
    elif method.upper()=="FOCK":
        method_string="Fock"
        method_function=Fock
        exec(calcS)
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
    elif method.upper()=="ORBITALET":
        method_string="Orbitalet"
        method_function=Orbitalet
        Ws=[-mo_mwfn.X_electric_dipole_moment_matrix,-mo_mwfn.Y_electric_dipole_moment_matrix,-mo_mwfn.Z_electric_dipole_moment_matrix]
        W2s=[-mo_mwfn.XX_electric_quadrupole_moment_matrix,-mo_mwfn.YY_electric_quadrupole_moment_matrix,-mo_mwfn.ZZ_electric_quadrupole_moment_matrix]
        conv=method_optn.get("conv",None)
        args=[(Ws,W2s,conv) for i in range(3)]
        exec(calcS)
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
        gamma=method_optn.get("gamma",0.7959)
        CC=method_optn.get("C",1000)
        conv=method_optn.get("conv",None)
        args=[(Ws,W2s,S,F[i],gamma,CC,conv) for i in range(3)]

    print(method_string+" localization:")
    if mo_mwfn.Wfntype==0 or mo_mwfn.Wfntype==1:
        for spin in ([0] if mo_mwfn.Wfntype==0 else [1,2]):
            print("Spin "+str(spin))
            C=mo_mwfn.getCoefficientMatrix(spin)
            arg=args[spin]
            nocc=round(mo_mwfn.getNumElec(spin)/(2 if spin==0 else 1))
            if space=="mix":
                print("Localizing occupied and virtual orbitals together")
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
