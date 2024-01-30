import numpy as np
import scipy.linalg as sl
import itertools as it
import copy as cp
from WaveFunction import MultiWaveFunction as mwfn


def generateNaturalHybridOrbital(basis_indices_by_center,D,S,maxnatoms,maxnorbitals,threshold):
    natoms=1
    norbitals=0
    Dresidual=cp.deepcopy(D)
    H=np.zeros_like(D)
    N=np.zeros(D.shape[0])
    finished=False
    while not finished:
        atomlists=list(it.combinations(range(len(basis_indices_by_center)),natoms))
        for atomlist in atomlists:
            indices=[]
            if atomlist is int:
                atomlist=[atomlist]
            for atom in atomlist:
                indices+=basis_indices_by_center[atom]
            Dblock=Dresidual[np.ix_(indices,indices)]
            Sblock=S[np.ix_(indices,indices)]
            Nblock,Hblock=sl.eigh(Dblock,Sblock)
            Nblock=Nblock[::-1]
            print(Nblock)
            Hblock=np.fliplr(Hblock)
            for i in range(Nblock.shape[0]):
                if Nblock[i]<1.-threshold:
                    Nblock[i]=0
                else:
                    N[norbitals]=Nblock[i]
                    H[np.ix_(indices,[norbitals])]=Hblock[:,i]
                    norbitals+=1
            Dresidual[np.ix_(indices,indices)]-=Hblock*Nblock*Hblock.conj().T
            finished= natoms==maxnatoms or norbitals==maxnorbitals
            if finished:
                break
        finished= natoms==maxnatoms or norbitals==maxnorbitals
        natoms+=1
    Ss=sl.sqrtm(S)
    H=Ss@H
    H=sl.orth(H) # wfwfwefwefwefwefw
    return N,H

def NaturalHybridOrbital(mwfn_obj,maxnatoms=-1,maxnorbitals=-1,threshold=0.0075):
    maxnatoms=mwfn_obj.getNumCenters() if maxnatoms==-1 else maxnatoms
    basis_indices_by_center=[]
    for center in mwfn_obj.Centers:
        basis_indices_by_center.append(mwfn_obj.getBasisIndexByCenter(center))
    S=mwfn_obj.Overlap_matrix
    result_mwfn_obj=cp.deepcopy(mwfn_obj)
    if mwfn_obj.Wfntype==0:
        maxnorbitals=mwfn_obj.Naelec if maxnorbitals==-1 else maxnorbitals
        D=mwfn_obj.Total_density_matrix/2.
        N,H=generateNaturalHybridOrbital(basis_indices_by_center,D,S,maxnatoms,maxnorbitals,threshold)
        for iorbital in range(len(N)):
            result_mwfn_obj.Orbitals[iorbital].Energy=0.
            result_mwfn_obj.Orbitals[iorbital].Occ=N[iorbital]
            result_mwfn_obj.Orbitals[iorbital].Coeff=H[:,iorbital]
        result_mwfn_obj.Total_density_matrix=2*H@H.conj().T
        result_mwfn_obj.Hamiltonian_matrix=np.zeros_like(S)
        result_mwfn_obj.Kinetic_energy_matrix=np.zeros_like(S)
        result_mwfn_obj.Potential_energy_matrix=np.zeros_like(S)
    return result_mwfn_obj

