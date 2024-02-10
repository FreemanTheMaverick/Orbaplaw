import numpy as np
import scipy.linalg as sl
import itertools as it
import copy as cp
from WaveFunction import MultiWaveFunction as mwfn


def generateNaturalBondOrbital(basis_indices_by_center,P,maxnatoms,maxnorbitals,threshold):
    natoms=1
    norbitals=0
    Presidual=cp.deepcopy(P)
    H=np.zeros_like(P)
    N=np.zeros(P.shape[0])
    finished=False
    while not finished:
        atomlists=list(it.combinations(range(len(basis_indices_by_center)),natoms))
        for atomlist in atomlists:
            indices=[]
            if atomlist is int:
                atomlist=[atomlist]
            for atom in atomlist:
                indices+=basis_indices_by_center[atom]
            Pblock=Presidual[np.ix_(indices,indices)]
            Nblock,Hblock=sl.eigh(Pblock)
            Nblock=Nblock[::-1]
            Hblock=np.fliplr(Hblock)
            for i in range(Nblock.shape[0]):
                if Nblock[i]<1.-threshold:
                    Nblock[i]=0
                else:
                    #print(atomlist,Nblock[i])
                    N[norbitals]=Nblock[i]
                    H[np.ix_(indices,[norbitals])]=np.array([Hblock[:,i]]).T
                    norbitals+=1
            Presidual[np.ix_(indices,indices)]-=Hblock@np.diag(Nblock)@Hblock.T
            finished= natoms==maxnatoms or norbitals==maxnorbitals
            if finished:
                break
        finished= natoms==maxnatoms or norbitals==maxnorbitals
        natoms+=1
        '''
    HH=sl.orth(H[:,:norbitals])
    H=np.zeros_like(H)
    H[:,:norbitals]=HH
    '''
    return N,H

def NaturalBondOrbital(nao_mwfn,maxnatoms=-1,maxnorbitals=-1,threshold=0.0075):
    maxnatoms=nao_mwfn.getNumCenters() if maxnatoms==-1 else maxnatoms
    basis_indices_by_center=nao_mwfn.getBasisIndexByCenter()
    #S=np.eye(nao_mwfn.getNumBasis())
    #S=nao_mwfn.Extra_info["NAO_overlap_matrix"]
    #S=nao_mwfn.Overlap_matrix
    nbo_mwfn=cp.deepcopy(nao_mwfn)
    if nao_mwfn.Wfntype==0:
        maxnorbitals=nao_mwfn.Naelec if maxnorbitals==-1 else maxnorbitals
        P=nao_mwfn.Extra_info["NAO_density_matrix"]
        #D=nao_mwfn.Total_density_matrix/2
        N,H=generateNaturalBondOrbital(basis_indices_by_center,P,maxnatoms,maxnorbitals,threshold)
        nbo_mwfn.setEnergy([0 for i in range(len(N))])
        nbo_mwfn.setOccupation(N)
        nbo_mwfn.setCoefficientMatrix(nbo_mwfn.getCoefficientMatrix()@H)
        #nbo_mwfn.setCoefficientMatrix(H)
        nbo_mwfn.Comment="Natural bond orbital."
    return nbo_mwfn

