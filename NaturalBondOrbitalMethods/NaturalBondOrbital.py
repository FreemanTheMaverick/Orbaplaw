import numpy as np
import scipy.linalg as sl
import itertools as it
import copy as cp


def generateNaturalBondOrbital(basis_indices_by_center,P,maxnfragments,maxnorbitals,threshold):
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
                    N[norbitals]=Nblock[i]
                    H[np.ix_(indices,[norbitals])]=np.array([Hblock[:,i]]).T
                    norbitals+=1
            Presidual[np.ix_(indices,indices)]-=Hblock@np.diag(Nblock)@Hblock.T
            finished= natoms==maxnfragments or norbitals==maxnorbitals
            if finished:
                break
        finished= natoms==maxnfragments or norbitals==maxnorbitals
        natoms+=1
    return N,H

def NaturalBondOrbital(nao_mwfn,fragments=[],maxnfragments=-1,maxnorbitals=-1,threshold=0.0075): # By default, every atom is a fragment, which is the case of NBO. By combining atoms into fragments one extends NBO to natural fragment bond orbital (NFBO).
    maxnfragments=nao_mwfn.getNumCenters() if maxnfragments==-1 else maxnfragments
    fragments=[[i] for i in range(nao_mwfn.getNumCenters())] if fragments==[] else fragments
    basis_indices_by_center=nao_mwfn.getBasisIndexByCenter()
    basis_indices_by_fragment=[]
    for fragment in fragments:
        basis_indices_this_fragment=[]
        for icenter in fragment:
            assert icenter<nao_mwfn.getNumCenters(),"Atom index out of range!"
            basis_indices_this_fragment.extend(basis_indices_by_center[icenter])
        basis_indices_by_fragment.append(basis_indices_this_fragment)
    nbo_mwfn=cp.deepcopy(nao_mwfn)
    if nao_mwfn.Wfntype==0:
        maxnorbitals=nao_mwfn.Naelec if maxnorbitals==-1 else maxnorbitals
        P=nao_mwfn.Extra_info["NAO_density_matrix"]
        N,H=generateNaturalBondOrbital(basis_indices_by_fragment,P,maxnfragments,maxnorbitals,threshold)
        nbo_mwfn.setEnergy([0 for i in range(len(N))])
        nbo_mwfn.setOccupation(N)
        nbo_mwfn.setCoefficientMatrix(nbo_mwfn.getCoefficientMatrix()@H)
        nbo_mwfn.Comment="Natural bond orbital."
    return nbo_mwfn

