import numpy as np
import scipy.linalg as sl
import itertools as it
import copy as cp


def generateNaturalBondOrbital(basis_indices_by_center,P,maxnfragments,maxnorbitals,threshold):
    natoms=1
    norbitals=0
    nnhos=0
    nho_indices=[]
    Presidual=cp.deepcopy(P)
    H=np.zeros_like(P) # Coefficient of primitive NBO in NAO basis
    NHO_nao=np.zeros_like(P) # Coefficient of primitive NHO in NAO basis
    N=np.zeros(P.shape[0])
    finished=False
    while not finished:
        atomlists=list(it.combinations(range(len(basis_indices_by_center)),natoms))
        for atomlist in atomlists:
            indices=[]
            for atom in atomlist:
                indices.append(basis_indices_by_center[atom])
            indices_all=sum(indices,[])
            Pblock=Presidual[np.ix_(indices_all,indices_all)]
            Nblock,Hblock=sl.eigh(Pblock)
            Nblock=Nblock[::-1]
            Hblock=np.fliplr(Hblock)
            for i in range(Nblock.shape[0]):
                if Nblock[i]<1.-threshold:
                    Nblock[i]=0
                else:
                    N[norbitals]=Nblock[i]
                    H[np.ix_(indices_all,[norbitals])]=np.array([Hblock[:,i]]).T
                    these_nho_indices=[]
                    for indices_atom in indices:
                        NHO_nao[np.ix_(indices_atom,[nnhos])]=H[np.ix_(indices_atom,[norbitals])]
                        these_nho_indices+=[nnhos]
                        nnhos+=1
                    nho_indices.append(these_nho_indices)
                    norbitals+=1
            Presidual[np.ix_(indices_all,indices_all)]-=Hblock@np.diag(Nblock)@Hblock.T
            finished= natoms==maxnfragments or norbitals==maxnorbitals
            if finished:
                break
        finished= natoms==maxnfragments or norbitals==maxnorbitals
        natoms+=1

    # Orthogonalization of natural hybrid orbital
    S_nho_nao=NHO_nao[:,:nnhos].T@NHO_nao[:,:nnhos]
    NHO_O_nao=sl.sqrtm(np.linalg.inv(S_nho_nao))
    NHO_nao[:,:nnhos]=NHO_nao[:,:nnhos]@NHO_O_nao
    P_nho=NHO_nao.T@P@NHO_nao

    # Generation of natural bond orbital
    NBO_nho=np.zeros_like(P)
    NBO_n=np.zeros(P.shape[0])
    nnbos=0
    for these_nho_indices in nho_indices:
        n,nbo=np.linalg.eigh(P_nho[np.ix_(these_nho_indices,these_nho_indices)])
        for i in range(len(n)):
            NBO_n[nnbos]=n[i]
            NBO_nho[np.ix_(these_nho_indices,[nnbos])]=np.array([nbo[:,i]]).T
            nnbos+=1
    NBO_nao=NHO_nao@NBO_nho
    #print(sl.norm(NBO_nao.T@NBO_nao-np.diag(np.diag(NBO_nao.T@NBO_nao))))
    #print(np.diag(NBO_nao.T@NBO_nao))
    return NHO_nao,NBO_n,NBO_nao

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
    nho_mwfn=cp.deepcopy(nao_mwfn)
    nbo_mwfn=cp.deepcopy(nao_mwfn)
    if nao_mwfn.Wfntype==0:
        nbasis=nao_mwfn.getNumIndBasis()
        maxnorbitals=nao_mwfn.Naelec if maxnorbitals==-1 else maxnorbitals
        P=nao_mwfn.Extra_info["NAO_density_matrix"]
        NHO_nao,NBO_n,NBO_nao=generateNaturalBondOrbital(basis_indices_by_fragment,P,maxnfragments,maxnorbitals,threshold)
        nho_mwfn.setEnergy([0 for i in range(nbasis)])
        nho_mwfn.setOccupation([0 for i in range(nbasis)])
        nho_mwfn.setCoefficientMatrix(nho_mwfn.getCoefficientMatrix()@NHO_nao)
        nbo_mwfn.setEnergy([0 for i in range(nbasis)])
        nbo_mwfn.setOccupation(NBO_n)
        nbo_mwfn.setCoefficientMatrix(nbo_mwfn.getCoefficientMatrix()@NBO_nao)
        nbo_mwfn.Comment="Natural bond orbital."
    return nho_mwfn,nbo_mwfn

