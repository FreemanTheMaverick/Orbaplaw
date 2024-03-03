import numpy as np
import scipy.linalg as sl
import itertools as it
import copy as cp
from Orbaplaw import Localization as loc


def generateNaturalBondOrbital(basis_indices_by_frag,P,C,S,maxnfrags,maxnnbos,threshold):
    nfrags=1
    nnbos=0
    nnhos=0
    nho_indices=[] # [combinations] -> [[NBO]] -> [[[NHO]]]
    NHO_INDICES=[] # [combinations] -> [[fragments]] -> [[[NHO]]]
    combinations=[] # List of combinations of fragments
    info=""
    Pr=cp.deepcopy(P) # Residual of matrix P
    H=np.zeros_like(P) # Coefficient matrix of NBO in NAO basis
    I=np.zeros_like(P) # Coefficient matrix of NHO in NAO basis
    N=np.zeros(P.shape[0]) # Occupation of NBO
    finished=False # Whether search is finished
    while not finished:
        fraglists=list(it.combinations(range(len(basis_indices_by_frag)),nfrags)) # Combinations of fragment indices, [combinations] -> [[fragments]] -> [[[basis]]]
        for fraglist in fraglists: # Looping over combinations
            combinations.append(fraglist)
            basis_indices_comb=[] # Basis indices of this combination, [fragments] -> [[basis]]
            for frag in fraglist: # Looping over fragments
                basis_indices_comb.append(basis_indices_by_frag[frag])
            bic=sum(basis_indices_comb,[]) # Basis indices of this combination, [fragments -> basis]
            Pblock=Pr[np.ix_(bic,bic)] # Block of P matrix belonging to this combination
            Nblock,Hblock=sl.eigh(Pblock) # The related occupation and orbitals
            Nblock=Nblock[::-1] # Decreasing order
            Hblock=np.fliplr(Hblock)
            nho_indices_comb=[] # NHO indices of this combination, [NBO] -> [[NHO]]
            for i in range(Nblock.shape[0]):
                if Nblock[i]<1.-threshold: # Ignoring all small occupation
                    Nblock[i]=0
                else:
                    N[nnbos]=Nblock[i] # Copying the occupation number of this combination into the holistic occupation array
                    H[bic,nnbos]=Hblock[:,i] # Copying the NBO into the holistic NBO coefficient matrix
                    nho_indices_nbo=[] # NHO indices of this fragment, [NHO]
                    for basis_indices_frag in basis_indices_comb: # Basis indices of this fragment, [basis]
                        I[basis_indices_frag,nnhos]=H[basis_indices_frag,nnbos] # Copying the elements of NBO belonging to each fragment NHO into the holistic NHO coefficient matrix
                        nho_indices_nbo.append(nnhos) # Adding the index of this NHO to the list
                        nnhos+=1 # One more NHO is found.
                    nho_indices_comb.append(nho_indices_nbo) # Adding the list to the list of lists
                    nnbos+=1 # One more NBO is found.
            nho_indices.append(nho_indices_comb) # Adding the list of lists to the list of lists of lists
            Pr[np.ix_(bic,bic)]-=Hblock@np.diag(Nblock)@Hblock.T # Removing the density of the found NBO from the total P
            finished= nfrags==maxnfrags or nnbos==maxnnbos # Exiting the loop when the max number of interacting fragments or NBOs are found
            if finished:
                break
        finished= nfrags==maxnfrags or nnbos==maxnnbos
        nfrags+=1 # One more interacting fragment is considered.
    #return I,N,H,""

    # Orthogonalization of natural hybrid orbital
    Snho=I[:,:nnhos].T@I[:,:nnhos] # Overlap matrix of pre-orthogonalized NHO
    Onho=sl.sqrtm(np.linalg.inv(Snho)) # Transformation matrix for symmetric orthogonalization

    I[:,:nnhos]=I[:,:nnhos]@Onho.real # Transforming the coefficient matrix of NHO to that of orthogonal NHO

    # Localization of natural hybrid orbital
    for nho_indices_comb in nho_indices: # NHO indices of this combination, [NBO] -> [[NHO]]
        if nho_indices_comb!=[]:
            if len(nho_indices_comb[0])>1:
                NHO_INDICES_COMB=[[] for i in range(len(nho_indices_comb[0]))] # NHO indices of this combination, [fragments] -> [[NHO]]; Capital letters correspond to NHO_INDICES, [combinations] -> [[fragments]] -> [[[NHO]]]
                for nho_indices_nbo in nho_indices_comb:
                    for ifrag in range(len(nho_indices_nbo)):
                        NHO_INDICES_COMB[ifrag].append(nho_indices_nbo[ifrag])
                NHO_INDICES.append(NHO_INDICES_COMB)
                for NHO_INDICES_FRAG in NHO_INDICES_COMB:
                    U=loc.generatePipekMezey(C@I[:,NHO_INDICES_FRAG],S,basis_indices_by_frag,loc.LowdinCharge,0) # C - NAO in AO basis; C@I - degenerate NHO in AO basis
                    I[:,NHO_INDICES_FRAG]=I[:,NHO_INDICES_FRAG]@U
    Pnho=I.T@P@I # Transforming the density matrix of NHO to that of orthogonal and localized NHO
    J=np.diag(Pnho) # Population of NHO
    #return I,N,H,""

    # Generation of natural bond orbital
    H=np.zeros_like(P) # Reseting NBO data, coefficient matrix of NBO in NHO basis
    N=np.zeros(P.shape[0]) # Occupation of NBO
    nnbos=0
    icomb=0
    for nho_indices_comb in nho_indices: # NHO indices of this combination, [NBO] -> [[NHO]]
        info+="Fragment combination "+str(combinations[icomb])+":\n"
        icomb+=1
        nic=sum(nho_indices_comb,[]) # NHO indices of this combination, [NBO -> NHO]
        Nblock,Hblock=np.linalg.eigh(Pnho[np.ix_(nic,nic)])
        degenerate=[114514,[]]
        for i in range(len(Nblock)+1):
            if ( (not np.isclose(degenerate[0],Nblock[i])) if i<len(Nblock) else i==len(Nblock) ) and degenerate[1]!=[]: # Finding degenerate NBOs.
                G=np.zeros([P.shape[0],len(degenerate[1])]) # Degenerate eigenvectors
                for j,k in zip(degenerate[1],range(len(degenerate[1]))):
                    G[nic,k]=Hblock[:,j]
                if len(nho_indices_comb[0])>1: # Partial localization of degenerate NBOs.
                    U=loc.generatePipekMezey(C@I@G,S,basis_indices_by_frag,loc.LowdinCharge,0) # C - NAO in AO basis; C@I@G - degenerate NBO in AO basis
                    G=G@U
                for j,k in zip(degenerate[1],range(len(degenerate[1]))):
                    Hblock[:,j]=G[nic,k]
                degenerate[1]=[]
            if i<len(Nblock):
                degenerate[0]=Nblock[i]
                degenerate[1].append(i)
        Nblock=Nblock[::-1] # Decreasing order
        Hblock=np.fliplr(Hblock)
        for i in range(len(Nblock)):
            N[nnbos]=Nblock[i] # Copying the occupation number of this combination into the holistic occupation array
            H[nic,nnbos]=Hblock[:,i] # Copying the NBO into the holistic NBO coefficient matrix
            info+="NBO_"+str(nnbos)+" ="
            H_nic=H[nic,nnbos]
            square=H_nic**2 # Putting the NBO_NHO coefficients in decreasing order
            order=np.argsort(square)[::-1]
            nic_=np.array(nic)[np.ix_(order)]
            H_nic=H_nic[np.ix_(order)]
            square_sum=0
            for j in nic_:
                if square_sum<0.99:
                    info+="  "+str(round(H[j,nnbos],3))+" * NHO_"+str(j)
                square_sum+=H[j,nnbos]**2
            info+="\n"
            nnbos+=1
 
    H=I@H # Transforming the coefficient matrix of NBO from NHO basis to NAO basis
    #print(sl.norm(NBO_nao.T@NBO_nao-np.diag(np.diag(NBO_nao.T@NBO_nao))))
    #print(np.diag(NBO_nao.T@NBO_nao))

    return J,I,N,H,info


def NaturalBondOrbital(nao_mwfn,frags=[],maxnfrags=-1,maxnnbos=-1,threshold=0.0075): # By default, every atom is a fragment, which is the case of NBO. By combining atoms into fragments one extends NBO to natural fragment bond orbital (NFBO).
    maxnfrags=nao_mwfn.getNumCenters() if maxnfrags==-1 else maxnfrags
    frags=[[i] for i in range(nao_mwfn.getNumCenters())] if frags==[] else frags
    basis_indices_by_center=nao_mwfn.getBasisIndexByCenter()
    basis_indices_by_frag=[]
    for frag in frags:
        basis_indices_this_fragment=[]
        for icenter in frag:
            assert icenter<nao_mwfn.getNumCenters(),"Atom index out of range!"
            basis_indices_this_fragment.extend(basis_indices_by_center[icenter])
        basis_indices_by_frag.append(basis_indices_this_fragment)
    nho_mwfn=cp.deepcopy(nao_mwfn)
    nbo_mwfn=cp.deepcopy(nao_mwfn)
    if nao_mwfn.Wfntype==0:
        nbasis=nao_mwfn.getNumIndBasis()
        C=nao_mwfn.getCoefficientMatrix()
        S=nao_mwfn.Overlap_matrix
        maxnnbos=nao_mwfn.Naelec if maxnnbos==-1 else maxnnbos
        P=nao_mwfn.Extra_info["NAO_density_matrix"]
        NHO_n,NHO_nao,NBO_n,NBO_nao,info=generateNaturalBondOrbital(basis_indices_by_frag,P,C,S,maxnfrags,maxnnbos,threshold)
        print(info)
        nho_mwfn.setEnergy([0 for i in range(nbasis)])
        nho_mwfn.setOccupation(2*NHO_n)
        nho_mwfn.setCoefficientMatrix(nho_mwfn.getCoefficientMatrix()@NHO_nao)
        nho_mwfn.Comment="Natural hybrid orbital."
        nbo_mwfn.setEnergy([0 for i in range(nbasis)])
        nbo_mwfn.setOccupation(2*NBO_n)
        nbo_mwfn.setCoefficientMatrix(nbo_mwfn.getCoefficientMatrix()@NBO_nao)
        nbo_mwfn.Comment="Natural bond orbital."
    return nho_mwfn,nbo_mwfn
