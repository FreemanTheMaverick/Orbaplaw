import numpy as np
import scipy.linalg as sl
import itertools as it
import copy as cp
from Orbaplaw import Localization as loc


def generateNaturalBondOrbital(basis_indices_by_frag,P,C,S,maxnfrags,maxnnbos,threshold):
    # P - The NAO-based density matrix.
    # C - The NAO coefficient matrix in AO basis set.
    # S - The AO-based overlap matrix.

    # Searching for pNBOs and pNHOs.
    pnbo=0 # Number of pNBOs.
    pnho=0 # Number of pNHOs.
    pnbo_info=[] # Information of pNBOs.
    pnho_info=[] # Information of pNHOs.
    Pr=cp.deepcopy(P) # Residual of matrix P.
    H=np.zeros_like(P) # Coefficient matrix of (p-)NBO in NAO basis.
    I=np.zeros_like(P) # Coefficient matrix of (p-)NHO in NAO basis.
    N=np.zeros(P.shape[0]) # Occupation of (p-)NBOs.
    J=np.zeros(P.shape[0]) # Occupation of (p-)NHOs.
    combs=[] # Combinations of fragments.
    for nfrags in range(1,maxnfrags+1): # Generating all possible combinations.
        combs+=list(it.combinations(range(len(basis_indices_by_frag)),nfrags))
    for icomb in range(len(combs)): # Looping over combinations.
        basis_indices_comb=[] # Basis indices of this combination, [fragments] -> [[basis]].
        nfrags=len(combs[icomb]) # Number of fragments.
        for jfrag in range(nfrags): # Looping over fragments of this combination to find all its basis indices.
            frag=combs[icomb][jfrag] # Index of this fragment.
            basis_indices_comb.append(basis_indices_by_frag[frag])
        bic=sum(basis_indices_comb,[]) # Basis indices of this combination, [fragments -> basis].
        Pblock=Pr[np.ix_(bic,bic)] # Block of P matrix belonging to this combination.
        Nblock,Hblock=sl.eigh(Pblock) # The related occupation and orbitals.
        Nblock=Nblock[::-1] # Decreasing order.
        Hblock=np.fliplr(Hblock)
        for jpnbo in range(len(bic)):
            if Nblock[jpnbo]<1.-threshold: # Ignoring all little occupation.
                Nblock[jpnbo]=0
            else:
                N[pnbo]=Nblock[jpnbo] # Copying the occupation number of this combination into the pNBO occupation array.
                H[bic,pnbo]=Hblock[:,jpnbo] # Copying the pNBO into the pNBO coefficient matrix.
                pnho_indices_pnbo=[]
                for kfrag in range(nfrags): # Looping over the fragments of this pNBO to find its pNHOs.
                    basis_indices_frag=basis_indices_comb[kfrag] # Basis indices of this fragment.
                    pnho_indices_pnbo.append(pnho) # Adding the index of this pNHO to the list.
                    pnho_info.append({ # Recording the information about this pNHO.
                        "comb":icomb, # Index of its combination.
                        "pnbo":pnbo, # Index of its pNBO.
                        "frag":combs[icomb][kfrag]}) # Index of its fragment.
                    pnho+=1 # A new pNHO is found.
                pnbo_info.append({ # Recording the information about this NBO.
                    "comb":icomb, # Index of its combination.
                    "pnhos":pnho_indices_pnbo}) # Indices of its pNHOs.
                pnbo+=1 # A new pNBO is found.
        Pr[np.ix_(bic,bic)]-=Hblock@np.diag(Nblock)@Hblock.T # Removing the density of the found pNBOs from the total P.
        if pnbo==maxnnbos: # Exiting the loop when the max number of NBOs are found.
            combs=combs[:icomb+1] # Removing all the unused combinations.
            break

    # Finding degenerate pNBOs and partially localizing them.
    pnbo_indices=[[] for i in range(len(combs))] # pNBO indices of each combination, [combination] -> [[pNBO]].
    for jpnbo in range(pnbo): # Looping over all pNBOs.
        pnbo_indices[pnbo_info[jpnbo]["comb"]].append(jpnbo)
    eigenlists=[] # List of indices of degenerate pNBOs, [degenerate group] -> [[pNBO]].
    for icomb in range(len(combs)):
        if len(combs[icomb])==1: # One-fragment non-bonding orbitals do not need localizing.
            continue
        eigenvalue=114514
        eigenlist=[] # List of indices of degenerate pNBOs in this degenerate group.
        for jpnbo in pnbo_indices[icomb]:
            if not np.isclose(eigenvalue,N[jpnbo]) and eigenlist!=[]: # One group of degenerate pNBOs is found if the eigenvalue of the current pNBO is different from that of last one.
                eigenlists.append(eigenlist) # Recording this degenerate group.
                eigenlist=[] # Clean the degenerate pNBO list.
            eigenvalue=N[jpnbo] # Recording the eigenvalue and the index of the current pNBO.
            eigenlist.append(jpnbo)
            if jpnbo==pnbo_indices[icomb][-1]: # One group of degenerate pNBOs is found if the current pNBO is the last one in this combination.
                eigenlists.append(eigenlist)
    for keigen in range(len(eigenlists)): # Looping over degenerate groups.
        nbid=eigenlists[keigen] # pNBO indices of this degenerate group.
        if len(nbid)==1: # Degenerate groups of only one pNBO, namely non-degenerate-pNBOs, do not need localizing.
            continue
        U=loc.generatePipekMezey(C@H[:,nbid],S,basis_indices_by_frag,loc.LowdinCharge,0) # C - NAO in AO basis; C@H - degenerate pNBO in AO basis
        H[:,nbid]=H[:,nbid]@U # Partially localizing the degenerate pNBOs.

    # Dividing pNBOs into pNHOs.
    for ipnbo in range(pnbo): # Looping over pNBOs.
        for jpnho in pnbo_info[ipnbo]["pnhos"]: # Looping over pNHOs of this pNBO.
            bif=basis_indices_by_frag[pnho_info[jpnho]["frag"]]
            I[bif,jpnho]=H[bif,ipnbo] # Copying the elements of this pNBO belonging to each fragment into the pNHO coefficient matrix.

    # Orthogonalization of pNHOs and generalization of NHOs.
    pnho_indices=[[] for i in range(len(basis_indices_by_frag))] # pNHO indices of each fragment, [fragment] -> [[pNHO]]; len(basis_indices_by_frag) is the total number of fragments.
    for ipnho in range(pnho): # Looping over all pNHOs.
        pnho_indices[pnho_info[ipnho]["frag"]].append(ipnho)
    for ifrag in range(len(pnho_indices)): # Looping over fragments.
        Iblock=I[:,pnho_indices[ifrag]] # The coefficients of the pNHOs of this fragment.
        Sblock=Iblock.T@Iblock # The overlap matrix of pNHOs in NAO basis set. NAOs are mutually orthonormal.
        Jblock=np.diag(np.diag(Iblock.T@P@Iblock)) # Population of pNHOs derived from the fragment density matrix in pNHO basis set.
        Oblock=Jblock@sl.sqrtm(np.linalg.inv(Jblock@Sblock@Jblock)) # Occupancy-weighted symmetric orthogonalization.
        Iblock=Iblock@Oblock
        I[:,pnho_indices[ifrag]]=Iblock
        J[pnho_indices[ifrag]]=np.diag(Iblock.T@P@Iblock) # Occupation numbers of NHOs.

    # Generation of NBO by diagonalization of the NHO-based density matrix.
    nbo=0 # Number of NBOs.
    nbo_info=[] # Information of NBOs.
    N=np.zeros_like(N) # Occupation of NBOs.
    H=np.zeros_like(H) # Temporarily the NBO coefficient matrix is in NHO basis. It will be transformed to NAO basis later.
    for ipnbo in range(pnbo): # Looping over all pNBOs.
        Iblock=I[:,pnbo_info[ipnbo]["pnhos"]] # The NHOs of this pNBO.
        Pblock=Iblock.T@P@Iblock # The density matrix expressed by these NHOs.
        Nblock,Hblock=np.linalg.eigh(Pblock) # NBOs.
        Nblock=Nblock[::-1] # Decreasing order.
        Hblock=np.fliplr(Hblock)
        for jpnho in range(len(pnbo_info[ipnbo]["pnhos"])): # Looping over the component NHOs.
            N[nbo]=Nblock[jpnho]
            H[pnbo_info[ipnbo]["pnhos"],nbo]=Hblock[:,jpnho]
            nbo_info.append({ # Recording the information about this NBO, which is the same as that about its primitive counterpart, except that pNHO is replaced by NHO.
                "comb":pnbo_info[ipnbo]["comb"],
                "nhos":pnbo_info[ipnbo]["pnhos"]})
            nbo+=1 # A new NBO is found.
    '''
    # Checking the orthonormality of NHOs and NBOs
    for A in [I,H]:
        print(sl.norm(I.T@I-np.diag(np.diag(I.T@I))))
        print(np.diag(I.T@I))
    '''

    # Printing NBO and NHO information
    info=""
    for icomb in range(len(combs)):
        info+="Fragment combination "+str(combs[icomb])+"\n"
        for jnbo in range(nbo):
            if nbo_info[jnbo]["comb"]==icomb:
                info+="NBO_"+str(jnbo)+" ("+str(round(N[jnbo],3))+")  ="
                for knho in nbo_info[jnbo]["nhos"]:
                    info+="  "+str(round(H[knho,jnbo],3))+" * NHO_"+str(knho)+" ("+str(round(J[knho],3))+", F_"+str(pnho_info[knho]["frag"])+")"
                info+="\n"
    H=I@H # Transforming the coefficient matrix of NBO from NHO basis to NAO basis
    return J,I,N,H,info


def NaturalBondOrbital(nao_mwfn,frags=[],maxnfrags=-1,maxnnbos=-1,threshold=0.05): # By default, every atom is a fragment, which is the case of NBO. By combining atoms into fragments one extends NBO to natural fragment bond orbital (NFBO).
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
