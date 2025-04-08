import numpy as np
import scipy.linalg as sl
import itertools as it
import copy as cp
from Orbaplaw import Localization as loc


class pNBO:
    Combination=None
    pNHOs=None
    Occupancy=None
    Vector=None
    pNMEB=None

class pNHO:
    Fragment=None
    pNBO=None
    Multi=None
    Vector=None
    Siblings=None

class NBO:
    Combination=None
    NHOs=None
    Occupancy=None
    Vector=None
    Multi=None

class NHO:
    Fragment=None
    NBO=None
    Occupancy=None
    Vector=None

def getMatrix(pnbos): # Or pnhos.
    matrix=np.zeros((len(pnbos[0].Vector),len(pnbos)))
    for ipnbo,pnbo in zip(range(len(pnbos)),pnbos):
        matrix[:,ipnbo]=pnbo.Vector
    return matrix

def setMatrix(pnbos,matrix): # Or pnhos.
    for ipnbo,pnbo in zip(range(len(pnbos)),pnbos):
        pnbo.Vector=matrix[:,ipnbo]

def SortpNBO(pnbos,feature): # Or pnhos.
    pnbos_features=dict() # Sorting pNBOs by features, [feature] -> [[pNBO]]. "feature" can be "Combination", "Fragment", etc.
    pnbo_feature=None
    for pnbo in pnbos:
        match feature:
            case "Combination":
                pnbo_feature=pnbo.Combination
            case "Fragment":
                pnbo_feature=pnbo.Fragment
        if pnbo_feature in pnbos_features.keys():
            pnbos_features[pnbo_feature].append(pnbo)
        else:
            pnbos_features[pnbo_feature]=[pnbo]
    return pnbos_features


def generateNaturalBondOrbital(basis_indices_by_frag,P,C,S,maxnfrags,maxnnbos,occ_thres,multi_thres,pdeg_thres,deg_thres):
    # P - The NAO-based density matrix.
    # C - The NAO coefficient matrix in AO basis set.
    # S - The AO-based overlap matrix.

    nbasis=P.shape[0]
    Pr=cp.deepcopy(P) # Residual of matrix P.
    pnbos=[]
    pnhos=[]

    # Searching for pNBOs and pNHOs.
    for nfrags in range(maxnfrags+1):
        combs=list(it.combinations(range(len(basis_indices_by_frag)),nfrags)) # Combinations of fragments.
        for comb in combs:
            bic=[] # Basis indices of this combination, [fragment] -> [[basis]].
            for frag in comb:
                bic.append(basis_indices_by_frag[frag])
            bicf=sum(bic,[]) # Flattened basis indices of this combination, [fragment -> basis]
            Pblock=Pr[np.ix_(bicf,bicf)] # Block of P matrix belonging to this combination
            Nblock,Hblock=np.linalg.eigh(Pblock) # The related occupation and orbitals.
            Nblock=Nblock[::-1] # Putting the eigenvalues and eigenvectors in decreasing order.
            Hblock=np.fliplr(Hblock)
            for ieigen in range(len(bicf)):
                if Nblock[ieigen]<occ_thres: # Ignoring all small occupancies.
                    Nblock[ieigen]=0
                else:
                    pnbo=pNBO()
                    pnbo.Combination=comb
                    pnbo.pNHOs=[]
                    pnbo.Occupancy=Nblock[ieigen]
                    pnbo.Vector=np.zeros(nbasis)
                    pnbo.Vector[bicf]=Hblock[:,ieigen]
                    for frag,bif in zip(comb,bic):
                        pnho=pNHO()
                        pnho.Fragment=frag
                        pnho.pNBO=pnbo
                        pnho.Multi=0
                        pnho.Vector=np.zeros(nbasis)
                        pnho.Siblings=[]
                        pnho.Explored=False
                        pnbo.pNHOs.append(pnho)
                        pnhos.append(pnho)
                    pnbos.append(pnbo)
            Pr[np.ix_(bicf,bicf)]-=Hblock@np.diag(Nblock)@Hblock.T
            if len(pnbos)>=maxnnbos:
                break
        
    # Finding degenerate pNBOs and partially localizing them.
    pnbos_combs=SortpNBO(pnbos,"Combination") # Sorting pNBOs by combinations, [combination : pNBO].
    eigenlists=[] # List of indices of degenerate pNBOs, [degenerate group] -> [[pNBO]].
    for comb,pnbos_comb in pnbos_combs.items():
        if len(comb)==1:
            continue
        eigenvalue=114514
        eigenlist=[] # List of indices of degenerate pNBOs in this degenerate group, [pNBO].
        for pnbo in pnbos_comb:
            if not np.isclose(eigenvalue,pnbo.Occupancy,atol=pdeg_thres) and eigenlist!=[]:
                eigenlists.append(eigenlist) # Recording this degenerate group.
                eigenlist=[] # Clean the degenerate pNBO list.
            eigenvalue=pnbo.Occupancy # Recording the eigenvalue and the index of the current pNBO.
            eigenlist.append(pnbo)
            if pnbo is pnbos_comb[-1]: # One group of degenerate pNBOs is found if the current pNBO is the last one in this combination.
                eigenlists.append(eigenlist)
    for eigenlist in eigenlists: # Looping over degenerate groups.
        if len(eigenlist)==1: # Degenerate groups of only one pNBO, namely non-degenerate-pNBOs, do not need localizing.
            continue
        H=getMatrix(eigenlist)
        U=loc.oldPipekMezey(C@H,S,basis_indices_by_frag,"Lowdin",None) # C - NAO in AO basis; C@H - degenerate pNBO in AO basis
        H=H@U # Partially localizing the degenerate pNBOs.
        setMatrix(eigenlist,H)

    # Dividing pNBOs into pNHOs.
    for pnbo in pnbos: # Looping over pNBOs.
        for pnho in pnbo.pNHOs: # Looping over pNHOs of this pNBO.
            bif=basis_indices_by_frag[pnho.Fragment]
            pnho.Vector=np.zeros_like(pnbo.Vector)
            pnho.Vector[bif]=pnbo.Vector[bif] # Assigning the elements of this pNBO belonging to each fragment to the corresponding pNHO.

    # Identifying multi-electron bonds
    pnhos_frags=SortpNBO(pnhos,"Fragment") # Sorting pNHOs by fragments, [fragment : pNHO].
    for frag,pnhos_frag in pnhos_frags.items():
        for ipnho in range(len(pnhos_frag)):
            pnhoi=pnhos_frag[ipnho]
            if pnhoi.Multi==-1:
                continue
            vectori=pnhoi.Vector/np.linalg.norm(pnhoi.Vector)
            for jpnho in range(ipnho+1,len(pnhos_frag)):
                pnhoj=pnhos_frag[jpnho]
                if pnhoj.Multi==-1:
                    continue
                vectorj=pnhoj.Vector/np.linalg.norm(pnhoj.Vector)
                if abs(np.inner(vectori,vectorj))>multi_thres:
                    pnhoi.Multi=1
                    pnhoi.Siblings.append(pnhoj)
                    pnhoj.Multi=-1
                    pnhoj.Siblings.append(pnhoi)
    for pnhoi in pnhos: # Averaging over sibling pNHOs by population.
        if pnhoi.Multi==1:
            for pnhoj in pnhoi.Siblings:
                pnhoi.Vector+=pnhoj.Vector*np.sign(np.inner(pnhoi.Vector,pnhoj.Vector))
    n_two_e_pnbos=len(pnbos)
    for ipnbo in range(n_two_e_pnbos):
        pnboi=pnbos[ipnbo]
        for jpnbo in range(ipnbo+1,n_two_e_pnbos):
            pnboj=pnbos[jpnbo]
            combined=False
            for pnhok in pnboi.pNHOs:
                for pnhol in pnboj.pNHOs:
                    if pnhok in pnhol.Siblings:
                        combined=True
                        break
            if combined:
                pnmeb=None
                if pnboi.pNMEB!=None and pnboj.pNMEB!=None:
                    continue
                elif pnboi.pNMEB!=None and pnboj.pNMEB==None:
                    pnmeb=pnboi.pNMEB
                    pnboj.pNMEB=pnmeb
                    pnmeb.pNHOs+=[pnho for pnho in pnboj.pNHOs if pnho.Multi!=-1]
                elif pnboi.pNMEB==None and pnboj.pNMEB!=None:
                    pnmeb=pnboj.pNMEB
                    pnboi.pNMEB=pnmeb
                    pnmeb.pNHOs+=[pnho for pnho in pnboi.pNHOs if pnho.Multi!=-1]
                else:
                    pnmeb=pNBO()
                    pnmeb.pNHOs=[pnho for pnho in pnboi.pNHOs+pnboj.pNHOs if pnho.Multi!=-1]
                    pnboi.pNMEB=pnmeb
                    pnboj.pNMEB=pnmeb
                    pnbos.append(pnmeb)

    # Removing redundant pNHOs and pNBOs.
    delete_list=[ipnho for ipnho in range(len(pnhos)) if pnhos[ipnho].Multi==-1]
    '''
    for pnho in pnhos:
        if pnho.Multi==1:
            print(pnhos.index(pnho)+1,pnbos.index(pnho.pNBO)+1,end=' ')
            for pnhoj in pnho.Siblings:
                print(pnhos.index(pnhoj)+1,pnbos.index(pnhoj.pNBO)+1,end=' ')
            print('')
    '''
    delete_list.reverse()
    for ipnho in delete_list:
        pnhos.pop(ipnho)
    delete_list=[ipnbo for ipnbo in range(len(pnbos)) if pnbos[ipnbo].pNMEB!=None]
    delete_list.reverse()
    for ipnbo in delete_list:
        pnbos.pop(ipnbo)
    for pnbo in pnbos:
        if pnbo.Combination==None:
            comb_set=set()
            for pnho in pnbo.pNHOs:
                comb_set.add(pnho.Fragment)
            comb=tuple(comb_set)
            pnbo.Combination=comb

    # Orthogonalization of pNHOs and generation of NHOs.
    pnhos_frags=SortpNBO(pnhos,"Fragment") # Sorting pNHOs by fragments, [fragment : pNHO].
    for frag,pnhos_frag in pnhos_frags.items(): # Looping over fragments.
        Iblock=getMatrix(pnhos_frag) # The coefficients of the pNHOs of this fragment.
        Sblock=Iblock.T@Iblock # The overlap matrix of pNHOs in NAO basis set. NAOs are mutually orthonormal.
        Jblock=np.zeros([len(pnhos_frag),len(pnhos_frag)])
        for ipnho,pnho in zip(range(len(pnhos_frag)),pnhos_frag):
            Jblock[ipnho,ipnho]=np.linalg.norm(pnho.Vector) # Normalization factors of each pNHO.
        Oblock=Jblock@sl.sqrtm(np.linalg.inv(Jblock@Sblock@Jblock)) # Normalization-factor-weighted symmetric orthogonalization.
        Iblock=Iblock@Oblock
        setMatrix(pnhos_frag,Iblock)
        Jblock=np.diag(Iblock.T@P@Iblock)
        for ipnho in range(len(pnhos_frag)):
            pnhos_frag[ipnho].Occupancy=Jblock[ipnho] # Occupation numbers of NHOs.

    # Generation of NBO by diagonalization of the NHO-based density matrix.
    nhos=[]
    nbos=[]
    pnbos_combs=SortpNBO(pnbos,"Combination") # Sorting pNBOs by combinations, [combination : pNBO].
    for comb,pnbos_comb in pnbos_combs.items():
        for pnbo in pnbos_comb:
            nho_index_head=len(nhos)
            nhos_pnho=[]
            for pnho in pnbo.pNHOs:
                nho=NHO()
                nho.Fragment=pnho.Fragment
                nho.Occupancy=pnho.Occupancy
                nho.Vector=pnho.Vector
                nhos.append(nho)
                nhos_pnho.append(nho)
            nho_index_tail=len(nhos)
            nho_indices_pnbo=[i for i in range(nho_index_head,nho_index_tail)]
            Iblock=getMatrix(pnbo.pNHOs) # The NHOs of this pNBO.
            Pblock=Iblock.T@P@Iblock # The density matrix expressed by these NHOs.
            Nblock,Hblock=np.linalg.eigh(Pblock) # NBOs.
            Nblock=Nblock[::-1]
            Hblock=np.fliplr(Hblock)
            eigenlists=[]
            eigenlist=[]
            eigenvalue=114514
            for i in range(len(Nblock)):
                if not np.isclose(Nblock[i],eigenvalue,atol=deg_thres):
                    eigenlists.append(eigenlist)
                    eigenlist=[]
                eigenvalue=Nblock[i]
                eigenlist.append(i)
                if i+1==len(Nblock):
                    eigenlists.append(eigenlist)
            for eigenlist in eigenlists:
                if len(eigenlist)<2:
                    continue
                H=Hblock[:,eigenlist]
                U=loc.PipekMezey(C@Iblock@H,S,basis_indices_by_frag,"Lowdin",None)
                Hblock[:,eigenlist]=H@U
            for i in range(len(Nblock)):
                nbo=NBO()
                nbo.Combination=pnbo.Combination
                nbo.NHOs=nhos_pnho
                nbo.Occupancy=Nblock[i]
                nbo.Vector=np.zeros(nbasis)
                nbo.Vector[nho_indices_pnbo]=Hblock[:,i]
                nbo.Multi=any([pnho.Multi==1 for pnho in pnbo.pNHOs])
                nbos.append(nbo)

    return nhos,nbos


def NaturalBondOrbital(nao_mwfn,frags=[],maxnfrags=-1,maxnnbos=-1,occ_thres=0.95,multi_thres=1,pdeg_thres=1e-5,deg_thres=0): # By default, every atom is a fragment, which is the case of NBO. By combining atoms into fragments one extends NBO to natural fragment bond orbital (NFBO).
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
    if nao_mwfn.Overlap_matrix is None:
        nao_mwfn.calcOverlap()
    S=nao_mwfn.Overlap_matrix
    nbasis=nao_mwfn.getNumBasis()
    nho_mwfn=cp.deepcopy(nao_mwfn)
    nbo_mwfn=cp.deepcopy(nao_mwfn)
    print("Natural (fragment) bond orbitals:")
    if nao_mwfn.Wfntype==0 or nao_mwfn.Wfntype==1:
        for spin in ([0] if nao_mwfn.Wfntype==0 else [1,2]):
            C=nao_mwfn.getCoefficientMatrix(spin)
            if maxnnbos==-1:
                maxnnbos=round(nao_mwfn.getNumElec(1) if spin==1 else nao_mwfn.getNumElec(2))
            P=None
            match spin:
                case 0:
                    P=nao_mwfn.Extra_info["NAO_density_matrix"]
                case 1:
                    P=nao_mwfn.Extra_info["NAO_alpha_density_matrix"]
                case 2:
                    P=nao_mwfn.Extra_info["NAO_beta_density_matrix"]
            if nao_mwfn.Wfntype==0:
                occ_thres*=2
                deg_thres*=2
            nhos,nbos=generateNaturalBondOrbital(basis_indices_by_frag,P,C,S,maxnfrags,maxnnbos,occ_thres,multi_thres,pdeg_thres,deg_thres)
            output="Spin "+str(spin)+"\n" # Printing NBO and NHO information
            nbos_combs=SortpNBO(nbos,"Combination") # Sorting NBOs by combinations, [combination : NBO].
            for comb,nbos_comb in nbos_combs.items():
                output+="Fragment combination "+str(comb)+"\n"
                for nbo in nbos_comb:
                    inbo=nbos.index(nbo)
                    occ=nbo.Occupancy.real
                    nbo_mwfn.Orbitals[inbo+ (nbasis if spin==2 else 0)].Occ=occ
                    output+="NBO_"+str(inbo+ (nbasis if spin==2 else 0))+" ("+str(round(occ,3))+")  ="
                    for nho in nbo.NHOs:
                        jnho=nhos.index(nho)
                        frag=nho.Fragment
                        occ=nho.Occupancy.real
                        nho_mwfn.Orbitals[jnho+ (nbasis if spin==2 else 0)].Occ=occ
                        coef=nbo.Vector[jnho]
                        output+="  "+str(round(coef,3))+" * NHO_"+str(jnho+ (nbasis if spin==2 else 0))+" ("+str(round(occ,3))+", F_"+str(frag)+")"
                    output+="\n"
            print(output)
            I=np.zeros([nbasis,nbasis])
            I[:,:len(nhos)]=getMatrix(nhos)
            H=np.zeros([nbasis,nbasis])
            H[:,:len(nbos)]=getMatrix(nbos)
            nho_mwfn.setEnergy(spin,[0 for i in range(nbasis)])
            nho_mwfn.setCoefficientMatrix(spin,C@I)
            nbo_mwfn.setEnergy(spin,[0 for i in range(nbasis)])
            nbo_mwfn.setCoefficientMatrix(spin,C@I@H)
            nbo_mwfn.GramSchmidt(spin,len(nbos))
    nho_mwfn.Comment="Natural hybrid orbital."
    nbo_mwfn.Comment="Natural bond orbital."
    return nho_mwfn,nbo_mwfn
