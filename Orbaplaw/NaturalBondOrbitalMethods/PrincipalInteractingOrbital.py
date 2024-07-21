import numpy as np
import copy as cp


def generatePrincipalInteractingOrbital(basis_indices_by_frag,P):

    # Obtaining NAOs in fragments
    bfA=basis_indices_by_frag[0]
    bfB=basis_indices_by_frag[1]

    # Principal interacting orbital
    Pab=P[np.ix_(bfA,bfB)]
    U,Sigma,VT=np.linalg.svd(Pab)
    T=np.zeros_like(P)
    T[np.ix_(bfA,bfA)]=U # NAO -> PIO
    T[np.ix_(bfB,bfB)]=VT.T
    I=np.zeros(P.shape[0])
    I[:len(Sigma)]=Sigma**2 # PIO-based bond index

    # Ordering PIOs in pairs
    TT=cp.deepcopy(T)
    II=cp.deepcopy(I)
    pair_info=[]
    pio=0
    nopair=set([i for i in range(P.shape[0])])
    for ipair in range(len(bfA)+len(bfB)):
        if ipair<len(bfA) and ipair<len(bfB):
            T[:,pio]=TT[:,bfA[ipair]]
            I[pio]=II[ipair]
            T[:,pio+1]=TT[:,bfB[ipair]]
            I[pio+1]=II[ipair]
            nopair-={bfA[ipair],bfB[ipair]}
            pair_info.append({"pimos":[pio,pio+1],"pios":[pio,pio+1]})
            pio+=2
        elif len(nopair)>0:
            ipio=nopair.pop()
            T[:,pio]=TT[:,ipio]
            I[pio]=II[ipio]
            pair_info.append({"pimos":[pio],"pios":[pio]})
            pio+=1
    Ppio=T.T@P@T # PIO-based density matrix
    N=np.diag(Ppio) # PIO population

    # Principal interacting molecular orbital
    Y=np.zeros_like(T) # PIO -> PIMO
    O=np.zeros_like(N) # Occupation of PIMO
    if len(basis_indices_by_frag)==2:
        for ipair in range(len(pair_info)):
            pips=pair_info[ipair]["pios"]
            Pblock=Ppio[np.ix_(pips,pips)]
            Oblock,Yblock=np.linalg.eigh(Pblock)
            O[pips]=Oblock[::-1]
            Y[np.ix_(pips,pips)]=np.fliplr(Yblock)
    
    return I,N,T,O,Y,pair_info

def PrincipalInteractingOrbital(nao_mwfn,frags):
    basis_indices_by_center=nao_mwfn.getBasisIndexByCenter()
    basis_indices_by_frag=[]
    for frag in frags:
        basis_indices_this_fragment=[]
        for icenter in frag:
            assert icenter<nao_mwfn.getNumCenters(),"Atom index out of range!"
            basis_indices_this_fragment.extend(basis_indices_by_center[icenter])
        basis_indices_by_frag.append(basis_indices_this_fragment)
    pio_mwfn=cp.deepcopy(nao_mwfn)
    pimo_mwfn=cp.deepcopy(nao_mwfn)
    nbasis=nao_mwfn.getNumBasis()
    print("Principal interacting orbitals:")
    if nao_mwfn.Wfntype==0 or nao_mwfn.Wfntype==1:
        for spin in ([0] if nao_mwfn.Wfntype==0 else [1,2]):
            C=nao_mwfn.getCoefficientMatrix(spin)
            P=None
            match spin:
                case 0:
                    P=nao_mwfn.Extra_info["NAO_density_matrix"]
                case 1:
                    P=nao_mwfn.Extra_info["NAO_alpha_density_matrix"]
                case 2:
                    P=nao_mwfn.Extra_info["NAO_beta_density_matrix"]
            I,N,T,O,Y,pair_info=generatePrincipalInteractingOrbital(basis_indices_by_frag,P)
            I*=1 if nao_mwfn.Wfntype==0 else 2
            output="Spin "+str(spin)+"\n"
            if len(basis_indices_by_frag)==2: # Printing PIMO components
                for ipair in range(len(pair_info)):
                    pimos=pair_info[ipair]["pimos"]
                    pios=pair_info[ipair]["pios"]
                    for pimo in pimos:
                        output+="PIMO_"+str(pimo+ (nao_mwfn.getNumIndBasis() if spin==2 else 0))+" ("+str(round(O[pimo],3))+", "+str(round(I[pimo],3))+") ="
                        for pio in pios:
                            output+="  "+str(round(Y[pio,pimo],3))+" * PIO_"+str(pio+ (nao_mwfn.getNumIndBasis() if spin==2 else 0))+" ("+str(round(N[pio],3))+")"
                        output+="\n"
            print(output)
            pio_mwfn.setOccupation(spin,N)
            pio_mwfn.setEnergy(spin,I)
            pio_mwfn.setCoefficientMatrix(spin,C@T)
            pimo_mwfn.setOccupation(spin,O)
            pimo_mwfn.setEnergy(spin,I)
            pimo_mwfn.setCoefficientMatrix(spin,C@T@Y)
    return pio_mwfn,pimo_mwfn
