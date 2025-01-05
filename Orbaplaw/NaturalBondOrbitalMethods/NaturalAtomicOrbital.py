import numpy as np
import scipy.linalg as sl
import copy as cp


def generateNaturalAtomicOrbital(shell_indices_by_center,basis_indices_by_shell,basis_indices_by_center,angulars,D,S,minimal_shells):
    P=S@D@S

    # 2: Intraatomic orthogonalization

    # 2a: Transformation from Cartesian to pure d, f, g AOs.
    # I do not want to code this. Tell users to use pure AOs only.

    # 2b: Partitioning and symmetry averaging of P and S.
    W=np.zeros(P.shape[0])
    N=np.zeros_like(P)
    for shell_indices_this_center in shell_indices_by_center:
        angulars_this_center=[angulars[shell_index] for shell_index in shell_indices_this_center]
        angulars_this_center_set=set(angulars_this_center)
        for angular in angulars_this_center_set:
            shell_indices=[shell_index for shell_index in shell_indices_this_center if angulars[shell_index]==angular]
            basis_index_heads=[basis_indices_by_shell[shell_index][0] for shell_index in shell_indices]
            PAl=np.zeros([len(basis_index_heads),len(basis_index_heads)])
            SAl=np.zeros([len(basis_index_heads),len(basis_index_heads)])
            for m in range(2*abs(angular)+1):
                basis_indices=[basis_index+m for basis_index in basis_index_heads]
                PAl+=P[np.ix_(basis_indices,basis_indices)]
                SAl+=S[np.ix_(basis_indices,basis_indices)]
            PAl/=2*abs(angular)+1
            SAl/=2*abs(angular)+1

    # 2c: Formation of pre-NAOs.
            OL=np.linalg.inv(sl.sqrtm(SAl))
            PAlL=OL.T@PAl@OL
            w,NL=sl.eigh(PAlL)
            for m in range(2*abs(angular)+1):
                basis_indices=[basis_index+m for basis_index in basis_index_heads]
                W[np.ix_(basis_indices)]=w
                N[np.ix_(basis_indices,basis_indices)]=OL@NL

    # 3: Initial division between orthogonal valence and Rydberg AO spaces.
    
    # 3a: Selection of NMB orbitals.
    WW=cp.deepcopy(W)
    NN=cp.deepcopy(N)
    W=np.array([])
    N=[]
    natoms=len(minimal_shells)
    basis_indices_nmb=[]
    basis_indices_nrb=[]
    jbasis=0
    for iatom in range(natoms):
        W_nmb_atom=np.array([])
        N_nmb_atom=[]
        old_shell_indices_nmb=[]
        old_basis_indices_nmb=[]
        W_nrb_atom=np.array([])
        N_nrb_atom=[]
        old_shell_indices_nrb=[]
        old_basis_indices_nrb=[]
        for l in range(len(minimal_shells[iatom])):
            shell_indices=[]
            ws=[]
            for shell_index in shell_indices_by_center[iatom]:
                if abs(angulars[shell_index])==l:
                    shell_indices.append(shell_index)
                    ws.append(WW[basis_indices_by_shell[shell_index][0]])
            args=np.argsort(ws)[::-1]
            for p in range(len(args)):
                if p<minimal_shells[iatom][l]: # NMB
                    for m in range(2*l+1):
                        W_nmb_atom=np.append(W_nmb_atom,WW[basis_indices_by_shell[shell_indices[args[p]]][0]+m])
                        N_nmb_atom.append(NN[:,basis_indices_by_shell[shell_indices[args[p]]][0]+m])
                    old_shell_indices_nmb.append(shell_indices[args[p]])
                else: # NRB
                    for m in range(2*l+1):
                        W_nrb_atom=np.append(W_nrb_atom,WW[basis_indices_by_shell[shell_indices[args[p]]][0]+m])
                        N_nrb_atom.append(NN[:,basis_indices_by_shell[shell_indices[args[p]]][0]+m])
                    old_shell_indices_nrb.append(shell_indices[args[p]])
        basis_indices_nmb_this_atom=[]
        basis_indices_nrb_this_atom=[]
        for old_shell_index in old_shell_indices_nmb:
            basis_indices_nmb_this_atom_this_shell=[]
            for m in range(2*abs(angulars[old_shell_index])+1):
                basis_indices_nmb_this_atom_this_shell.append(jbasis)
                jbasis+=1
            basis_indices_nmb_this_atom.append(basis_indices_nmb_this_atom_this_shell)
        for old_shell_index in old_shell_indices_nrb:
            basis_indices_nrb_this_atom_this_shell=[]
            for m in range(2*abs(angulars[old_shell_index])+1):
                basis_indices_nrb_this_atom_this_shell.append(jbasis)
                jbasis+=1
            basis_indices_nrb_this_atom.append(basis_indices_nrb_this_atom_this_shell)
        basis_indices_nmb.append(basis_indices_nmb_this_atom)
        basis_indices_nrb.append(basis_indices_nrb_this_atom)
        W_atom=np.append(W_nmb_atom,W_nrb_atom)
        W=np.append(W,W_atom)
        N_nmb_atom=np.array(N_nmb_atom)
        N_nrb_atom=np.array(N_nrb_atom)
        N.extend(N_nmb_atom)
        N.extend(N_nrb_atom)
    N=np.array(N).T

    # 4a: Weighted interatomic orthogonalization within NMB.
    all_basis_indices_nmb=sum(sum(basis_indices_nmb,[]),[]) # Flattening a list of lists of lists.
    all_basis_indices_nrb=sum(sum(basis_indices_nrb,[]),[])
    Spnao=N.T@S@N
    Snmb=Spnao[np.ix_(all_basis_indices_nmb,all_basis_indices_nmb)]
    Wnmb=np.diag(W[np.ix_(all_basis_indices_nmb)])
    Ownmb=Wnmb@np.linalg.inv(sl.sqrtm(Wnmb@Snmb@Wnmb))
    N[:,np.ix_(all_basis_indices_nmb)]=N[:,np.ix_(all_basis_indices_nmb)]@Ownmb
    Spnao=N.T@S@N

    # 3b: Schmidt interatomic orthogonalization of NRB to NMB orbitals.
    for a in range(Spnao.shape[0]):
        for r in all_basis_indices_nrb:
            for m in all_basis_indices_nmb:
                N[a,r]-=N[a,m]*Spnao[m,r]
    Spnao=N.T@S@N

    # 3c: Restoration of natural character of the NRB.
    W.setflags(write=True)
    Ppnao=N.T@P@N
    Nryd=np.zeros_like(N)
    Nryd[np.ix_(all_basis_indices_nmb,all_basis_indices_nmb)]=np.eye(len(all_basis_indices_nmb))
    for basis_indices_this_atom in basis_indices_nrb:
        lmax=int((max(map(len,basis_indices_this_atom))-1)/2) if len(basis_indices_this_atom)>0 else -2 # Skipping this part if a minimal basis is used.
        for l in range(lmax+1):
            basis_index_heads=[]
            for basis_index_this_atom_this_shell in basis_indices_this_atom:
                if 2*l+1==len(basis_index_this_atom_this_shell):
                    basis_index_heads.append(basis_index_this_atom_this_shell[0])
            PAl=np.zeros([len(basis_index_heads),len(basis_index_heads)])
            SAl=np.zeros([len(basis_index_heads),len(basis_index_heads)])
            for m in range(2*l+1):
                basis_indices=[basis_index+m for basis_index in basis_index_heads]
                PAl+=Ppnao[np.ix_(basis_indices,basis_indices)]
                SAl+=Spnao[np.ix_(basis_indices,basis_indices)]
            PAl/=2*l+1
            SAl/=2*l+1
            OL=np.linalg.inv(sl.sqrtm(SAl))
            PAlL=OL.T@PAl@OL
            w,NL=sl.eigh(PAlL)
            for m in range(2*l+1):
                basis_indices=[basis_index+m for basis_index in basis_index_heads]
                W[np.ix_(basis_indices)]=w[::-1]
                Nryd[np.ix_(basis_indices,basis_indices)]=(OL@NL)[:,::-1]
    N=N@Nryd
    Spnao=N.T@S@N


    # 4: Formation of the final NAO set
    # 4a: Weighted interatomic orthogonalization within NRB.
    W=np.diag(W)
    Ow=np.real(W@np.linalg.inv(sl.sqrtm(W@Spnao@W)))
    N=N@Ow
    Spnao=N.T@S@N
    W=np.diag(W)
    '''
    print(sl.norm(Spnao[np.ix_(all_basis_indices_nmb,all_basis_indices_nmb)]-np.eye(len(all_basis_indices_nmb))))
    print(np.diag(Spnao[np.ix_(all_basis_indices_nrb,all_basis_indices_nrb)]))
    print(sl.norm(Spnao[np.ix_(all_basis_indices_nrb,all_basis_indices_nrb)]-np.eye(len(all_basis_indices_nrb))))
    print(sl.norm(Spnao[np.ix_(all_basis_indices_nmb,all_basis_indices_nrb)]))
    '''

    # 4b: Restoration of natural character of the NAOs.
    W.setflags(write=True)
    Ppnao=N.T@P@N
    Nred=np.zeros_like(N)
    basis_indices_both=[]
    for basis_indices_nmb_this_atom,basis_indices_nrb_this_atom in zip(basis_indices_nmb,basis_indices_nrb):
        basis_indices_both.append(basis_indices_nmb_this_atom+basis_indices_nrb_this_atom)
    for basis_indices_this_atom in basis_indices_both:
        lmax=int((max(map(len,basis_indices_this_atom))-1)/2)
        for l in range(lmax+1):
            basis_index_heads=[]
            for basis_index_this_atom_this_shell in basis_indices_this_atom:
                if l==(len(basis_index_this_atom_this_shell)-1)/2:
                    basis_index_heads.append(basis_index_this_atom_this_shell[0])
            PAl=np.zeros([len(basis_index_heads),len(basis_index_heads)])
            SAl=np.zeros([len(basis_index_heads),len(basis_index_heads)])
            for m in range(2*l+1):
                basis_indices=[basis_index+m for basis_index in basis_index_heads]
                PAl+=Ppnao[np.ix_(basis_indices,basis_indices)]
                SAl+=Spnao[np.ix_(basis_indices,basis_indices)]
            PAl/=2*l+1
            SAl/=2*l+1
            OL=np.linalg.inv(sl.sqrtm(SAl))
            PAlL=OL.T@PAl@OL
            w,NL=sl.eigh(PAlL)
            for m in range(2*l+1):
                basis_indices=[basis_index+m for basis_index in basis_index_heads]
                W[np.ix_(basis_indices)]=w[::-1]
                Nred[np.ix_(basis_indices,basis_indices)]=(OL@NL)[:,::-1]
    N=N@Nred

    # Computing population of NAOs
    P_nao=N.T@P@N
    W=np.diag(P_nao)

    return basis_indices_nmb,basis_indices_nrb,W,N

def MinimalShells(an,nc): # an - Atomic number; nc - Nuclear charge
    match an:
        case 1:
            return [1,0,0,0,0,0,0] #  H
        case 2:
            return [1,0,0,0,0,0,0] # He
        case 3:
            return [2,0,0,0,0,0,0] # Li
        case 4:
            return [2,0,0,0,0,0,0] # Be
        case 5:
            return [2,1,0,0,0,0,0] #  B
        case 6:
            return [2,1,0,0,0,0,0] #  C
        case 7:
            return [2,1,0,0,0,0,0] #  N
        case 8:
            return [2,1,0,0,0,0,0] #  O
        case 9:
            return [2,1,0,0,0,0,0] #  F
        case 10:
            return [2,1,0,0,0,0,0] # Ne
        case 11:
            return [3,1,0,0,0,0,0] # Na
        case 12:
            return [3,1,0,0,0,0,0] # Mg
        case 13:
            return [3,2,0,0,0,0,0] # Al
        case 14:
            return [3,2,0,0,0,0,0] # Si
        case 15:
            match nc:
                case 5:
                    return [1,1,0,0,0,0]
            return [3,2,0,0,0,0,0] #  P
        case 16:
            return [3,2,0,0,0,0,0] #  S
        case 17:
            return [3,2,0,0,0,0,0] # Cl
        case 18:
            return [3,2,0,0,0,0,0] # Ar
        case 19:
            match nc:
                case 9:
                    return [2,1,0,0,0,0,0]
            return [4,2,0,0,0,0,0] #  K
        case 20:
            return [4,2,0,0,0,0,0] # Ca
        case 21:
            match nc:
                case 11:
                    return [2,1,1,0,0,0,0]
            return [4,2,1,0,0,0,0] # Sc
        case 22:
            match nc:
                case 12:
                    return [2,1,1,0,0,0,0]
            return [4,2,1,0,0,0,0] # Ti
        case 23:
            match nc:
                case 13:
                    return [2,1,1,0,0,0,0]
            return [4,2,1,0,0,0,0] #  V
        case 24:
            match nc:
                case 14:
                    return [2,1,1,0,0,0,0]
            return [4,2,1,0,0,0,0] # Cr
        case 25:
            match nc:
                case 15:
                    return [2,1,1,0,0,0,0]
            return [4,2,1,0,0,0,0] # Mn
        case 26:
            match nc:
                case 16:
                    return [2,1,1,0,0,0,0]
            return [4,2,1,0,0,0,0] # Fe
        case 27:
            match nc:
                case 17:
                    return [2,1,1,0,0,0,0]
            return [4,2,1,0,0,0,0] # Co
        case 28:
            match nc:
                case 18:
                    return [2,1,1,0,0,0,0]
            return [4,2,1,0,0,0,0] # Ni
        case 29:
            match nc:
                case 19:
                    return [2,1,1,0,0,0,0]
            return [4,2,1,0,0,0,0] # Cu
        case 30:
            match nc:
                case 20:
                    return [2,1,1,0,0,0,0]
            return [4,2,1,0,0,0,0] # Zn
        case 31:
            return [4,3,1,0,0,0,0] # Ga
        case 32:
            return [4,3,1,0,0,0,0] # Ge
        case 33:
            match nc:
                case 5:
                    return [1,1,0,0,0,0,0]
            return [4,3,1,0,0,0,0] # As
        case 34:
            return [4,3,1,0,0,0,0] # Se
        case 35:
            return [4,3,1,0,0,0,0] # Br
        case 36:
            return [4,3,1,0,0,0,0] # Kr
        case 37:
            match nc:
                case 9:
                    return [2,1,0,0,0,0,0]
            return [5,3,1,0,0,0,0] # Rb
        case 38:
            match nc:
                case 10:
                    return [2,1,0,0,0,0,0]
            return [5,3,1,0,0,0,0] # Sr
        case 39:
            return [5,3,2,0,0,0,0] #  Y
        case 40:
            return [5,3,2,0,0,0,0] # Zr
        case 41:
            return [5,3,2,0,0,0,0] # Nb
        case 42:
            return [5,3,2,0,0,0,0] # Mo
        case 43:
            return [5,3,2,0,0,0,0] # Tc
        case 44:
            match nc:
                case 16:
                    return [2,1,1,0,0,0,0]
            return [5,3,2,0,0,0,0] # Ru
        case 45:
            return [5,3,2,0,0,0,0] # Rh
        case 46:
            match nc:
                case 18:
                    return [2,1,1,0,0,0,0]
            return [5,3,2,0,0,0,0] # Pd
        case 47:
            return [5,3,2,0,0,0,0] # Ag
        case 48:
            return [5,3,2,0,0,0,0] # Cd
        case 49:
            return [5,4,2,0,0,0,0] # In
        case 50:
            match nc:
                case 4:
                    return [1,1,0,0,0,0,0]
                case 22:
                    return [2,2,1,0,0,0,0]
            return [5,4,2,0,0,0,0] # Sn
        case 51:
            match nc:
                case 5:
                    return [1,1,0,0,0,0,0]
                case 23:
                    return [2,2,1,0,0,0,0]
            return [5,4,2,0,0,0,0] # Sb
        case 52:
            return [5,4,2,0,0,0,0] # Te
        case 53:
            match nc:
                case 25:
                    return [2,2,1,0,0,0,0]
            return [5,4,2,0,0,0,0] #  I
        case 54:
            return [5,4,2,0,0,0,0] # Xe
        case 55:
            return [6,4,2,0,0,0,0] # Cs
        case 56:
            return [6,4,2,0,0,0,0] # Ba
        case 57:
            return [6,4,3,0,0,0,0] # La
        case 58:
            return [6,4,3,1,0,0,0] # Ce
        case 59:
            return [6,4,2,1,0,0,0] # Pr
        case 60:
            return [6,4,2,1,0,0,0] # Nd
        case 61:
            return [6,4,2,1,0,0,0] # Pm
        case 62:
            return [6,4,2,1,0,0,0] # Sm
        case 63:
            return [6,4,2,1,0,0,0] # Eu
        case 64:
            return [6,4,3,1,0,0,0] # Gd
        case 65:
            return [6,4,2,1,0,0,0] # Tb
        case 66:
            return [6,4,2,1,0,0,0] # Dy
        case 67:
            return [6,4,2,1,0,0,0] # Ho
        case 68:
            return [6,4,2,1,0,0,0] # Er
        case 69:
            return [6,4,2,1,0,0,0] # Tm
        case 70:
            return [6,4,2,1,0,0,0] # Yb
        case 71:
            match nc:
                case 11:
                    return [2,1,1,0,0,0,0]
            return [6,4,3,1,0,0,0] # Lu
        case 72:
            match nc:
                case 12:
                    return [2,1,1,0,0,0,0]
            return [6,4,3,1,0,0,0] # Hf
        case 73:
            match nc:
                case 13:
                    return [2,1,1,0,0,0,0]
            return [6,4,3,1,0,0,0] # Ta
        case 74:
            match nc:
                case 14:
                    return [2,1,1,0,0,0,0]
            return [6,4,3,1,0,0,0] #  W
        case 75:
            match nc:
                case 15:
                    return [2,1,1,0,0,0,0]
            return [6,4,3,1,0,0,0] # Re
        case 76:
            match nc:
                case 16:
                    return [2,1,1,0,0,0,0]
            return [6,4,3,1,0,0,0] # Os
        case 77:
            match nc:
                case 17:
                    return [2,1,1,0,0,0,0]
            return [6,4,3,1,0,0,0] # Ir
        case 78:
            match nc:
                case 18:
                    return [2,1,1,0,0,0,0]
            return [6,4,3,1,0,0,0] # Pt
        case 79:
            match nc:
                case 19:
                    return [2,1,1,0,0,0,0]
            return [6,4,3,1,0,0,0] # Au
        case 80:
            match nc:
                case 20:
                    return [2,1,1,0,0,0,0]
            return [6,4,3,1,0,0,0] # Hg
        case 81:
            return [6,5,3,1,0,0,0] # Tl
        case 82:
            return [6,5,3,1,0,0,0] # Pb
        case 83:
            match nc:
                case 5:
                    return [1,1,0,0,0,0,0]
                case 23:
                    return [2,2,1,0,0,0,0]
            return [6,5,3,1,0,0,0] # Bi
        case 84:
            return [6,5,3,1,0,0,0] # Po
        case 85:
            return [6,5,3,1,0,0,0] # At
        case 86:
            return [6,5,3,1,0,0,0] # Rn
        case 87:
            return [7,5,3,1,0,0,0] # Fr
        case 88:
            return [7,5,3,1,0,0,0] # Ra
        case 89:
            return [7,5,4,1,0,0,0] # Ac
        case 90:
            match nc:
                case 12:
                    return [2,1,1,1,0,0,0]
                case 30:
                    return [3,2,2,1,0,0,0]
            return [7,5,4,1,0,0,0] # Th
        case 91:
            return [7,5,4,2,0,0,0] # Pa
        case 92:
            match nc:
                case 14:
                    return [2,1,1,2,0,0,0]
                case 32:
                    return [3,2,2,2,0,0,0]
            return [7,5,4,2,0,0,0] # U


def NaturalAtomicOrbital(mo_mwfn,modify_minimal_shells=[]):
    nao_mwfn=cp.deepcopy(mo_mwfn)
    basis_indices_by_center=nao_mwfn.getBasisIndexByCenter()
    shell_indices_by_center=nao_mwfn.getShellIndexByCenter()
    basis_indices_by_shell=nao_mwfn.getBasisIndexByShell()
    basis_indices_nmb,basis_indices_nrb=None,None
    if nao_mwfn.Overlap_matrix is None:
        nao_mwfn.calcOverlap()
    S=nao_mwfn.Overlap_matrix
    angulars=[shell.Type for shell in nao_mwfn.Shells]
    minimal_shells=[MinimalShells(center.Index,int(center.Nuclear_charge)) for center in mo_mwfn.Centers]
    for icenter, minimal_shell in modify_minimal_shells:
        minimal_shells[icenter]=minimal_shell
    print("Natural atomic orbitals:")
    if mo_mwfn.Wfntype==0 or mo_mwfn.Wfntype==1:
        for spin in ([0] if mo_mwfn.Wfntype==0 else [1,2]):
            D=None
            match spin:
                case 0:
                    D=nao_mwfn.Total_density_matrix
                case 1:
                    D=nao_mwfn.Alpha_density_matrix
                case 2:
                    D=nao_mwfn.Beta_density_matrix
            basis_indices_nmb,basis_indices_nrb,W,N=generateNaturalAtomicOrbital(shell_indices_by_center,basis_indices_by_shell,basis_indices_by_center,angulars,D,S,minimal_shells)
            nao_mwfn.setOccupation(spin,W)
            nao_mwfn.setCoefficientMatrix(spin,N)
            nao_mwfn.setEnergy(spin,[0 for i in range(nao_mwfn.getNumBasis())])
            match spin:
                case 0:
                    nao_mwfn.Extra_info["NAO_density_matrix"]=N.T@S@D@S@N
                case 1:
                    nao_mwfn.Extra_info["NAO_alpha_density_matrix"]=N.T@S@D@S@N
                case 2:
                    nao_mwfn.Extra_info["NAO_beta_density_matrix"]=N.T@S@D@S@N
    nao_mwfn.Extra_info["NAO_minimal_basis_indices"]=basis_indices_nmb
    nao_mwfn.Extra_info["NAO_Rydberg_basis_indices"]=basis_indices_nrb
    nao_mwfn.Comment="Natural atomic orbitals."
    return nao_mwfn
