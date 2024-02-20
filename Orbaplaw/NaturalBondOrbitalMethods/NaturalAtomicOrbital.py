import numpy as np
import scipy.optimize as so
import scipy.special as ss
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
        lmax=int((max(map(len,basis_indices_this_atom))-1)/2)
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

    return basis_indices_nmb,basis_indices_nrb,W,N


MinimalShells=[
        [0,0,0,0,0,0,0],
        [1,0,0,0,0,0,0], #  H
        [1,0,0,0,0,0,0], # He
        [2,0,0,0,0,0,0], # Li
        [2,0,0,0,0,0,0], # Be
        [2,1,0,0,0,0,0], #  B
        [2,1,0,0,0,0,0], #  C
        [2,1,0,0,0,0,0], #  N
        [2,1,0,0,0,0,0], #  O
        [2,1,0,0,0,0,0], #  F
        [2,1,0,0,0,0,0], # Ne
        [3,1,0,0,0,0,0], # Na
        [3,1,0,0,0,0,0], # Mg
        [3,2,0,0,0,0,0], # Al
        [3,2,0,0,0,0,0], # Si
        [3,2,0,0,0,0,0], #  P
        [3,2,0,0,0,0,0], #  S
        [3,2,0,0,0,0,0], # Cl
        [3,2,0,0,0,0,0], # Ar
        [4,2,0,0,0,0,0], #  K
        [4,2,0,0,0,0,0], # Ca
        [4,2,1,0,0,0,0], # Sc
        [4,2,1,0,0,0,0], # Ti
        [4,2,1,0,0,0,0], #  V
        [4,2,1,0,0,0,0], # Cr
        [4,2,1,0,0,0,0], # Mn
        [4,2,1,0,0,0,0], # Fe
        [4,2,1,0,0,0,0], # Co
        [4,2,1,0,0,0,0], # Ni
        [4,2,1,0,0,0,0]] # Cu

def NaturalAtomicOrbital(mwfn_obj):
    result_mwfn_obj=cp.deepcopy(mwfn_obj)
    basis_indices_by_center=result_mwfn_obj.getBasisIndexByCenter()
    shell_indices_by_center=result_mwfn_obj.getShellIndexByCenter()
    basis_indices_by_shell=result_mwfn_obj.getBasisIndexByShell()
    D=result_mwfn_obj.Total_density_matrix/2
    S=result_mwfn_obj.Overlap_matrix
    angulars=[shell.Type for shell in result_mwfn_obj.Shells]
    minimal_shells=[MinimalShells[int(center.Nuclear_charge)] for center in mwfn_obj.Centers]
    basis_indices_nmb,basis_indices_nrb,W,N=generateNaturalAtomicOrbital(shell_indices_by_center,basis_indices_by_shell,basis_indices_by_center,angulars,D,S,minimal_shells)
    result_mwfn_obj.setOccupation(W)
    result_mwfn_obj.setCoefficientMatrix(N)
    result_mwfn_obj.setEnergy([0 for iorbital in result_mwfn_obj.Orbitals])
    result_mwfn_obj.Extra_info["NAO_density_matrix"]=N.T@S@D@S@N
    result_mwfn_obj.Extra_info["NAO_minimal_basis_indices"]=basis_indices_nmb
    result_mwfn_obj.Extra_info["NAO_Rydberg_basis_indices"]=basis_indices_nrb
    result_mwfn_obj.Comment="Natural atomic orbitals. The Total_density_matrix is NAO-based P."
    return result_mwfn_obj

