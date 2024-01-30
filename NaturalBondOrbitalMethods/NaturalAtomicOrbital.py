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
    for iatom in range(natoms):
        W_nmb_atom=np.array([])
        N_nmb_atom=[]
        W_nrb_atom=np.array([])
        N_nrb_atom=[]
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
                else: # NRB
                    for m in range(2*l+1):
                        W_nrb_atom=np.append(W_nrb_atom,WW[basis_indices_by_shell[shell_indices[args[p]]][0]+m])
                        N_nrb_atom.append(NN[:,basis_indices_by_shell[shell_indices[args[p]]][0]+m])
        W_atom=np.append(W_nmb_atom,W_nrb_atom)
        W=np.append(W,W_atom)
        N_nmb_atom=np.array(N_nmb_atom)
        N_nrb_atom=np.array(N_nrb_atom)
        N.extend(N_nmb_atom)
        N.extend(N_nrb_atom)
    N=np.array(N).T
    return W,N

    # 3b: Schmidt interatomic orthogonalization of NRB to NMB orbitals.
    # 3c: Restoration of natural character of the NRB.
    # 4: Formation of the final NAO set
    # 4a: Weighted interatomic orthogonalization within NMB and within NRB.
    # 4b: Resroration of natural character of the NAOs.




