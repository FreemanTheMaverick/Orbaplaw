import numpy as np

def OrbitalAlignment(A,B,epsilon,diagmat=True,diagmis=True):
    # A - the coefficient matrix of the whole molecule
    # B - the coefficient matrix of the fragment
    # epsilon - the orbital energies of the whole molecule
    # diagmat - whether to diagonalize the matched space
    # diagmis - whether to diagonalize the mismatched space
    m=A.shape[0] # The number of A's basis functions
    n=A.shape[1] # The number of A's occupied orbitals
    p=B.shape[0] # The number of B's basis functions
    q=B.shape[1] # The number of B's occupied orbitals
    # m > n, p > q.
    # In the case of the same cutoff and the same lattice parameters, m = p and n > q.
    ''' # The overlap matrix is identity for plane wave basis sets.
    S=np.zeros([p,m])
    for i in range(min(p,m)):
        S[i,i]=1
    '''
    I=np.zeros([q,n])
    for i in range(min(q,n)):
        I[i,i]=1.
    U,Sigma,VH=np.linalg.svd((B.conj().T@A).conj().T@I)
    X=U@VH

    E=np.zeros(n)
    Y=np.eye(n,dtype='complex128')
    Fo=X.conj().T@(epsilon if epsilon.ndim==2 else np.diag(epsilon))@X
    if diagmat:
        Fmat=Fo[0:q,0:q]
        E[0:q],Y[0:q,0:q]=np.linalg.eigh(Fmat)
    if diagmis:
        Fmis=Fo[q:n,q:n]
        E[q:n],Y[q:n,q:n]=np.linalg.eigh(Fmis)
    return A@X@Y,Sigma,E
