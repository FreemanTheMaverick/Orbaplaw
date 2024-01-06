import numpy as np

def OrbitalAlignment(A,B): # A - the whole molecule; B - the fragment
    m=A.shape[0] # The number of A's basis functions
    n=A.shape[1] # The number of A's occupied orbitals
    p=B.shape[0] # The number of B's basis functions
    q=B.shape[1] # The number of B's occupied orbitals
    # m > n, p > q.
    # In the case of the same cutoff and the same lattice parameters, m = p and n > q
    ''' # The overlap matrix is identity for plane wave basis sets.
    S=np.zeros([p,m])
    for i in range(min(p,m)):
        S[i,i]=1
    '''
    I=np.zeros([q,n])
    for i in range(min(q,n)):
        I[i,i]=1.
    U,Sigma,VH=np.linalg.svd((B.conj().T@A).T@I)
    return A@U@VH,Sigma
