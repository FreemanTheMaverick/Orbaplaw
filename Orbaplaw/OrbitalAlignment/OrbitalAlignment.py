import numpy as np

def OrbitalAlignment(A, B, S, epsilon, diagmat, diagmis):
	# A - the coefficient matrix of the whole molecule
	# B - the coefficient matrices of the fragments
	# S - the overlap matrix of the two basis sets, row -> basis of A, column -> basis of B
	# epsilon - the orbital energies of the whole molecule
	# diagmat - whether to diagonalize the matched space
	# diagmis - whether to diagonalize the mismatched space
	m = A.shape[0] # The number of A's basis functions
	n = A.shape[1] # The number of A's occupied orbitals
	p = B.shape[0] # The number of B's basis functions
	q = B.shape[1] # The number of B's occupied orbitals
	# m > n, p > q.
	# In the case of the same cutoff and the same lattice parameters, m = p and n > q.
	I = np.zeros([q, n])
	for i in range(min(q, n)):
		I[i, i] = 1.
		U, Sigma_or_E, VH = np.linalg.svd((A.conj().T @ B @ I) if S is None else (A.conj().T @ S @ B @ I))
	X = U @ VH
	C = A @ X
	Sigma_or_E = np.zeros(n)
	SS = C.conj().T @ B if S is None else C.conj().T @ S @ B
	for i in range(q):
		Sigma_or_E[i] = np.max(SS[i, :])

	Y = np.eye(n)
	if epsilon is not None:
		Fo = X.conj().T @ np.diag(epsilon) @ X
		if diagmat:
			Fmat = Fo[0 : q, 0 : q]
			Sigma_or_E[0 : q], Y[0 : q, 0 : q] = np.linalg.eigh(Fmat)
		if diagmis:
			Fmis = Fo[q : n, q : n]
			Sigma_or_E[q : n], Y[q : n, q : n] = np.linalg.eigh(Fmis)
	return A @ X @ Y, Sigma_or_E
