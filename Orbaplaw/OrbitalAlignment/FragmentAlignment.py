import numpy as np
import copy as cp
from Orbaplaw import Integrals as eint
from . import OrbitalAlignment


def FragmentAlignment(big_mwfn, small_mwfn_list, diagmat = False, diagmis = True):
	famo_mwfn = cp.deepcopy(big_mwfn)
	overlap_matrices = []
	for small_mwfn in small_mwfn_list:
		overlap_matrices.append(eint.PyscfOverlap(big_mwfn,small_mwfn))
	S = np.hstack(overlap_matrices)
	for spin in ([0] if big_mwfn.Wfntype == 0 else [1, 2]):
		noccA = round(big_mwfn.getNumElec(spin) / (2 if big_mwfn.Wfntype == 0 else 1))
		A = big_mwfn.getCoefficientMatrix(spin)[:, :noccA]
		epsilon = big_mwfn.getEnergy(spin)[:noccA]
		nbasisB = 0
		B_basis_heads = []
		B_basis_tails = []
		noccB = 0
		B_orbital_heads = []
		B_orbital_tails = []
		for small_mwfn in small_mwfn_list:
			B_basis_heads.append(nbasisB)
			nbasisB += small_mwfn.getNumBasis()
			B_basis_tails.append(nbasisB)
			B_orbital_heads.append(noccB)
			noccB += round(small_mwfn.getNumElec(0 if small_mwfn.Wfntype == 0 else spin) / (2 if small_mwfn.Wfntype == 0 else 1))
			B_orbital_tails.append(noccB)
		B = np.zeros([nbasisB, noccB])
		for i, small_mwfn in enumerate(small_mwfn_list):
			B_basis_head = B_basis_heads[i]
			B_basis_tail = B_basis_tails[i]
			B_orbital_head = B_orbital_heads[i]
			B_orbital_tail = B_orbital_tails[i]
			noccB = B_orbital_tail - B_orbital_head
			B[B_basis_head : B_basis_tail, B_orbital_head : B_orbital_tail] = small_mwfn.getCoefficientMatrix(0 if small_mwfn.Wfntype == 0 else spin)[:, :noccB]
		C = big_mwfn.getCoefficientMatrix(spin)
		E = big_mwfn.getEnergy(spin)
		C[:, :noccA], E[:noccA] = OrbitalAlignment(A, B, S, epsilon, diagmat, diagmis)
		famo_mwfn.setCoefficientMatrix(spin, C)
		famo_mwfn.setEnergy(spin, E)
	return famo_mwfn
