import numpy as np
import copy as cp
from Orbaplaw import Integrals as eint
from . import OrbitalAlignment


def SpinAlignment(mo_mwfn, diagmat = False, diagmis = True):
	sno_mwfn = cp.deepcopy(mo_mwfn)
	if mo_mwfn.Overlap_matrix is None:
		mo_mwfn.calcOverlap()
	assert(sno_mwfn.Wfntype == 1)
	S = mo_mwfn.Overlap_matrix
	noccA = round(sno_mwfn.getNumElec(1))
	A = sno_mwfn.getCoefficientMatrix(1)
	eA = sno_mwfn.getEnergy(1)
	noccB = round(sno_mwfn.getNumElec(2))
	B = sno_mwfn.getCoefficientMatrix(2)
	C = np.zeros_like(sno_mwfn.getCoefficientMatrix(1))
	A[:, :noccA], eA[:noccA] = OrbitalAlignment(A[:, :noccA], B[:, :noccB], S, eA[:noccA], diagmat, diagmis)
	A[noccA:, :] = 0
	eA[noccA:] = [0] * (len(eA) - noccA)
	sno_mwfn.setCoefficientMatrix(1, A)
	sno_mwfn.setEnergy(1, eA)
	return sno_mwfn
