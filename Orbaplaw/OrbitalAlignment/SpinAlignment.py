import numpy as np
import copy as cp
from Orbaplaw import Integrals as eint
from . import OrbitalAlignment


def SpinAlignment(mo_mwfn, diagmat = False, diagmis = True):
	if mo_mwfn.Overlap.shape != tuple([mo_mwfn.getNumBasis()] * 2):
		mo_mwfn.Overlap = eint.PyscfOverlap(mo_mwfn, mo_mwfn)
	sno_mwfn = mo_mwfn.Clone()
	if sno_mwfn.Wfntype != 1:
		raise RuntimeError("Spin alignment can only be done on spin-unrestricted wavefunctions!")
	S = sno_mwfn.Overlap
	noccA = round(sno_mwfn.getNumElec(1))
	A = sno_mwfn.getCoefficientMatrix(1)
	eA = sno_mwfn.getEnergy(1)
	noccB = round(sno_mwfn.getNumElec(2))
	B = sno_mwfn.getCoefficientMatrix(2)
	C = np.zeros_like(sno_mwfn.getCoefficientMatrix(1))
	if noccA < noccB:
		raise RuntimeError("The first spin type must have no fewer electrons than the second!")
	A[:, :noccA], eA[:noccA] = OrbitalAlignment(A[:, :noccA], B[:, :noccB], S, eA[:noccA], diagmat, diagmis)
	A[noccA:, :] = 0
	eA[noccA:] = [0] * (len(eA) - noccA)
	sno_mwfn.setCoefficientMatrix(A, 1)
	sno_mwfn.setEnergy(eA, 1)
	return sno_mwfn
