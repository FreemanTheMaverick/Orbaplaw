import numpy as np
import copy as cp
from Orbaplaw import Integrals as eint
from . import PipekMezey
from . import FosterBoys
from . import Fock
from . import Orbitalet

calcS='''
if mo_mwfn.Overlap_matrix is None:
	mo_mwfn.calcOverlap()
'''

def Localizer(mo_mwfn, method = "PipekMezey-Lowdin", space = "occ"):

	loc_mwfn = cp.deepcopy(mo_mwfn)
	norbitals = loc_mwfn.getNumIndBasis()
	for spin in ([0] if loc_mwfn.Wfntype == 0 else [1, 2]):
		nocc = round(loc_mwfn.getNumElec(spin))
		if spin == 0:
			nocc //= 2
		orbital_range = []
		if type(space) is str:
			if space.upper() == "OCC":
				orbital_range = list(range(nocc))
			if space.upper() == "VIR":
				orbital_range = list(range(nocc, norbitals))
			if space.upper() == "MIX":
				orbital_range = list(range(norbitals))
		elif type(space) is list:
			orbital_range = space
		Cold = loc_mwfn.getCoefficientMatrix(spin)[:, orbital_range]
		Cnew = Cold.copy()

		if "PM" in method.upper() or "PIPEK" in method.upper() or "MEZEY" in method.upper():
			charge_type = ""
			if "LOWDIN" in method.upper():
				charge_type = "Lowdin"
			elif "MULLIKEN" in method.upper():
				charge_type = "Mulliken"
			else:
				raise RuntimeError("Unrecognized charge type!")
			if space.upper() == "MIX":
				raise RuntimeError("Fractionally occupied orbitals and mixing occupied and virtual orbitals in Pipek-Mezey localization is not supported!")
			print("Pipek-Mezey localization (%s) on Spin %d Orbitals %d - %d:" % ( charge_type, spin, orbital_range[0], orbital_range[-1] ))
			if mo_mwfn.Overlap_matrix is None:
				mo_mwfn.calcOverlap()
			S = mo_mwfn.Overlap_matrix
			basis_indices_by_center = loc_mwfn.getBasisIndexByCenter()
			Cnew = Cold @ PipekMezey(Cold, S, basis_indices_by_center, charge_type)

		elif "FB" in method.upper() or "FOSTER" in method.upper() or "BOYS" in method.upper():
			print("Foster-Boys localization on Spin %d Orbitals %d - %d:" % (spin, orbital_range[0], orbital_range[-1]))
			X, Y, Z = eint.PyscfDipole(loc_mwfn, loc_mwfn)
			XX, _, _, YY, _, ZZ = eint.PyscfQuadrupole(loc_mwfn, loc_mwfn)
			Waos = [ X, Y, Z ]
			W2aoSum = - XX - YY - ZZ
			Cnew = Cold @ FosterBoys(Cold, Waos, W2aoSum)

		C = loc_mwfn.getCoefficientMatrix(spin)
		C[:, orbital_range] = Cnew
		loc_mwfn.setCoefficientMatrix(spin, C)

	return loc_mwfn
