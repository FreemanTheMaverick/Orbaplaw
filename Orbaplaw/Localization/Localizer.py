import numpy as np
import copy as cp
from Orbaplaw import Integrals as eint
from Orbaplaw import Population as pop
from Orbaplaw import Miscellany as mis
from . import PipekMezey
from . import FosterBoys
from . import Fock
from . import Orbitalet


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
		format_range = mis.FormatRange(orbital_range)
		Cold = loc_mwfn.getCoefficientMatrix(spin)[:, orbital_range]
		Cnew = Cold.copy()

		if "PM" in method.upper() or "PIPEK" in method.upper() or "MEZEY" in method.upper():
			if mo_mwfn.Overlap_matrix is None:
				mo_mwfn.calcOverlap()
			S = mo_mwfn.Overlap_matrix
			basis_indices_by_center = loc_mwfn.getBasisIndexByCenter()
			charge_type = ""
			Qrefs = []
			if "LOWDIN" in method.upper():
				charge_type = "Lowdin"
				Qrefs = pop.Lowdin(Cold, S, basis_indices_by_center)
			elif "MULLIKEN" in method.upper():
				charge_type = "Mulliken"
				Qrefs = pop.Mulliken(Cold, S, basis_indices_by_center)
			else:
				raise RuntimeError("Unrecognized charge type!")
			print("Pipek-Mezey localization (%s) on Spin %d Orbitals %s:" % ( charge_type, spin, format_range) )
			occ_all = mo_mwfn.getOccupation(spin)
			occ = [occ_all[i] for i in orbital_range]
			mix = len(set(occ)) > 1
			if mix:
				raise RuntimeError("Fractionally occupied orbitals and mixing occupied and virtual orbitals in Pipek-Mezey localization is not supported!")
			Cnew = Cold @ PipekMezey(Qrefs)

		elif "FB" in method.upper() or "FOSTER" in method.upper() or "BOYS" in method.upper():
			print("Foster-Boys localization on Spin %d Orbitals %s:" % (spin, format_range))
			X, Y, Z = eint.PyscfDipole(loc_mwfn, loc_mwfn)
			XX, _, _, YY, _, ZZ = eint.PyscfQuadrupole(loc_mwfn, loc_mwfn)
			Waos = [ X, Y, Z ]
			Wrefs = [ Cold.T @ Wao @ Cold for Wao in Waos ]
			Wao2Sum = - XX - YY - ZZ
			W2refSum = Cold.T @ Wao2Sum @ Cold
			Cnew = Cold @ FosterBoys(Wrefs, W2refSum)

		elif "FOCK" in method.upper(): # Must start from CMOs
			print("Fock localization on Spin %d Orbitals %s:" % (spin, format_range))
			Eref = np.array(mo_mwfn.getEnergy(spin))[orbital_range]
			Cnew = Cold @ Fock(Eref)

		elif "OL" in method.upper() or "ORBITALET" in method.upper(): # Must start from CMOs
			gamma_e = 0.7959
			for word in method.split():
				try:
					gamma_e = float(word)
					break
				except ValueError:
					continue
			if gamma_e <= 0 or gamma_e >= 1:
				raise RuntimeError("Gamma_e must be in (0, 1)!")
			print("Orbitalet localization (Gamma_e = %f) on Spin %d Orbitals %s:" % (gamma_e, spin, format_range))
			X, Y, Z = eint.PyscfDipole(loc_mwfn, loc_mwfn)
			XX, _, _, YY, _, ZZ = eint.PyscfQuadrupole(loc_mwfn, loc_mwfn)
			Waos = [ X, Y, Z ]
			Wrefs = [ Cold.T @ Wao @ Cold for Wao in Waos ]
			Wao2Sum = - XX - YY - ZZ
			W2refSum = Cold.T @ Wao2Sum @ Cold
			Eref = np.array(mo_mwfn.getEnergy(spin))[orbital_range]
			Cnew = Cold @ Orbitalet(Wrefs, W2refSum, Eref, gamma_e)

		C = loc_mwfn.getCoefficientMatrix(spin)
		C[:, orbital_range] = Cnew
		loc_mwfn.setCoefficientMatrix(spin, C)

	return loc_mwfn
