import numpy as np
from Orbaplaw import Miscellany as mis
from . import Lowdin
#from . import Mulliken

def PopulationAnalyzer(mo_mwfn, method = "Lowdin", space = "occ"):
	natoms = mo_mwfn.getNumCenters()
	if mo_mwfn.Overlap_matrix is None:
		mo_mwfn.calcOverlap()
	S = mo_mwfn.Overlap_matrix
	basis_indices_by_center = mo_mwfn.getBasisIndexByCenter()
	for spin in ([0] if mo_mwfn.Wfntype == 0 else [1, 2]):
		nocc = round(mo_mwfn.getNumElec(spin))
		if spin == 0:
			nocc //= 2
		orbital_range = []
		if type(space) is str:
			if space.upper() == "OCC":
				orbital_range = list(range(nocc))
		elif type(space) is list:
			orbital_range = space
		format_range = mis.FormatRange(orbital_range)
		C = mo_mwfn.getCoefficientMatrix(spin)[:, orbital_range]

		charge_type = ""
		Qs = []
		if "LOWDIN" in method.upper():
			charge_type = "Lowdin"
			Qs = Lowdin(C, S, basis_indices_by_center)
		elif "MULLIKEN" in method.upper():
			charge_type = "Mulliken"
			Qs = Mulliken(C, S, basis_indices_by_center)
		else:
			raise RuntimeError("Unrecognized charge type!")

		print("%s population on Spin %d Orbitals %s:" % ( charge_type, spin, format_range ))
		print("--------------------------------------------")
		print("| Index | Symbol | Population | Net Charge |")
		for iatom in range(natoms):
			symbol = mo_mwfn.Centers[iatom].Symbol
			population = np.sum(np.diag(Qs[iatom]))
			if spin == 0:
				population *= 2
			charge = mo_mwfn.Centers[iatom].Nuclear_charge - population
			print("| %5d | %6s | % 10f | % 10f |" % ( iatom, symbol, population, charge))
		print("--------------------------------------------")
