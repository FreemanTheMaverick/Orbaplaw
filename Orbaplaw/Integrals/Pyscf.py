import numpy as np
from pyscf import gto


def PyscfDecompose(mwfn):
	atom_string = ""
	for iatom, atom in enumerate(mwfn.Centers):
		atom_string += atom.getSymbol() + str(iatom) + " " + str(atom.Coordinates[0]) + " " + str(atom.Coordinates[1]) + " " + str(atom.Coordinates[2]) + "\n"
	mol_list = []
	cart_list = []
	head_list = []
	tail_list = []
	center_list = mwfn.Shell2Atom()
	nbasis = 0
	charge = round(mwfn.getCharge())
	spin = round(mwfn.getNumElec(1) - mwfn.getNumElec(2))
	shell2atom = mwfn.Shell2Atom()
	for ishell in range(mwfn.getNumShells()):
		shell = mwfn.getShell(ishell)
		mol = gto.Mole()
		mol_list.append(mol)
		cart_list.append(shell.Type >= 2)
		head_list.append(nbasis)
		nbasis += shell.getSize()
		tail_list.append(nbasis)
		jcenter = shell2atom[ishell]
		mol.atom = atom_string
		mol.basis = {
				mwfn.Centers[jcenter].getSymbol() + str(jcenter):[ [abs(shell.Type)] + [(shell.Exponents[j], shell.Coefficients[j]) for j in range(shell.getNumPrims())] ]
		}
		mol.unit = 'B'
		mol.charge = charge
		mol.spin = spin
		mol.build()
	return mol_list, cart_list, head_list, tail_list

def PyscfOverlap(mwfn1, mwfn2):
	overlap = np.zeros([mwfn1.getNumBasis(), mwfn2.getNumBasis()])
	mol_list1, cart_list1, head_list1, tail_list1 = PyscfDecompose(mwfn1)
	mol_list2, cart_list2, head_list2, tail_list2 = PyscfDecompose(mwfn2)
	for ishell in range(len(mol_list1)):
		imol = mol_list1[ishell]
		icart = cart_list1[ishell]
		ihead = head_list1[ishell]
		itail = tail_list1[ishell]
		for jshell in range(len(mol_list2)):
			jmol = mol_list2[jshell]
			jcart = cart_list2[jshell]
			jhead = head_list2[jshell]
			jtail = tail_list2[jshell]
			if icart or jcart:
				raise RuntimeError("Currently the integrals among cartesian basis functions (l >= 2) are not supported!")
			overlap[ihead:itail, jhead:jtail] = gto.intor_cross("int1e_ovlp_sph", imol, jmol)
	return overlap

# Below is an equivalent but more flexible realization of the function above.

nagging='''
mol_list1, cart_list1, head_list1, tail_list1 = PyscfDecompose(mwfn1)
mol_list2, cart_list2, head_list2, tail_list2 = PyscfDecompose(mwfn2)
for ishell in range(len(mol_list1)):
	imol = mol_list1[ishell]
	icart = cart_list1[ishell]
	ihead = head_list1[ishell]
	itail = tail_list1[ishell]
	for jshell in range(len(mol_list2)):
		jmol = mol_list2[jshell]
		jcart = cart_list2[jshell]
		jhead = head_list2[jshell]
		jtail = tail_list2[jshell]
		if icart or jcart:
			raise RuntimeError("Currently the integrals among cartesian basis functions (l >= 2) are not supported!")
__replacement__
'''

def PyscfOverlap(mwfn1, mwfn2):
	overlap = np.zeros([mwfn1.getNumBasis(), mwfn2.getNumBasis()])
	replacement = '''
		overlap[ihead:itail, jhead:jtail] = gto.intor_cross("int1e_ovlp_sph", imol, jmol)
	'''
	exec(nagging.replace("__replacement__", replacement))
	return overlap

def PyscfDipole(mwfn1, mwfn2):
	X = np.zeros([mwfn1.getNumBasis(), mwfn2.getNumBasis()])
	Y = np.zeros([mwfn1.getNumBasis(), mwfn2.getNumBasis()])
	Z = np.zeros([mwfn1.getNumBasis(), mwfn2.getNumBasis()])
	replacement = '''
		X[ihead:itail, jhead:jtail], Y[ihead:itail, jhead:jtail], Z[ihead:itail, jhead:jtail] = gto.intor_cross("int1e_r_sph", imol, jmol, comp=3)
	'''
	exec(nagging.replace("__replacement__", replacement))
	return X, Y, Z

def PyscfQuadrupole(mwfn1, mwfn2):
	XX = np.zeros([mwfn1.getNumBasis(), mwfn2.getNumBasis()])
	XY = np.zeros([mwfn1.getNumBasis(), mwfn2.getNumBasis()])
	XZ = np.zeros([mwfn1.getNumBasis(), mwfn2.getNumBasis()])
	YY = np.zeros([mwfn1.getNumBasis(), mwfn2.getNumBasis()])
	YZ = np.zeros([mwfn1.getNumBasis(), mwfn2.getNumBasis()])
	ZZ = np.zeros([mwfn1.getNumBasis(), mwfn2.getNumBasis()])
	replacement = '''
		XX[ihead:itail, jhead:jtail], XY[ihead:itail, jhead:jtail], XZ[ihead:itail, jhead:jtail], _, YY[ihead:itail, jhead:jtail], YZ[ihead:itail, jhead:jtail], _, _, ZZ[ihead:itail, jhead:jtail] = gto.intor_cross("int1e_rr_sph", imol, jmol, comp=9)
	'''
	exec(nagging.replace("__replacement__", replacement))
	return XX, XY, XZ, YY, YZ, ZZ
