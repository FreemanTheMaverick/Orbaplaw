import numpy as np
from pyscf import gto


def PyscfDecompose(mwfn):
	atom_string = ""
	for iatom, atom in enumerate(mwfn.Centers):
		atom_string += atom.Symbol + str(iatom) + " " + np.array2string(atom.Coordinates)[1:-1] + "\n"
	mol_list = []
	cart_list = []
	head_list = []
	tail_list = []
	center_list = mwfn.getShellCenters()
	nbasis = 0
	charge = round(mwfn.getCharge())
	spin = round(mwfn.getSpin())
	for ishell in range(mwfn.getNumShells()):
		shell = mwfn.Shells[ishell]
		mol = gto.Mole()
		mol_list.append(mol)
		cart_list.append(shell.Type >= 2)
		head_list.append(nbasis)
		nbasis += shell.getSize()
		tail_list.append(nbasis)
		jcenter = center_list[ishell]
		mol.atom = atom_string
		mol.basis = {
				shell.Center.Symbol + str(jcenter):[ [abs(shell.Type)] + [(shell.Exponents[j], shell.Coefficients[j]) for j in range(shell.getNumPrims())] ]
		}
		mol.unit = 'B'
		mol.charge = charge
		mol.spin = spin
		mol.build()
	return mol_list, cart_list, head_list, tail_list

def PyscfOverlap(mwfn1,mwfn2):
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
			assert icart == jcart
			overlap[ihead:itail, jhead:jtail] = gto.intor_cross("int1e_ovlp_" + ("cart" if icart else "sph"), imol, jmol)
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
		assert icart == jcart
__replacement__
'''

def PyscfOverlap(mwfn1, mwfn2):
	overlap = np.zeros([mwfn1.getNumBasis(), mwfn2.getNumBasis()])
	replacement = '''
		overlap[ihead:itail, jhead:jtail] = gto.intor_cross("int1e_ovlp_" + ("cart" if icart else "sph"), imol, jmol)
	'''
	exec(nagging.replace("__replacement__", replacement))
	return overlap

def PyscfDipole(mwfn1, mwfn2):
	X = np.zeros([mwfn1.getNumBasis(), mwfn2.getNumBasis()])
	Y = np.zeros([mwfn1.getNumBasis(), mwfn2.getNumBasis()])
	Z = np.zeros([mwfn1.getNumBasis(), mwfn2.getNumBasis()])
	replacement = '''
		X[ihead:itail, jhead:jtail], Y[ihead:itail, jhead:jtail], Z[ihead:itail, jhead:jtail] = gto.intor_cross("int1e_r_" + ("cart" if icart else "sph"), imol, jmol, comp=3)
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
		XX[ihead:itail, jhead:jtail], XY[ihead:itail, jhead:jtail], XZ[ihead:itail, jhead:jtail], _, YY[ihead:itail, jhead:jtail], YZ[ihead:itail, jhead:jtail], _, _, ZZ[ihead:itail, jhead:jtail] = gto.intor_cross("int1e_rr_" + ("cart" if icart else "sph"), imol, jmol, comp=9)
	'''
	exec(nagging.replace("__replacement__", replacement))
	return XX, XY, XZ, YY, YZ, ZZ
