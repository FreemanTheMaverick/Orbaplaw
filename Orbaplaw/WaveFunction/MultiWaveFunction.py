import numpy as np
import scipy.linalg as sl
from Orbaplaw import Integrals as eint


__angstrom2bohr__ = 1.8897259886

class MwfnCenter:
	Symbol = None # String
	Index = None # Int
	Nuclear_charge = None # Float
	Coordinates = None # np.array (1,3)
	Shells = None # MwfnShell
	def getNumShells(self):
		return len(self.Shells)

class MwfnShell:
	Type = None # int
	Center = None # MwfnCenter
	Exponents = None # List
	Coefficients = None # List
	def getSize(self):
		if self.Type >= 0:
			return ( self.Type + 1 ) * ( self.Type + 2 ) // 2
		else:
			return 2 * abs(self.Type) + 1
	def getNumPrims(self):
		return len(self.Exponents)
	
class MwfnOrbital:
	Type = None
	Energy = None
	Occ = None
	Sym = None
	Coeff = None # np.array (Nshell,1)

class MultiWaveFunction:

	Comment=""

	# Field 1
	Wfntype=None
	E_tot=None
	VT_ratio=None

	# Field 2
	Centers=None # List of MwfnCenter

	# Field 3
	Shells=None # List of MwfnShell

	# Field 4
	Orbitals=None # List of np.array

	# Field 5
	Total_density_matrix=None # np.array
	Alpha_density_matrix=None
	Beta_density_matrix=None
	Hamiltonian_matrix=None
	Alpha_Hamiltonian_matrix=None
	Beta_Hamiltonian_matrix=None
	Overlap_matrix=None
	Kinetic_energy_matrix=None
	Potential_energy_matrix=None
	X_electric_dipole_moment_matrix=None
	Y_electric_dipole_moment_matrix=None
	Z_electric_dipole_moment_matrix=None
	XX_electric_quadrupole_moment_matrix=None
	YY_electric_quadrupole_moment_matrix=None
	ZZ_electric_quadrupole_moment_matrix=None
	XY_electric_quadrupole_moment_matrix=None
	YZ_electric_quadrupole_moment_matrix=None
	XZ_electric_quadrupole_moment_matrix=None

	Extra_info={}

	def MatrixTransform(self):
		SPDFGHI=dict()
		SPDFGHI[0]=np.array([[1]])
		SPDFGHI[1]=np.eye(3)
		SPDFGHI[-2]=np.array([
			[0,0,1,0,0],
			[0,0,0,1,0],
			[0,1,0,0,0],
			[0,0,0,0,1],
			[1,0,0,0,0]])
		SPDFGHI[-3]=np.array([
			[0,0,0,1,0,0,0],
			[0,0,0,0,1,0,0],
			[0,0,1,0,0,0,0],
			[0,0,0,0,0,1,0],
			[0,1,0,0,0,0,0],
			[0,0,0,0,0,0,1],
			[1,0,0,0,0,0,0]])
		SPDFGHI[-4]=np.array([
			[0,0,0,0,1,0,0,0,0],
			[0,0,0,0,0,1,0,0,0],
			[0,0,0,1,0,0,0,0,0],
			[0,0,0,0,0,0,1,0,0],
			[0,0,1,0,0,0,0,0,0],
			[0,0,0,0,0,0,0,1,0],
			[0,1,0,0,0,0,0,0,0],
			[0,0,0,0,0,0,0,0,1],
			[1,0,0,0,0,0,0,0,0]])
		SPDFGHI[-5]=np.array([
			[0,0,0,0,0,1,0,0,0,0,0],
			[0,0,0,0,0,0,1,0,0,0,0],
			[0,0,0,0,1,0,0,0,0,0,0],
			[0,0,0,0,0,0,0,1,0,0,0],
			[0,0,0,1,0,0,0,0,0,0,0],
			[0,0,0,0,0,0,0,0,1,0,0],
			[0,0,1,0,0,0,0,0,0,0,0],
			[0,0,0,0,0,0,0,0,0,1,0],
			[0,1,0,0,0,0,0,0,0,0,0],
			[0,0,0,0,0,0,0,0,0,0,1],
			[1,0,0,0,0,0,0,0,0,0,0]])
		SPDFGHI[-6]=np.array([
			[0,0,0,0,0,0,1,0,0,0,0,0,0],
			[0,0,0,0,0,0,0,1,0,0,0,0,0],
			[0,0,0,0,0,1,0,0,0,0,0,0,0],
			[0,0,0,0,0,0,0,0,1,0,0,0,0],
			[0,0,0,0,1,0,0,0,0,0,0,0,0],
			[0,0,0,0,0,0,0,0,0,1,0,0,0],
			[0,0,0,1,0,0,0,0,0,0,0,0,0],
			[0,0,0,0,0,0,0,0,0,0,1,0,0],
			[0,0,1,0,0,0,0,0,0,0,0,0,0],
			[0,0,0,0,0,0,0,0,0,0,0,1,0],
			[0,1,0,0,0,0,0,0,0,0,0,0,0],
			[0,0,0,0,0,0,0,0,0,0,0,0,1],
			[1,0,0,0,0,0,0,0,0,0,0,0,0]])
		nbasis=self.getNumBasis()
		mwfntransform=np.zeros([nbasis,nbasis])
		jbasis=0
		for center in self.Centers:
			for shell in center.Shells:
				l=shell.Type
				num=shell.getSize()
				mwfntransform[jbasis:jbasis+num,jbasis:jbasis+num]=SPDFGHI[l]
				jbasis+=num
		return mwfntransform

	def __init__(self,filename):
		def ReadMatrix(f,nrows,ncols,lower):
			matrix=np.zeros([nrows,ncols])
			assert nrows==ncols
			nelements=((1+nrows)*nrows//2) if lower else nrows*ncols
			elements=["114514" for i in range(nelements)]
			ielement=0
			finished=False
			while not finished:
				newline=f.readline()
				if not newline:
					break
				newvalues=newline.split()
				for newvalue in newvalues:
					try:
						elements[ielement]=float(newvalue)
					except ValueError:
						elements[ielement]=0.
					ielement+=1
					if ielement==nelements:
						finished=True
						break
			matrix=np.zeros([nrows,ncols])
			ielement=0
			for irow in range(nrows):
				for jcol in range(irow+1 if lower else ncols):
					matrix[irow,jcol]=elements[ielement]
					if lower:
						matrix[jcol,irow]=elements[ielement]
					ielement+=1
			return matrix

		with open(filename,'r') as f:
			nbasis=-114
			nindbasis=-114
			mwfntransform=None
			tmp_int=None
			while True:
				line=f.readline()
				if not line:
					break
				values=[]
				value=None
				if len(line.split())>1:
					values=line.split()
					value=values[1]

				if "# Comment: " in line:
					self.Comment=line[11:]
				
				# Field 1
				elif "Wfntype=" in line:
					self.Wfntype=int(value)
				elif "E_tot=" in line:
					self.E_tot=float(value)
				elif "VT_ratio=" in line:
					self.VT_ratio=float(value)

				# Field 2
				elif "Ncenter=" in line:
					self.Centers=[MwfnCenter() for i in range(int(value))]
				elif "$Centers" in line:
					for center in self.Centers:
						newline=f.readline()
						newvalues=newline.split()
						center.Symbol=newvalues[1]
						center.Index=int(newvalues[2])
						center.Nuclear_charge=float(newvalues[3])
						center.Coordinates=np.array(newvalues[4:7],dtype="float")*__angstrom2bohr__

				# Field 3
				elif "Nbasis=" in line:
					nbasis=int(value)
					if nindbasis!=114:
						self.Orbitals=[MwfnOrbital() for i in range(nindbasis*(1 if self.Wfntype==0 else 2))]
						for orbital in self.Orbitals:
							orbital.Coeff=np.zeros(nbasis)
				elif "Nindbasis=" in line:
					nindbasis=int(value)
					if nbasis!=114:
						self.Orbitals=[MwfnOrbital() for i in range(nindbasis*(1 if self.Wfntype==0 else 2))]
						for orbital in self.Orbitals:
							orbital.Coeff=np.zeros(nbasis)
				elif "Nshell=" in line:
					self.Shells=[MwfnShell() for i in range(int(value))]
				elif "$Shell types" in line:
					ishell=0
					finished=False
					while not finished:
						newline=f.readline()
						if not newline:
							break
						newvalues=newline.split()
						for newvalue in newvalues:
							self.Shells[ishell].Type=int(newvalue)
							ishell+=1
							if ishell==len(self.Shells):
								finished=True
								break
				elif "$Shell centers" in line:
					ishell=0
					finished=False
					while not finished:
						newline=f.readline()
						if not newline:
							break
						newvalues=newline.split()
						for newvalue in newvalues:
							self.Shells[ishell].Center=self.Centers[int(newvalue)-1]
							ishell+=1
							if ishell==len(self.Shells):
								finished=True
								break
					for center in self.Centers:
						center.Shells=[]
					for shell in self.Shells:
						shell.Center.Shells.append(shell)
				elif "$Shell contraction degrees" in line:
					ishell=0
					finished=False
					while not finished:
						newline=f.readline()
						if not newline:
							break
						newvalues=newline.split()
						for newvalue in newvalues:
							self.Shells[ishell].Exponents=[0 for i in range(int(newvalue))]
							self.Shells[ishell].Coefficients=[0 for i in range(int(newvalue))]
							ishell+=1
							if ishell==len(self.Shells):
								finished=True
								break
				elif "$Primitive exponents" in line:
					ishell=0
					jprim=0
					finished=False
					while not finished:
						newline=f.readline()
						if not newline:
							break
						newvalues=newline.split()
						for newvalue in newvalues:
							self.Shells[ishell].Exponents[jprim]=float(newvalue)
							jprim+=1
							if jprim==len(self.Shells[ishell].Exponents):
								ishell+=1
								jprim=0
							if ishell==len(self.Shells):
								finished=True
								break
				elif "$Contraction coefficients" in line:
					ishell=0
					jprim=0
					finished=False
					while not finished:
						newline=f.readline()
						if not newline:
							break
						newvalues=newline.split()
						for newvalue in newvalues:
							self.Shells[ishell].Coefficients[jprim]=float(newvalue)
							jprim+=1
							if jprim==len(self.Shells[ishell].Coefficients):
								ishell+=1
								jprim=0
							if ishell==len(self.Shells):
								finished=True
								break

				# Field 4
				elif "Index=" in line:
					if mwfntransform is None:
						mwfntransform=self.MatrixTransform()
					tmp_int=int(value)-1
				elif "Type=" in line:
					self.Orbitals[tmp_int].Type=int(value)
				elif "Energy=" in line:
					self.Orbitals[tmp_int].Energy=float(value)
				elif "Occ=" in line:
					self.Orbitals[tmp_int].Occ=float(value)
				elif "Sym=" in line:
					self.Orbitals[tmp_int].Sym=value
				elif "$Coeff" in line:
					ibasis=0
					finished=False
					while not finished:
						newline=f.readline()
						if not newline:
							break
						newvalues=newline.split()
						for newvalue in newvalues:
							try:
								self.Orbitals[tmp_int].Coeff[ibasis]=float(newvalue)
							except ValueError:
								self.Orbitals[tmp_int].Coeff[ibasis]=0.
							ibasis+=1
							if ibasis==self.getNumBasis():
								finished=True
								break

				# Field 5
				elif "$Total density matrix" in line:
					nrows=int(values[4])
					ncols=int(values[5])
					lower=int(values[7])
					self.Total_density_matrix=mwfntransform.T@ReadMatrix(f,nrows,ncols,lower)@mwfntransform
				elif "$Alpha density matrix" in line:
					nrows=int(values[4])
					ncols=int(values[5])
					lower=int(values[7])
					self.Alpha_density_matrix=mwfntransform.T@ReadMatrix(f,nrows,ncols,lower)@mwfntransform
				elif "$Beta density matrix" in line:
					nrows=int(values[4])
					ncols=int(values[5])
					lower=int(values[7])
					self.Beta_density_matrix=mwfntransform.T@ReadMatrix(f,nrows,ncols,lower)@mwfntransform
				elif "$1-e Hamiltonian matrix" in line:
					nrows=int(values[4])
					ncols=int(values[5])
					lower=int(values[7])
					self.Hamiltonian_matrix=mwfntransform.T@ReadMatrix(f,nrows,ncols,lower)@mwfntransform
				elif "$Alpha 1-e Hamiltonian matrix" in line:
					nrows=int(values[4])
					ncols=int(values[5])
					lower=int(values[7])
					self.Alpha_Hamiltonian_matrix=mwfntransform.T@ReadMatrix(f,nrows,ncols,lower)@mwfntransform
				elif "$Beta 1-e Hamiltonian matrix" in line:
					nrows=int(values[4])
					ncols=int(values[5])
					lower=int(values[7])
					self.Beta_Hamiltonian_matrix=mwfntransform.T@ReadMatrix(f,nrows,ncols,lower)@mwfntransform
				elif "$Overlap matrix" in line:
					nrows=int(values[3])
					ncols=int(values[4])
					lower=int(values[6])
					self.Overlap_matrix=mwfntransform.T@ReadMatrix(f,nrows,ncols,lower)@mwfntransform
				elif "$Kinetic energy matrix" in line:
					nrows=int(values[4])
					ncols=int(values[5])
					lower=int(values[7])
					self.Kinetic_energy_matrix=mwfntransform.T@ReadMatrix(f,nrows,ncols,lower)@mwfntransform
				elif "$Potential energy matrix" in line:
					nrows=int(values[4])
					ncols=int(values[5])
					lower=int(values[7])
					self.Potential_energy_matrix=mwfntransform.T@ReadMatrix(f,nrows,ncols,lower)@mwfntransform
				elif "$X electric dipole moment matrix" in line:
					nrows=int(values[6])
					ncols=int(values[7])
					lower=int(values[9])
					self.X_electric_dipole_moment_matrix=mwfntransform.T@ReadMatrix(f,nrows,ncols,lower)@mwfntransform
				elif "$Y electric dipole moment matrix" in line:
					nrows=int(values[6])
					ncols=int(values[7])
					lower=int(values[9])
					self.Y_electric_dipole_moment_matrix=mwfntransform.T@ReadMatrix(f,nrows,ncols,lower)@mwfntransform
				elif "$Z electric dipole moment matrix" in line:
					nrows=int(values[6])
					ncols=int(values[7])
					lower=int(values[9])
					self.Z_electric_dipole_moment_matrix=mwfntransform.T@ReadMatrix(f,nrows,ncols,lower)@mwfntransform
				elif "$XX electric quadrupole moment matrix" in line:
					nrows=int(values[6])
					ncols=int(values[7])
					lower=int(values[9])
					self.XX_electric_quadrupole_moment_matrix=mwfntransform.T@ReadMatrix(f,nrows,ncols,lower)@mwfntransform
				elif "$YY electric quadrupole moment matrix" in line:
					nrows=int(values[6])
					ncols=int(values[7])
					lower=int(values[9])
					self.YY_electric_quadrupole_moment_matrix=mwfntransform.T@ReadMatrix(f,nrows,ncols,lower)@mwfntransform
				elif "$ZZ electric quadrupole moment matrix" in line:
					nrows=int(values[6])
					ncols=int(values[7])
					lower=int(values[9])
					self.ZZ_electric_quadrupole_moment_matrix=mwfntransform.T@ReadMatrix(f,nrows,ncols,lower)@mwfntransform
				elif "$XY electric quadrupole moment matrix" in line:
					nrows=int(values[6])
					ncols=int(values[7])
					lower=int(values[9])
					self.XY_electric_quadrupole_moment_matrix=mwfntransform.T@ReadMatrix(f,nrows,ncols,lower)@mwfntransform
				elif "$YZ electric quadrupole moment matrix" in line:
					nrows=int(values[6])
					ncols=int(values[7])
					lower=int(values[9])
					self.YZ_electric_quadrupole_moment_matrix=mwfntransform.T@ReadMatrix(f,nrows,ncols,lower)@mwfntransform
				elif "$XZ electric quadrupole moment matrix" in line:
					nrows=int(values[6])
					ncols=int(values[7])
					lower=int(values[9])
					self.XZ_electric_quadrupole_moment_matrix=mwfntransform.T@ReadMatrix(f,nrows,ncols,lower)@mwfntransform
		for spin in ([0] if self.Wfntype==0 else [1,2]):
			self.setCoefficientMatrix(spin,mwfntransform.T@self.getCoefficientMatrix(spin))
		Extra_info={}

	def getNumElec(self,spin):
		naelec=0.
		nbelec=0.
		for orbital in self.Orbitals:
			if orbital.Type==0:
				naelec+=orbital.Occ/2
				nbelec+=orbital.Occ/2
			elif orbital.Type==1:
				naelec+=orbital.Occ
			elif orbital.Type==2:
				nbelec+=orbital.Occ
		if spin==0:
			return naelec+nbelec
		elif spin==1:
			return naelec
		elif spin==2:
			return nbelec

	def getCharge(self):
		n=-self.getNumElec(0)
		for center in self.Centers:
			n+=center.Nuclear_charge
		return n

	def getSpin(self):
		return self.getNumElec(1)-self.getNumElec(2)

	def getNumCenters(self):
		return len(self.Centers)

	def getNumBasis(self):
		n=0
		for shell in self.Shells:
			n+=shell.getSize()
		return n

	def getNumIndBasis(self):
		return len(self.Orbitals)//(2 if self.Wfntype==1 else 1)

	def getNumPrims(self):
		n=0
		for shell in self.Shells:
			n+=(abs(shell.Type)+1)*(abs(shell.Type)+2)//2*len(shell.Exponents)
		return n

	def getNumShells(self):
		return len(self.Shells)

	def getNumPrimShells(self):
		n=0
		for shell in self.Shells:
			n+=len(shell.Exponents)
		return n

	def getShellIndexByCenter(self,centers=-1):
		indices=[]
		if centers==-1:
			for center in self.Centers:
				indices.append(self.getShellIndexByCenter(center))
			return indices
		elif type(centers) is list:
			for center in centers:
				indices.append(self.getShellIndexByCenter(center))
			return indices
		else:
			assert type(centers) is MwfnCenter or int
			thiscenter=None
			if type(centers) is MwfnCenter:
				assert centers in self.Centers
				thiscenter=centers
			else:
				thiscenter=self.Centers[centers]
			ishell=0
			for shell in self.Shells:
				if shell.Center is thiscenter:
					indices.append(ishell)
				ishell+=1
		return indices

	def getShellCenters(self):
		center_list=[]
		for shell in self.Shells:
			center_list.append(self.Centers.index(shell.Center))
		return center_list

	def getBasisIndexByCenter(self,centers=-1):
		indices=[]
		if centers==-1:
			for center in self.Centers:
				indices.append(self.getBasisIndexByCenter(center))
			return indices
		elif type(centers) is list:
			for center in centers:
				indices.append(self.getBasisIndexByCenter(center))
			return indices
		else:
			assert type(centers) is MwfnCenter or int
			thiscenter=None
			if type(centers) is MwfnCenter:
				assert centers in self.Centers
				thiscenter=centers
			else:
				thiscenter=self.Centers[centers]
			ibasis=0
			for shell in self.Shells:
				if shell.Center is thiscenter:
					for count in range(shell.getSize()):
						indices.append(ibasis)
						ibasis+=1
				else:
					ibasis+=shell.getSize()
		return indices

	def getBasisIndexByShell(self,shells=-1):
		indices=[]
		if shells==-1:
			for shell in self.Shells:
				indices.append(self.getBasisIndexByShell(shell))
			return indices
		elif type(shells) is list:
			for shell in shells:
				indices.append(self.getBasisIndexByShell(shell))
			return indices
		else:
			assert type(shells) is MwfnShell or int
			thisshell=None
			if type(shells) is MwfnShell:
				assert shells in self.Shells
				thisshell=shells
			else:
				thisshell=self.Shells[shells]
			ibasis=0
			for shell in self.Shells:
				for i in range(shell.getSize()):
					if thisshell==shell:
						indices.append(ibasis)
					ibasis+=1
			return indices

	def getCoefficientMatrix(self,type_):
		matrix=np.zeros([self.getNumBasis(),self.getNumIndBasis()])
		orbital_indices=[i for i in range(len(self.Orbitals)) if self.Orbitals[i].Type==type_]
		for i in range(len(orbital_indices)):
			matrix[:,i]=self.Orbitals[orbital_indices[i]].Coeff
		return matrix

	def setCoefficientMatrix(self,type_,matrix):
		orbital_indices=[i for i in range(len(self.Orbitals)) if self.Orbitals[i].Type==type_]
		assert matrix.shape==(self.getNumBasis(),len(orbital_indices))
		for i in range(len(orbital_indices)):
			self.Orbitals[orbital_indices[i]].Coeff=matrix[:,i]

	def getEnergy(self,type_):
		return [orbital.Energy for orbital in self.Orbitals if orbital.Type==type_]

	def setEnergy(self,type_,energies):
		orbital_indices=[i for i in range(len(self.Orbitals)) if self.Orbitals[i].Type==type_]
		assert len(energies)==len(orbital_indices)
		for i in range(len(orbital_indices)):
			self.Orbitals[orbital_indices[i]].Energy=energies[i]

	def calcOverlap(self):
		self.Overlap_matrix = eint.PyscfOverlap(self,self)

	def calcHamiltonian(self):
		S=self.Overlap_matrix
		if self.Wfntype==0:
			E=np.diag(self.getEnergy(0))
			C=self.getCoefficientMatrix(0)
			self.Hamiltonian_matrix=S@C@E@C.T@S
		elif self.Wfntype==1:
			Ea=np.diag(self.getEnergy(1))
			Ca=self.getCoefficientMatrix(1)
			self.Alpha_Hamiltonian_matrix=S@Ca@Ea@Ca.T@S
			Eb=np.diag(self.getEnergy(2))
			Cb=self.getCoefficientMatrix(2)
			self.Beta_Hamiltonian_matrix=S@Cb@Eb@Cb.T@S

	def getFock(self, spin):
		S = self.Overlap_matrix
		E = np.diag(self.getEnergy(spin))
		C = self.getCoefficientMatrix(spin)
		return S @ C @ E @ C.T @ S

	def getOccupation(self,type_):
		return [orbital.Occ for orbital in self.Orbitals if orbital.Type==type_]

	def setOccupation(self,type_,occupations):
		orbital_indices=[i for i in range(len(self.Orbitals)) if self.Orbitals[i].Type==type_]
		assert len(occupations)==len(orbital_indices)
		for i in range(len(orbital_indices)):
			self.Orbitals[orbital_indices[i]].Occ=occupations[i]

	def calcDensity(self):
		S=self.Overlap_matrix
		if self.Wfntype==0:
			N=np.diag(self.getOccupation(0))
			C=self.getCoefficientMatrix(0)
			self.Total_density_matrix=C@N@C.T
		elif self.Wfntype==1:
			Na=np.diag(self.getOccupation(1))
			Ca=self.getCoefficientMatrix(1)
			self.Alpha_density_matrix=Ca@Na@Ca.T
			Nb=np.diag(self.getOccupation(2))
			Cb=self.getCoefficientMatrix(2)
			self.Beta_density_matrix=Cb@Nb@Cb.T

	def GramSchmidt(self,spin,fix):
		S=self.Overlap_matrix
		C=self.getCoefficientMatrix(spin)
		Cfix=C[:,:fix]
		Ccng=C[:,fix:]
		np.random.seed(0)
		Ccng=np.random.rand(*Ccng.shape)
		Cnew=np.hstack([Cfix,Ccng])
		for i in range(fix,Cnew.shape[1]):
			vector=Cnew[:,i]*1
			for j in range(i):
				Cnew[:,i]-=(vector@S@Cnew[:,j])*Cnew[:,j]
			Cnew[:,i]/=np.sqrt(Cnew[:,i]@S@Cnew[:,i])
		self.setCoefficientMatrix(spin,Cnew)
				
	def Export(self,filename):

		def PrintMatrix(f,matrix,lower):
			for irow in range(matrix.shape[0]):
				for jcol in range(irow+1 if lower else matrix.shape[1]):
					f.write(" "+str(matrix[irow,jcol]))
				f.write("\n")

		with open(filename,'w') as f:
			f.write("# Generated by Orbaplaw\n")
			f.write("\n# Comment: "+self.Comment)

			# Field 1
			f.write("\n\n# Overview\n")
			f.write("Wfntype= "+str(self.Wfntype)+"\n")
			f.write("Charge= "+str(self.getCharge())+"\n")
			f.write("Naelec= "+str(self.getNumElec(1))+"\n")
			f.write("Nbelec= "+str(self.getNumElec(2))+"\n")
			f.write("E_tot= "+str(self.E_tot)+"\n")
			f.write("VT_ratio= "+str(self.VT_ratio)+"\n")

			# Field 2
			f.write("\n\n# Atoms\n")
			f.write("Ncenter= "+str(self.getNumCenters())+"\n")
			f.write("$Centers\n")
			for icenter,center in enumerate(self.Centers):
				f.write(str(icenter+1)+" ")
				f.write(str(center.Symbol)+" ")
				f.write(str(center.Index)+" ")
				f.write(str(center.Nuclear_charge)+" ")
				f.write(str(center.Coordinates[0]/__angstrom2bohr__)+" ")
				f.write(str(center.Coordinates[1]/__angstrom2bohr__)+" ")
				f.write(str(center.Coordinates[2]/__angstrom2bohr__)+"\n")

			# Field 3
			f.write("\n\n# Basis set\n")
			f.write("Nbasis= "+str(self.getNumBasis())+"\n")
			f.write("Nindbasis= "+str(self.getNumIndBasis())+"\n")
			f.write("Nprims= "+str(self.getNumPrims())+"\n")
			f.write("Nshell= "+str(self.getNumShells())+"\n")
			f.write("Nprimshell= "+str(self.getNumPrimShells())+"\n")
			f.write("$Shell types")
			thiscenter=None
			for shell in self.Shells:
				if thiscenter is not shell.Center:
					f.write("\n")
					thiscenter=shell.Center
				f.write(" "+str(shell.Type))
			f.write("\n$Shell centers")
			thiscenter=None
			for shell in self.Shells:
				if thiscenter is not shell.Center:
					f.write("\n")
					thiscenter=shell.Center
				f.write(" "+str(self.Centers.index(shell.Center)+1))
			f.write("\n$Shell contraction degrees")
			for shell in self.Shells:
				if thiscenter is not shell.Center:
					f.write("\n")
					thiscenter=shell.Center
				f.write(" "+str(shell.getNumPrims()))
			f.write("\n$Primitive exponents\n")
			for shell in self.Shells:
				for exponent in shell.Exponents:
					f.write(" "+str(exponent))
				f.write("\n")
			f.write("$Contraction coefficients\n")
			for shell in self.Shells:
				for coefficient in shell.Coefficients:
					f.write(" "+str(coefficient))
				f.write("\n")

			# Field 4
			mwfntransform=self.MatrixTransform()
			f.write("\n\n# Orbitals")
			for iorbital,orbital in zip(range(len(self.Orbitals)),self.Orbitals):
				f.write("\nIndex= %9d\n" % (iorbital+1))
				f.write("Type= "+str(orbital.Type)+"\n")
				f.write("Energy= "+str(orbital.Energy)+"\n")
				f.write("Occ= "+str(orbital.Occ)+"\n")
				f.write("Sym= "+str(orbital.Sym)+"\n")
				f.write("$Coeff\n")
				for element in mwfntransform@orbital.Coeff:
					f.write(" "+str(element))
				f.write("\n")

			# Field 5
			f.write("\n\n# Matrices\n")
			if self.Total_density_matrix is not None:
				f.write("$Total density matrix, dim= "+str(self.Total_density_matrix.shape[0])+" "+str(self.Total_density_matrix.shape[1])+" lower= 1\n")
				PrintMatrix(f,mwfntransform@self.Total_density_matrix@mwfntransform.T,True)
			if self.Alpha_density_matrix is not None:
				f.write("$Alpha density matrix, dim= "+str(self.Alpha_density_matrix.shape[0])+" "+str(self.Alpha_density_matrix.shape[1])+" lower= 1\n")
				PrintMatrix(f,mwfntransform@self.Alpha_density_matrix@mwfntransform.T,True)
			if self.Beta_density_matrix is not None:
				f.write("$Beta density matrix, dim= "+str(self.Beta_density_matrix.shape[0])+" "+str(self.Beta_density_matrix.shape[1])+" lower= 1\n")
				PrintMatrix(f,mwfntransform@self.Beta_density_matrix@mwfntransform.T,True)
			if self.Hamiltonian_matrix is not None:
				f.write("$1-e Hamiltonian matrix, dim= "+str(self.Hamiltonian_matrix.shape[0])+" "+str(self.Hamiltonian_matrix.shape[1])+" lower= 1\n")
				PrintMatrix(f,mwfntransform@self.Hamiltonian_matrix@mwfntransform.T,True)
			if self.Alpha_Hamiltonian_matrix is not None:
				f.write("$Alpha 1-e Hamiltonian matrix, dim= "+str(self.Alpha_Hamiltonian_matrix.shape[0])+" "+str(self.Alpha_Hamiltonian_matrix.shape[1])+" lower= 1\n")
				PrintMatrix(f,mwfntransform@self.Alpha_Hamiltonian_matrix@mwfntransform.T,True)
			if self.Beta_Hamiltonian_matrix is not None:
				f.write("$Beta 1-e Hamiltonian matrix, dim= "+str(self.Beta_Hamiltonian_matrix.shape[0])+" "+str(self.Beta_Hamiltonian_matrix.shape[1])+" lower= 1\n")
				PrintMatrix(f,mwfntransform@self.Beta_Hamiltonian_matrix@mwfntransform.T,True)
			if self.Overlap_matrix is not None:
				f.write("$Overlap matrix, dim= "+str(self.Overlap_matrix.shape[0])+" "+str(self.Overlap_matrix.shape[1])+" lower= 1\n")
				PrintMatrix(f,mwfntransform@self.Overlap_matrix@mwfntransform.T,True)
			if self.Kinetic_energy_matrix is not None:
				f.write("$Kinetic energy matrix, dim= "+str(self.Kinetic_energy_matrix.shape[0])+" "+str(self.Kinetic_energy_matrix.shape[1])+" lower= 1\n")
				PrintMatrix(f,mwfntransform@self.Kinetic_energy_matrix@mwfntransform.T,True)
			if self.Potential_energy_matrix is not None:
				f.write("$Potential energy matrix, dim= "+str(self.Potential_energy_matrix.shape[0])+" "+str(self.Potential_energy_matrix.shape[1])+" lower= 1\n")
				PrintMatrix(f,mwfntransform@self.Potential_energy_matrix@mwfntransform.T,True)
