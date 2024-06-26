import numpy as np


class MwfnCenter:
    Symbol=None # String
    Index=None # Int
    Nuclear_charge=None # Float
    Coordinates=None # np.array (1,3)
    Shells=None # MwfnShell
    def getNumShells(self):
        return len(self.Shells)

class MwfnShell:
    Type=None # int
    Center=None # MwfnCenter
    Exponents=None # List
    Coefficients=None # List
    def getSize(self):
        if self.Type>=0:
            return int((self.Type+1)*(self.Type+2)/2)
        else:
            return int(2*abs(self.Type)+1)
    def getNumPrims(self):
        return len(self.Exponents)
    
class MwfnOrbital:
    Type=None
    Energy=None
    Occ=None
    Sym=None
    Coeff=None # np.array (Nshell,1)

class MultiWaveFunction:

    Comment=""

    # Field 1
    Wfntype=None
    Naelec=None
    Nbelec=None
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

    Extra_info={}

    def __init__(self,filename):
        def ReadMatrix(f,nrows,ncols,lower):
            matrix=np.zeros([nrows,ncols])
            assert nrows==ncols
            nelements=int((1+nrows)*nrows/2 if lower else nrows*ncols)
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
            naorbitals=0
            nborbitals=0
            aorbitals=[]
            borbitals=[]
            this_orbital=None
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
                    self.Comment=line[8:-1]
                
                # Field 1
                elif "Wfntype=" in line:
                    self.Wfntype=int(value)
                elif "Naelec=" in line: # Charge is inferred from Naelec, Nbelec and nuclear charges of centers.
                    self.Naelec=int(float(value))
                elif "Nbelec=" in line:
                    self.Nbelec=int(float(value))
                elif "E_tot=" in line:
                    self.E_tot=float(value)
                elif "VT_ratio=" in line:
                    self.VT_ratio=float(value)

                # Field 2
                elif "Ncenter=" in line:
                    self.Centers=[MwfnCenter() for i in range(int(value))]
                elif "$Centers" in line:
                    icenter=0
                    while True:
                        newline=f.readline()
                        if not newline:
                            break
                        newvalues=newline.split()
                        self.Centers[icenter].Symbol=newvalues[1]
                        self.Centers[icenter].Index=int(newvalues[2])
                        self.Centers[icenter].Nuclear_charge=float(newvalues[3])
                        self.Centers[icenter].Coordinates=np.array(newvalues[4:7],dtype="float")
                        icenter+=1
                        if icenter==len(self.Centers):
                            break

                # Field 3
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
                    if len(aorbitals+borbitals)==0:
                        aorbitals=[MwfnOrbital() for i in range(self.getNumBasis())]
                        for orbital in aorbitals:
                            orbital.Type=0 if self.Wfntype==0 else 1
                            orbital.Energy=0
                            orbital.Occ=-1
                            orbital.Coeff=np.zeros(self.getNumBasis())
                        if self.Wfntype==1:
                            borbitals=[MwfnOrbital() for i in range(self.getNumBasis())]
                            for orbital in borbitals:
                                orbital.Type=2
                                orbital.Energy=0
                                orbital.Occ=-1
                                orbital.Coeff=np.zeros(self.getNumBasis())
                elif "Type=" in line:
                    if int(value)!=2:
                        this_orbital=aorbitals[naorbitals]
                        naorbitals+=1
                    else:
                        this_orbital=borbitals[nborbitals]
                        nborbitals+=1
                elif "Energy=" in line:
                    this_orbital.Energy=float(value)
                elif "Occ=" in line:
                    this_orbital.Occ=float(value)
                elif "Sym=" in line:
                    this_orbital.Sym=value
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
                                this_orbital.Coeff[ibasis]=float(newvalue)
                            except ValueError:
                                this_orbitalCoeff[ibasis]=0.
                            ibasis+=1
                            if ibasis==self.getNumBasis():
                                finished=True
                                break

                # Field 5
                elif "$Total density matrix" in line:
                    nrows=int(values[4])
                    ncols=int(values[5])
                    lower=int(values[7])
                    self.Total_density_matrix=ReadMatrix(f,nrows,ncols,lower)
                elif "$Alpha density matrix" in line:
                    nrows=int(values[4])
                    ncols=int(values[5])
                    lower=int(values[7])
                    self.Alpha_density_matrix=ReadMatrix(f,nrows,ncols,lower)
                elif "$Beta density matrix" in line:
                    nrows=int(values[4])
                    ncols=int(values[5])
                    lower=int(values[7])
                    self.Beta_density_matrix=ReadMatrix(f,nrows,ncols,lower)
                elif "$1-e Hamiltonian matrix" in line:
                    nrows=int(values[4])
                    ncols=int(values[5])
                    lower=int(values[7])
                    self.Hamiltonian_matrix=ReadMatrix(f,nrows,ncols,lower)
                elif "$Alpha 1-e Hamiltonian matrix" in line:
                    nrows=int(values[4])
                    ncols=int(values[5])
                    lower=int(values[7])
                    self.Alpha_Hamiltonian_matrix=ReadMatrix(f,nrows,ncols,lower)
                elif "$Beta 1-e Hamiltonian matrix" in line:
                    nrows=int(values[4])
                    ncols=int(values[5])
                    lower=int(values[7])
                    self.Beta_Hamiltonian_matrix=ReadMatrix(f,nrows,ncols,lower)
                elif "$Overlap matrix" in line:
                    nrows=int(values[3])
                    ncols=int(values[4])
                    lower=int(values[6])
                    self.Overlap_matrix=ReadMatrix(f,nrows,ncols,lower)
                elif "$Kinetic energy matrix" in line:
                    nrows=int(values[4])
                    ncols=int(values[5])
                    lower=int(values[7])
                    self.Kinetic_energy_matrix=ReadMatrix(f,nrows,ncols,lower)
                elif "$Potential energy matrix" in line:
                    nrows=int(values[4])
                    ncols=int(values[5])
                    lower=int(values[7])
                    self.Potential_energy_matrix=ReadMatrix(f,nrows,ncols,lower)
        Extra_info={}
        self.Orbitals=[MwfnOrbital() for i in range(self.getNumBasis()*(2 if self.Wfntype==1 else 1))]
        for iaorbital in range(self.getNumBasis()):
            self.Orbitals[iaorbital]=aorbitals[iaorbital]
            self.Orbitals[iaorbital].Index=iaorbital+1
        if self.Wfntype==1:
            for iborbital in range(self.getNumBasis()):
                self.Orbitals[iborbital+self.getNumBasis()]=borbitals[iborbital]
                self.Orbitals[iborbital+self.getNumBasis()].Index=iborbital+self.getNumBasis()+1

    def getCharge(self):
        neutral=0
        for center in self.Centers:
            neutral+=center.Nuclear_charge
        actual=self.Naelec+self.Nbelec
        return actual-neutral

    def getNumCenters(self):
        return len(self.Centers)

    def getNumBasis(self):
        n=0
        for shell in self.Shells:
            n+=shell.getSize()
        return int(n)

    def getNumIndBasis(self):
        numindbasis=0
        for orbital in self.Orbitals:
            if orbital.Occ!=-1:
                numindbasis+=1
        numindbasis//=(2 if self.Wfntype==1 else 1)
        return numindbasis

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
        matrix=np.zeros([self.getNumBasis(),self.getNumBasis()])
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

    def getOccupation(self,type_):
        return [orbital.Occ for orbital in self.Orbitals if orbital.Type==type_]

    def setOccupation(self,type_,occupations):
        orbital_indices=[i for i in range(len(self.Orbitals)) if self.Orbitals[i].Type==type_]
        assert len(occupations)==len(orbital_indices)
        for i in range(len(orbital_indices)):
            self.Orbitals[orbital_indices[i]].Occ=occupations[i]

    def Export(self,filename):

        def PrintMatrix(f,matrix,lower):
            for irow in range(matrix.shape[0]):
                for jcol in range(irow+1 if lower else matrix.shape[1]):
                    f.write(" "+str(matrix[irow,jcol]))
                f.write("\n")

        with open(filename,'w') as f:
            f.write("# Generated by Orbaplaw\n")
            f.write("\n# Comment: "+self.Comment+"\n")

            # Field 1
            f.write("\n\n# Overview\n")
            f.write("Wfntype= "+str(self.Wfntype)+"\n")
            f.write("Charge= "+str(self.getCharge())+"\n")
            f.write("Naelec= "+str(self.Naelec)+"\n")
            f.write("Nbelec= "+str(self.Nbelec)+"\n")
            f.write("E_tot= "+str(self.E_tot)+"\n")
            f.write("VT_ratio= "+str(self.VT_ratio)+"\n")

            # Field 2
            f.write("\n\n# Atoms\n")
            f.write("Ncenter= "+str(self.getNumCenters())+"\n")
            f.write("$Centers\n")
            icenter=0
            for icenter,center in zip(range(self.getNumCenters()),self.Centers):
                f.write(str(icenter+1)+" ")
                f.write(str(center.Symbol)+" ")
                f.write(str(center.Index)+" ")
                f.write(str(center.Nuclear_charge)+" ")
                f.write(str(center.Coordinates[0])+" ")
                f.write(str(center.Coordinates[1])+" ")
                f.write(str(center.Coordinates[2])+"\n")

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
            f.write("\n\n# Orbitals")
            for iorbital,orbital in zip(range(len(self.Orbitals)),self.Orbitals):
                f.write("\nIndex= %9d\n" % (iorbital+1))
                f.write("Type= "+str(orbital.Type)+"\n")
                f.write("Energy= "+str(orbital.Energy)+"\n")
                f.write("Occ= "+str(orbital.Occ)+"\n")
                f.write("Sym= "+str(orbital.Sym)+"\n")
                f.write("$Coeff\n")
                for element in orbital.Coeff:
                    f.write(" "+str(element))
                f.write("\n")

            # Field 5
            f.write("\n\n# Matrices\n")
            if self.Total_density_matrix is not None:
                f.write("$Total density matrix, dim= "+str(self.Total_density_matrix.shape[0])+" "+str(self.Total_density_matrix.shape[1])+" lower= 1\n")
                PrintMatrix(f,self.Total_density_matrix,True)
            if self.Alpha_density_matrix is not None:
                f.write("$Alpha density matrix, dim= "+str(self.Alpha_density_matrix.shape[0])+" "+str(self.Alpha_density_matrix.shape[1])+" lower= 1\n")
                PrintMatrix(f,self.Alpha_density_matrix,True)
            if self.Beta_density_matrix is not None:
                f.write("$Beta density matrix, dim= "+str(self.Beta_density_matrix.shape[0])+" "+str(self.Beta_density_matrix.shape[1])+" lower= 1\n")
                PrintMatrix(f,self.Beta_density_matrix,True)
            if self.Hamiltonian_matrix is not None:
                f.write("$1-e Hamiltonian matrix, dim= "+str(self.Hamiltonian_matrix.shape[0])+" "+str(self.Hamiltonian_matrix.shape[1])+" lower= 1\n")
                PrintMatrix(f,self.Hamiltonian_matrix,True)
            if self.Alpha_Hamiltonian_matrix is not None:
                f.write("$Alpha 1-e Hamiltonian matrix, dim= "+str(self.Alpha_Hamiltonian_matrix.shape[0])+" "+str(self.Alpha_Hamiltonian_matrix.shape[1])+" lower= 1\n")
                PrintMatrix(f,self.Alpha_Hamiltonian_matrix,True)
            if self.Beta_Hamiltonian_matrix is not None:
                f.write("$Beta 1-e Hamiltonian matrix, dim= "+str(self.Beta_Hamiltonian_matrix.shape[0])+" "+str(self.Beta_Hamiltonian_matrix.shape[1])+" lower= 1\n")
                PrintMatrix(f,self.Beta_Hamiltonian_matrix,True)
            if self.Overlap_matrix is not None:
                f.write("$Overlap matrix, dim= "+str(self.Overlap_matrix.shape[0])+" "+str(self.Overlap_matrix.shape[1])+" lower= 1\n")
                PrintMatrix(f,self.Overlap_matrix,True)
            if self.Kinetic_energy_matrix is not None:
                f.write("$Kinetic energy matrix, dim= "+str(self.Kinetic_energy_matrix.shape[0])+" "+str(self.Kinetic_energy_matrix.shape[1])+" lower= 1\n")
                PrintMatrix(f,self.Kinetic_energy_matrix,True)
            if self.Potential_energy_matrix is not None:
                f.write("$Potential energy matrix, dim= "+str(self.Potential_energy_matrix.shape[0])+" "+str(self.Potential_energy_matrix.shape[1])+" lower= 1\n")
                PrintMatrix(f,self.Potential_energy_matrix,True)
