import numpy as np
import xml.etree.ElementTree as et
import os
import shutil as st
import copy as cp

class QeKPoint:
    x=None
    e=None
    npws=None
    mill=None
    ngw=None
    igwx=None
    evc=None

class QeWaveFunction:
    # Wavefunction info produced by QE.
    # For more information, see
    # https://mattermodeling.stackexchange.com/questions/9149/how-to-read-qes-wfc-dat-files-with-python
    # and 
    # https://gitlab.com/QEF/q-e/-/wikis/Developers/Format-of-data-files
    # .

    xml=None
    chargedat=None
    ispin=None
    gamma_only=None
    scalef=None
    npol=None
    nbnd=None
    b=None
    evc=None
    kpoints=None

    def __init__(self,folder):
        # Initializing a QeWaveFunction object from a *.save folder.

        # data-file-schema.xml
        self.xml=et.parse(folder+"/data-file-schema.xml")
        self.kpoints=[]
        for i in self.xml.getroot().find("output").find("band_structure").findall("ks_energies"):
            self.kpoints.append(QeKPoint())
            self.kpoints[-1].npws=np.array(i.find("npw").text.split(),'int32')[0]
            self.kpoints[-1].e=np.array(i.find("eigenvalues").text.split(),'float64')

        # charge-density.dat
        with open(folder+"/charge-density.dat","rb") as f:
            self.chargedat=f.read()

        # wfc*.dat
        for ik,k in zip(range(len(self.kpoints)),self.kpoints):
            with open(folder+"/wfc"+str(ik+1)+".dat",'rb') as f:

                f.seek(4)
                np.fromfile(f,dtype='int32',count=1)[0]
                k.x=np.fromfile(f,dtype='float64',count=3)
                self.ispin=np.fromfile(f,dtype='int32',count=1)[0]
                self.gamma_only=bool(np.fromfile(f,dtype='int32',count=1)[0])
                self.scalef=np.fromfile(f,dtype='float64',count=1)[0]
                f.seek(4,1)

                f.seek(4,1)
                k.ngw=np.fromfile(f,dtype='int32',count=1)[0]
                k.igwx=np.fromfile(f,dtype='int32',count=1)[0]
                self.npol=np.fromfile(f,dtype='int32',count=1)[0]
                self.nbnd=np.fromfile(f,dtype='int32',count=1)[0]
                f.seek(4,1)

                f.seek(4,1)
                self.b=np.zeros([3,3])
                self.b[0,:]=np.fromfile(f,dtype='float64',count=3)
                self.b[1,:]=np.fromfile(f,dtype='float64',count=3)
                self.b[2,:]=np.fromfile(f,dtype='float64',count=3)
                f.seek(4,1)

                f.seek(4,1)
                k.mill=np.fromfile(f,dtype='int32',count=3*k.igwx)
                k.mill=k.mill.reshape([k.igwx,3]).T
                f.seek(4,1)

                k.evc=np.zeros((self.npol*k.igwx,self.nbnd),dtype="complex128")
                for i in range(self.nbnd):
                    f.seek(4,1)
                    k.evc[:,i]=np.fromfile(f,dtype='complex128',count=self.npol*k.igwx)
                    f.seek(4,1)

    def Export(self,folder):
        if (os.path.exists(folder)):
            st.rmtree(folder)
        os.mkdir(folder)

        # data-file-schema.xml
        newxml=cp.deepcopy(self.xml)
        for i,k in zip(
                newxml.getroot().find("output").find("band_structure").findall("ks_energies"),
                self.kpoints):
            i.find("eigenvalues").text=np.array2string(k.e)[1:-1]
        newxml.write(folder+"/data-file-schema.xml")

        # charge-density.dat
        with open(folder+"/charge-density.dat","wb") as f:
            f.write(self.chargedat)
        '''
        with open(folder+"/charge-density.dat","wb") as f:
            f.write(np.array([12],dtype='int32').tobytes('F'))
            f.write(np.array([self.gamma_only],dtype='int32').tobytes('F'))
            f.write(np.array([self.igwx],dtype='int32').tobytes('F'))
            f.write(np.array([self.ispin],dtype='int32').tobytes('F'))
            f.write(np.array([12],dtype='int32').tobytes('F'))

            f.write(np.array([self.b.size*8],dtype='int32').tobytes('F'))
            f.write(np.array(self.b[0,:],dtype='float64').tobytes('F'))
            f.write(np.array(self.b[1,:],dtype='float64').tobytes('F'))
            f.write(np.array(self.b[2,:],dtype='float64').tobytes('F'))
            f.write(np.array([self.b.size*8],dtype='int32').tobytes('F'))

            f.write(np.array([self.mill.size*4],dtype='int32').tobytes('F'))
            f.write(np.array(self.mill,dtype='int32').tobytes('F'))
            f.write(np.array([self.mill.size*4],dtype='int32').tobytes('F'))

            f.write(np.array([self.igwx*4],dtype='int32').tobytes('F'))
            for i in range(self.igwx):
                f.write(np.array([self.igwx*4],dtype='int32').tobytes('F'))
            f.write(np.array([self.igwx*4],dtype='int32').tobytes('F'))
        '''

        # wfc*.dat
        for ik,k in zip(range(len(self.kpoints)),self.kpoints):
            with open(folder+"/wfc"+str(ik+1)+".dat",'wb') as f:
                f.write(np.array([44],dtype='int32').tobytes('F')) # Fortran feature: The byte size of the data block is shown at its head and tail.
                f.write(np.array([ik+1],dtype='int32').tobytes('F')) # Using fortran byte format.
                f.write(np.array(k.x,dtype='float64').tobytes('F')) # Reiterating the dtype is necessary because, say, there are int8, int32 and int64 but QE recognizes int32 only.
                f.write(np.array([self.ispin],dtype='int32').tobytes('F'))
                f.write(np.array([self.gamma_only],dtype='int32').tobytes('F'))
                f.write(np.array([self.scalef],dtype='float64').tobytes('F'))
                f.write(np.array([44],dtype='int32').tobytes('F'))

                f.write(np.array([16],dtype='int32').tobytes('F'))
                f.write(np.array([k.ngw],dtype='int32').tobytes('F'))
                f.write(np.array([k.igwx],dtype='int32').tobytes('F'))
                f.write(np.array([self.npol],dtype='int32').tobytes('F'))
                f.write(np.array([self.nbnd],dtype='int32').tobytes('F'))
                f.write(np.array([16],dtype='int32').tobytes('F'))

                f.write(np.array([self.b.size*8],dtype='int32').tobytes('F'))
                f.write(np.array(self.b[0,:],dtype='float64').tobytes('F'))
                f.write(np.array(self.b[1,:],dtype='float64').tobytes('F'))
                f.write(np.array(self.b[2,:],dtype='float64').tobytes('F'))
                f.write(np.array([self.b.size*8],dtype='int32').tobytes('F'))

                f.write(np.array([k.mill.size*4],dtype='int32').tobytes('F'))
                f.write(np.array(k.mill,dtype='int32').tobytes('F'))
                f.write(np.array([k.mill.size*4],dtype='int32').tobytes('F'))

                for i in range(self.nbnd):
                    f.write(np.array([k.evc[:,i].size*16],dtype='int32').tobytes('F'))
                    f.write(np.array(k.evc[:,i],dtype='complex128').tobytes('F'))
                    f.write(np.array([k.evc[:,i].size*16],dtype='int32').tobytes('F'))



def CompareQeWF(wf1,wf2): # Comparing the two QeWaveFunction objects.
    print('ispin',wf1.ispin-wf2.ispin)
    print('gamma_only',wf1.gamma_only is wf2.gamma_only)
    print('scalef',wf1.scalef-wf2.scalef)
    print('npol',wf1.npol-wf2.npol)
    print('nbnd',wf1.nbnd-wf2.nbnd)
    print('b',wf1.b-wf2.b)
    x=e=npws=mill=ngw=igwx=evc=0
    for i,j in zip(wf1.kpoints,wf2.kpoints):
        x+=np.linalg.norm(i.x-j.x)
        e+=np.linalg.norm(i.e-j.e)
        npws+=np.linalg.norm(i.npws-j.npws)
        mill+=np.linalg.norm(i.mill-j.mill)
        ngw+=np.linalg.norm(i.ngw-j.ngw)
        igwx+=np.linalg.norm(i.igwx-j.igwx)
        evc+=np.linalg.norm(i.evc-j.evc)
    print("x",x)
    print("e",e)
    print("npws",npws)
    print("mill",mill)
    print("ngw",ngw)
    print("igwx",igwx)
    print("evc",evc)


if __name__=='__main__':
    # Checking whether parsing and exporting is correct.
    from sys import argv
    wf1=QeWaveFunction(argv[1]) # Creating a QeWaveFunction object from a *.save folder.
    wf1.Export(argv[2]) # Exporting it to another *.save folder.
    wf2=QeWaveFunction(argv[2]) # Reading the two wfc*.dat file and checking whether they are the same.
    CompareQeWF(wf1,wf2)
