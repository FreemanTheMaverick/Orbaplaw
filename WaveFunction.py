import numpy as np

class QeWaveFunction:
    # Wavefunction info produced by QE.
    # See https://gitlab.com/QEF/q-e/-/wikis/Developers/Format-of-data-files for more information.
    ik=114514
    xk=np.empty(3)
    ispin=1919810
    gamma_only=True
    scalef=1.14514
    ngw=114
    igwx=514
    npol=1919
    nbnd=810
    b=np.empty([3,3])
    mill=np.empty([3,100])
    evc=np.empty([100,100])

    def __init__(self,wfcdat):
        # Initializing a QeWaveFunction object from a wfc*.dat file.
        # See https://mattermodeling.stackexchange.com/questions/9149/how-to-read-qes-wfc-dat-files-with-python for more information.
        with open(wfcdat,'rb') as f:
            f.seek(4)
            self.ik=np.fromfile(f,dtype='int32',count=1)[0]
            self.xk=np.fromfile(f,dtype='float64',count=3)
            self.ispin=np.fromfile(f,dtype='int32',count=1)[0]
            self.gamma_only=bool(np.fromfile(f,dtype='int32',count=1)[0])
            self.scalef=np.fromfile(f,dtype='float64',count=1)[0]
            f.seek(4,1)

            f.seek(4,1)
            self.ngw=np.fromfile(f,dtype='int32',count=1)[0]
            self.igwx=np.fromfile(f,dtype='int32',count=1)[0]
            self.npol=np.fromfile(f,dtype='int32',count=1)[0]
            self.nbnd=np.fromfile(f,dtype='int32',count=1)[0]
            f.seek(4,1)

            f.seek(4,1)
            self.b[0,:]=np.fromfile(f,dtype='float64',count=3)
            self.b[1,:]=np.fromfile(f,dtype='float64',count=3)
            self.b[2,:]=np.fromfile(f,dtype='float64',count=3)
            f.seek(4,1)

            f.seek(4,1)
            self.mill=np.fromfile(f,dtype='int32',count=3*self.igwx)
            self.mill=self.mill.reshape((self.igwx,3)).T
            f.seek(4,1)

            self.evc=np.zeros((self.npol*self.igwx,self.nbnd),dtype="complex128")
            for i in range(self.nbnd):
                f.seek(4,1)
                self.evc[:,i]=np.fromfile(f,dtype='complex128',count=self.npol*self.igwx)
                f.seek(4,1)

    def Export(self,wfcdat):
        # Exporting the QeWaveFunction info to a wfc*.dat file.
        with open(wfcdat,'wb') as f:
            f.write(np.array([44],dtype='int32').tobytes('F')) # Fortran feature: The byte size of the data block is shown at its head and tail.
            f.write(np.array([self.ik],dtype='int32').tobytes('F')) # Using fortran byte format.
            f.write(np.array(self.xk,dtype='float64').tobytes('F')) # Reiterating the dtype is necessary because, say, there are int8, int32 and int64 but QE recognizes int32 only.
            f.write(np.array([self.ispin],dtype='int32').tobytes('F'))
            f.write(np.array([self.gamma_only],dtype='int32').tobytes('F'))
            f.write(np.array([self.scalef],dtype='float64').tobytes('F'))
            f.write(np.array([44],dtype='int32').tobytes('F'))

            f.write(np.array([16],dtype='int32').tobytes('F'))
            f.write(np.array([self.ngw],dtype='int32').tobytes('F'))
            f.write(np.array([self.igwx],dtype='int32').tobytes('F'))
            f.write(np.array([self.npol],dtype='int32').tobytes('F'))
            f.write(np.array([self.nbnd],dtype='int32').tobytes('F'))
            f.write(np.array([16],dtype='int32').tobytes('F'))

            f.write(np.array([self.b.size*8],dtype='int32').tobytes('F'))
            f.write(np.array(self.b[0,:],dtype='float64').tobytes('F'))
            f.write(np.array(self.b[1,:],dtype='float64').tobytes('F'))
            f.write(np.array(self.b[2,:],dtype='float64').tobytes('F'))
            f.write(np.array([self.b.size*8],dtype='int32').tobytes('F'))

            f.write(np.array([self.mill.size*4],dtype='int32').tobytes('F'))
            f.write(np.array(self.mill,dtype='int32').tobytes('F'))
            f.write(np.array([self.mill.size*4],dtype='int32').tobytes('F'))

            for i in range(self.nbnd):
                f.write(np.array([self.evc[:,i].size*16],dtype='int32').tobytes('F'))
                f.write(np.array(self.evc[:,i],dtype='complex128').tobytes('F'))
                f.write(np.array([self.evc[:,i].size*16],dtype='int32').tobytes('F'))



def CompareQeWF(wf1,wf2): # Comparing the two QeWaveFunction objects.
    print('ik',wf1.ik-wf2.ik)
    print('xk',wf1.xk-wf2.xk)
    print('ispin',wf1.ispin-wf2.ispin)
    print('gamma_only',wf1.gamma_only is wf2.gamma_only)
    print('scalef',wf1.scalef-wf2.scalef)
    print('ngw',wf1.ngw-wf2.ngw)
    print('igwx',wf1.igwx-wf2.igwx)
    print('npol',wf1.npol-wf2.npol)
    print('nbnd',wf1.nbnd-wf2.nbnd)
    print('b',wf1.b-wf2.b)
    print('mill',np.linalg.norm(wf1.mill-wf2.mill))
    print('evc',np.linalg.norm(wf1.evc-wf2.evc))



if __name__=='__main__':
    # Checking whether parsing and exporting is correct.
    from sys import argv
    wf1=QeWaveFunction(argv[1]) # Creating a QeWaveFunction object from a wfc*.dat file.
    wf1.Export(argv[2]) # Exporting it to another wfc*.dat file.
    wf2=QeWaveFunction(argv[2]) # Reading the two wfc*.dat file and checking whether they are the same.
    CompareQeWF(wf1,wf2)
    


