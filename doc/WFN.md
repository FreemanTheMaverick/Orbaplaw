# Wavefunction processing (Gaussian-based)

Here is a typical procedure to deal with wavefunctions in `Orbaplaw`.

## *Ab initio* quantum chemistry computation
A wavefunction file given by an *ab initio* calculation is needed.
This can be done by a variety of computational chemistry packages.
Here we use `Gaussian 16` as an example.
```
%nprocshared=40
%mem=60GB
%chk=job.chk
# b3lyp 6-31g(d) 5d 7f

Title Card Required

0 1
[Geometry]
```

## Converting wavefunction file
The resultant wavefunction is stored in `job.chk`, an `chk` format file which is not supported by `Orbaplaw`.
We need to transform `job.chk` to `fchk` with `formchk` and then to `mwfn` format[^mwfn] with `Multiwfn`[^multiwfn].

```
$ formchk job.chk # Now we have job.fchk
$ Multiwfn job.fchk
100 # Other functions (Part 1)
2 # Export various files (mwfn/pdb/xyz/wfn/wfx/molden/fch/47/mkl...) or generate input file of quantum chemistry programs
32 # Output current wavefunction as .mwfn file
# Default name: job.mwfn
2 # Export wavefunction, density matrix and overlap matrix
0 # Return
q # Exit Multiwfn.
```
Now we have `job.mwfn`.

## Wavefunction processing

+ Loading the `numpy` package and the `WaveFunction` module.
```
import numpy as np
from Orbaplaw import WaveFunction as wfn
```

+ Reading and exporting wavefunction information from `job.mwfn` to `job_new.mwfn`.
```
job_mwfn=wfn.MultiWaveFunction("job.mwfn")
job_mwfn.Export("job_new.mwfn")
```

+ Getting the numbers of electrons.
```
n=mwfn.getNumElec(0) # The total number of electrons.
n_alpha=mwfn.getNumElec(1) # The number of alpha electrons.
n_beta=mwfn.getNumElec(2) # The number of beta electrons.
charge=mwfn.getCharge()
spin=mwfn.getSpin() # equals n_alpha - n_beta.
```
Note that the returned values are all floats. You may convert them to integers with the `round()` function.

+ Getting and setting the occupation numbers.
```
homo=round(mwfn.getNumElec(0)/2)-1 # The indices of the HOMO and LUMO.
lumo=round(mwfn.getNumElec(0)/2)
N=mwfn.getOccupation(0) # 0 for spin-restricted, 1 and 2 for alpha and beta in spin-unrestricted.
N[homo],N[lumo]=N[lumo],N[homo] # Exchanging the occupation numbers of the HOMO and LUMO.
mwfn.setOccupation(0,N)
# or equivalently
N=[orbital.Occ for orbital in mwfn.Orbitals if orbital.Type==0]
mwfn.Orbitals[homo].Occ,mwfn.Orbitals[lumo].Occ=mwfn.Orbitals[lumo].Occ,mwfn.Orbitals[homo].Occ
```

+ Getting and setting the orbital energies.
```
E=mwfn.getEnergy(0) # 0 for spin-restricted, 1 and 2 for alpha and beta in spin-unrestricted.
mwfn.setEnergy(0,np.zeros_like(E) # Setting the orbital energies to zeros.
# or equivalently
E=[orbital.Energy for orbital in mwfn.Orbitals if orbital.Type==0]
for orbital in mwfn.Orbitals:
	if orbital.Type==0:
		orbital.Energy=0.
```

+ Getting and setting the coefficient matrix.
```
C=mwfn.getCoefficientMatrix(0) # 0 for spin-restricted, 1 and 2 for alpha and beta in spin-unrestricted.
C_0=C[:,0] # The coefficient vector of the first orbital.
C_1=C[:,1] # The coefficient vector of the second orbital.
C[:,0]=(C_0+C_1)/1.414213562 # Mixing the first two orbitals.
C[:,1]=(C_0-C_1)/1.414213562
mwfn.setCoefficientMatrix(0,C)
```

## Caution
+ Currently only pure basis functions (5d, 7f, 9g, etc.) are supported.
+ Ordering of basis functions.
The ordering of basis functions of $l \ge 2$ in `Orbaplaw` is [d-2, d-1, d0, d+1, d+2] and [f-3, f-2, f-1, f0, f+1, f+2, f+3], etc., different from [d0, d+1, d-1, d+2, d-2] and [f0, f+1, f-1, f+2, f-2, f+3, f-3, f+4, f-4], etc., in the `mwfn` file format.
Upon I/O of `mwfn` file, a matrix transformation is applied to accommodate this difference.


[^multiwfn]: Lu, T.; Chen, F. Multiwfn: A multifunctional wavefunction analyzer. *J. Comput. Chem.* **2012**, *33*, 580-592. [link](https://doi.org/10.1002/jcc.22885)
[^mwfn]: Lu, T.; Chen, Q. mwfn: A strict, concise and extensible format for electronic wavefunction storage and exchange. *ChemRxiv.* **2022**. This content is a preprint and has not been peer-reviewed. [link](https://doi.org/10.26434/chemrxiv-2021-lt04f-v6)
