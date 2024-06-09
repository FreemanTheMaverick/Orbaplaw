# Orbaplaw
Orbital alignment analysis for plane wave basis sets

## Functions
`Orbaplaw` can be used to perform
+ Principal interacting orbital (PIO) analysis,
+ Natual fragment bond orbital (NFBO) analysis,
+ Fragment-aligned molecular orbital (FAMO) analysis (in development).

## Installation
### Prerequisites
`Orbaplaw` is written in `Python 3` with several commonly used scientific computation packages.
Therefore, a recent distribution of `Anaconda` is all `Orbaplaw` needs.
### Download
+ `$ cd [Installation directory]`
+ `$ git clone https://github.com/FreemanTheMaverick/Orbaplaw.git`
### Environment variables
+ `$ vim ~/.bashrc`
+ Input the following texts to the end
  ```
  # Orbaplaw
  export PYTHONPATH=[Installation directory]/Orbaplaw:$PYTHONPATH
  ```
+ `$ source ~/.bashrc`
  
## Usage
Here is a typical procedure to perform NFBO analysis with `Orbaplaw`.
### *Ab initio* quantum chemistry computation
A single-configuration wavefunction given by an HF/DFT calculation is needed for NFBO analysis.
This can be done by a variety of computational chemistry packages.
Here we use `Gaussian 16` as an example.
```
%nprocshared=40
%mem=60GB
%chk=job.chk
# b3lyp 6-31g(d) 5d

Title Card Required

0 1
......
```
NFBO analysis requires a basis set of natural atomic orbitals (NAOs).
In `Orbaplaw`, NAO construction is supported only for pure spherical gaussian basis functions, so the keyword `5d` is necessary in the route section.
### Converting wavefunction file
The resultant wavefunction is store in `job.chk`, an `chk` format file which is not supported by `Orbaplaw`.
We need to transform `job.chk` to `fchk` with `formchk` and then to `mwfn` format with `Multiwfn`.
```
$ formchk job.chk # Now we have job.fchk
$ Multiwfn job.fchk
100 # Other functions (Part 1)
2 # Export various files (mwfn/pdb/xyz/wfn/wfx/molden/fch/47/mkl...) or generate input file of quantum chemistry programs
32 # Output current wavefunction as .mwfn file
# Default name: job.mwfn
2 # Export wavefunction, density matrix and overlap matrix
0 # Return
q # Exit Multiwfn. Now we have job.mwfn
```
## Gallery
