# Orbaplaw
Orbital alignment analysis for plane wave basis sets

## Functions
`Orbaplaw` can be used to perform
+ Principal interacting orbital (PIO) analysis,
+ Natual fragment bond orbital (NFBO) analysis,
+ Fragment-aligned molecular orbital (FAMO) analysis (under development).

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
The resultant wavefunction is stored in `job.chk`, an `chk` format file which is not supported by `Orbaplaw`.
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
### Making a `Python` script for NFBO analysis
Here is an example `Python` script for NFBO analysis, named `job.py`.
```
from Orbaplaw import WaveFunction as wfn # Orbaplaw.WaveFunction is a subpackage for wavefunction file I/O.
from Orbaplaw import NaturalBondOrbitalMethods as nbo # Orbaplaw.NaturalBondOrbitalMethods is a subpackage for PIO and N(F)BO analysis.

frag1=[i for i in range(43)] # NFBO analysis needs the user to manually divide the molecule into fragments. In this example, we simply divide the molecule into two fragments.
frag2=[i for i in range(43,58)] # The first fragment covers Atoms 1-43 and the second 43-58. Note that indices start from 0 in Python.
mo=wfn.MultiWaveFunction("job.mwfn") # Reading the wavefunction stored in the file "job.mwfn" into a MultiWaveFunction object.
nao=nbo.NaturalAtomicOrbital(mo) # Transforming the wavefunction into NAO basis set.
nao.Export("job_nao.mwfn") # Exporting the NAOs to the file "job_nao.mwfn".
nfho,nfbo=nbo.NaturalBondOrbital(nao,frags=[frag1,frag2]) # Performing NFBO analysis based on the NAO-based wavefunction. Both NFHOs and NFBOs will be returned.
nfho.Export("job_nfho.mwfn") # Exporting the NFHOs and NFBOs into ".mwfn" file.
nfbo.Export("job_nfbo.mwfn")
```
Run the command `$ python job.py` and you will find the NFBO information printed on the screen, including the orbital indices, the population, the coefficients of NFHOs contributing to NFBOs and the fragments they belong to.
```
Fragment combination (0, 1)
NBO_109 (2.0)  =  0.265 * NHO_109 (0.141, F_0)  0.964 * NHO_110 (1.859, F_1)
NBO_110 (0.0)  =  -0.964 * NHO_109 (0.141, F_0)  0.265 * NHO_110 (1.859, F_1)
NBO_111 (2.0)  =  -0.899 * NHO_111 (1.615, F_0)  -0.438 * NHO_112 (0.385, F_1)
NBO_112 (0.0)  =  0.438 * NHO_111 (1.615, F_0)  -0.899 * NHO_112 (0.385, F_1)
```
You can view the orbitals (NAOs, NFHOs and NFBOs) in the ".mwfn" files with Multiwfn.

## Gallery
### Transition state of Diels-Alder reaction between dodecahexaene and ethene
See [^nfbo].
### The delocalized $\sigma$ bonds in [Zn<sub>3</sub>Cp<sub>3</sub>]<sup>+</sup> and [{Th(C<sub>8</sub>H<sub>8</sub>)Cl<sub>2</sub>}<sub>3</sub>]<sup>2-</sup>
See [^zn3cp3] [^th1] [^th2] [^th3] and [^nfbo].
### The B<sub>19</sub><sup>-</sup> cluster: A 2-D three-layer model
See [^wankel1] [^wankel2] and [^nfbo].

[^nfbo]: [This papar](https://doi.org/10.26434/chemrxiv-2024-rt585) elaborates on the concept of NFBO. It is written in a way as pedagogical as possible.
[^zn3cp3]: Freitag, K.; Gemel, C.; Jerabek, P.; Oppel, I. M.; Seidel, R. W.; Frenking, G.; Banh, H.;Dilchert, K.; Fischer, R. A. The σ-aromatic clusters [Zn<sub>3</sub>]<sup>+</sup> and [Zn<sub>2</sub>Cu]: Embryonic brass. *Angew. Chem. Int. Ed.* **2015**, *54*, 4370–4374. [link](https://onlinelibrary.wiley.com/doi/10.1002/anie.201410737)
[^th1]: Boronski, J. T.; Seed, J. A.; Hunger, D.; Woodward, A. W.; van Slageren, J.; Wooles, A. J.; Natrajan, L. S.; Kaltsoyannis, N.; Liddle, S. T. A crystalline tri-thorium cluster with σ-aromatic metal–metal bonding. *Nature* **2021**, *598*, 72–75. [link](https://doi.org/10.1038/s41586-021-03888-3)
[^th2]: Cuyacot, B. J. R.; Foroutan-Nejad, C. [{Th(C<sub>8</sub>H<sub>8</sub>)Cl<sub>2</sub>}<sub>3</sub>]<sup>2-</sup> is stable but not aromatic. *Nature* **2022**, *603*, E18–E20. [link](https://doi.org/10.1038/s41586-021-04319-z)
[^th3]: Boronski, J. T.; Seed, J. A.; Hunger, D.; Woodward, A. W.; van Slageren, J.; Wooles, A. J.; Natrajan, L. S.; Kaltsoyannis, N.; Liddle, S. T. Reply to: [{Th(C<sub>8</sub>H<sub>8</sub>)Cl<sub>2</sub>}<sub>3</sub>]<sup>2-</sup> is stable but not aromatic. *Nature* **2022**, *603*, E21–E22. [link](https://doi.org/10.1038/s41586-021-04320-6)
[^wankel1]: Huang, W.; Sergeeva, A. P.; Zhai, H.-J.; Averkiev, B. B.; Wang, L.-S.; Boldyrev, A. I. A concentric planar doubly π-aromatic B<sub>19</sub><sup>-</sup> cluster. *Nat. Chem.* **2010**, *2*, 202–206. [link](https://doi.org/10.1038/nchem.534)
[^wankel2]: Li, R.; You, X.-R.; Guo, J.-C.; Zhai, H.-J. Concentric inner 2π/6σ and outer 10π/14σ aromaticity underlies the dynamic structural fluxionality of planar B<sub>19</sub><sup>-</sup> Wankel motor cluster. *J. Phys. Chem. A* **2021**, *125*, 5022–5030. [link](https://doi.org/10.1021/acs.jpca.1c02764)
