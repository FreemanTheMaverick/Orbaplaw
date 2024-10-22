# Natural atomic orbital

## Contents
+ [Theory](#theory)
+ [Usage](#usage)
+ [Gallery](#gallery)


## Theory
[^nao]

## Usage

+ *Ab initio* computation.

In `Orbaplaw`, NAO construction is supported only for **pure** spherical gaussian basis functions, so the keyword `5d 7f` is necessary in the route section of the `Gaussian` input file.

+ Loading necessary packages.
```
from Orbaplaw import WaveFunction as wfn
from Orbaplaw import NaturalBondOrbitalMethods as nbo # Orbaplaw.NaturalBondOrbitalMethods contains the NAO function.
```

+ Loading the wavefunction file.
```
job_mwfn=wfn.MultiWaveFunction("job.mwfn")
```

+ Calculating the density matrix, which is necessary for generation of NAOs.
```
job_mwfn.calcDensity()
```

+ Conducting NAO population analysis.
```
job_nao_mwfn=nbo.NaturalAtomicOrbital(job_mwfn)
```

+ Exporting the NAOs to `job_nao.mwfn`.
```
job_nao_mwfn.Export("job_nao.mwfn")
```

## Gallery


[^nao]: Reed, A. E.; Weinstock, R. B.; Weinhold, F. Natural population analysis. *J. Chem. Phys.* **1985**, *83*, 735–746. [link](https://doi.org/10.1063/1.449486)
