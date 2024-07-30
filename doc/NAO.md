# Natural atomic orbital

## Contents
+ [Theory](#theory)
+ [Usage](#usage)
+ [Gallery](#gallery)


## Theory

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

+ Conducting NAO population analysis.
```
job_nao_mwfn=nbo.NaturalAtomicOrbital(job_mwfn)
```

+ Exporting the NAO to `job_nao.mwfn`.
```
job_nao_mwfn.Export("job_nao.mwfn")
```

## Gallery