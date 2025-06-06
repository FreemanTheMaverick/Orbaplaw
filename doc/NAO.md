# Natural atomic orbital

## Contents
+ [Theory](#theory)
+ [Usage](#usage)
+ [Gallery](#gallery)


## Theory
[^nao]

## Usage

### *Ab initio* computation

In **Orbaplaw**, NAO construction is supported only for **pure** spherical gaussian basis functions, so the keyword `5d 7f` is necessary in the route section of the **Gaussian** input file.

### Command-line tool

The following commands stores the NAOs generated from `job.mwfn` into a new `job_nao.mwfn`.
```shell
$ orbaplaw nao -h
usage: orbaplaw nao [-h] -i INPUT -o OUTPUT

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Original mwfn file (Required)
  -o OUTPUT, --output OUTPUT
                        Mwfn file for exported NAOs (Required)
$ nao -i job.mwfn -o job_nao.mwfn
```

### Script

You may also write a script to do the same thing if you want to go into details.
+ Loading necessary packages.
```
from Orbaplaw import WaveFunction as wfn
from Orbaplaw import NaturalBondOrbitalMethods as nbo # Orbaplaw.NaturalBondOrbitalMethods contains the NAO function.
```

+ Loading the wavefunction file.
```python
job_mwfn = wfn.MultiWaveFunction("job.mwfn")
```

+ Conducting NAO population analysis.
```python
job_nao_mwfn, job_nao_info = nbo.NaturalAtomicOrbital(job_mwfn)
```

+ Exporting the NAOs to `job_nao.mwfn`.
```python
job_nao_mwfn.Export("job_nao.mwfn")
```

## Gallery


[^nao]: Reed, A. E.; Weinstock, R. B.; Weinhold, F. Natural population analysis. *J. Chem. Phys.* **1985**, *83*, 735–746. [link](https://doi.org/10.1063/1.449486)
