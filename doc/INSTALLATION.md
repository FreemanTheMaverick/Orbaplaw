# Installation

## Prerequisites
+ **Python 3**
+ **pip**
+ A C++ compiler that supports C++20 standard (if any of the dependencies need compiling)

## Dependencies (Automatically resolved by **pip**)
+ **NumPy** and **SciPy**
+ **libmwfn** for `mwfn` file format handling
+ **Maniverse** for orbital localization on orthogonal manifolds
+ **PySCF** for electron integrals

## Related tools
+ `Multiwfn` for wavefunction file format conversion and visualization of orbitals

## Download and build
```shell
$ pip install Orbaplaw
```

## Tryout
+ Command-line tool
```shell
$ orbaplaw -h
usage: orbaplaw [-h] {pop,loc,famo,sno,nao,pio,nbo} ...

Command-line tool for Orbaplaw

options:
  -h, --help            show this help message and exit

Job types:
  {pop,loc,famo,sno,nao,pio,nbo}
    pop                 Population Analysis
    loc                 Orbital localization
    famo                Fragment aligned molecular orbital
    sno                 Spin natural orbital
    nao                 Natural atomic orbital
    pio                 Principal interacting orbital
    nbo                 Natural (fragment) bond orbital
```
+ Script
```python
>>> from Orbaplaw import Population as pop
```
