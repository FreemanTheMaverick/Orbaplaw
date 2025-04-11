# Installation

## Prerequisites
+ `Python 3` (`Anaconda` is recommended)
+ `Numpy`
+ `Scipy`
+ `Pyscf` for electron integrals
+ `Maniverse` for orbital localization on orthogonal manifolds
+ `Multiwfn` for visualization of orbitals

## Download - Hard 
+ `$ wget https://github.com/FreemanTheMaverick/Orbaplaw/archive/refs/tags/v1.0.zip`
+ `$ unzip -d [Installation directory] v1.0.zip`
+ Set the environment variable
```
$ export PYTHONPATH=[Installation directory]/Orbaplaw:$PYTHONPATH
```
+ Handle dependencies by yourself.
+ + `pip install numpy, scipy, pyscf`
+ + Install `Maniverse` from [here](https://github.com/FreemanTheMaverick/Maniverse.git) and set its environment variable.

## Download - Easy
+ If you have not installed `Maniverse` via `pip`, set the following environment variable.
```
$ export PYTHON3=[The path where you can find "Python.h".] # You may check this by the command "locate Python.h".
```
+ `pip install Orbaplaw`
Usually `pip` installs packages to a `lib/` directory that is already in `$PYTHONPATH`, so you do not need to set the environment variable for Orbaplaw.
