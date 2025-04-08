# Installation
## Prerequisites
+ `Python 3` (`Anaconda` is recommended)
+ `Numpy`
+ `Scipy`
+ `Pyscf` for electron integrals
+ `Multiwfn` for visualization of orbitals

## Download
+ `$ wget https://github.com/FreemanTheMaverick/Orbaplaw/archive/refs/tags/v1.0.zip`
+ `$ unzip -d [Installation directory] v1.0.zip

## Environment variables
+ `$ vim ~/.bashrc`
+ Input the following commands
  ```
  # Orbaplaw
  export PYTHONPATH=[Installation directory]/Orbaplaw:$PYTHONPATH
  ```
+ `$ source ~/.bashrc`
