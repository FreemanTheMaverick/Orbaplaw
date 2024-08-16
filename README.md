# Orbaplaw
> **Orb**ital **a**lignment analysis for **pla**ne **w**ave basis sets
>
> But now we are at the interlude of gaussian basis sets.

## Functions

`Orbaplaw` can be used to perform

+ Population analysis
  + Lowin Population
+ Orbital localization
  + Pipek-Mezey localization
  + Foster-Boys localization
  + Fock space localization (Inefficiently implemented)
  + Localized orbitalet (Inefficiently implemented)
+ Inter-fragment bonding analysis
  + Principal interacting orbital (PIO) analysis
  + Natual fragment bond orbital (NFBO) analysis
  + Fragment-aligned molecular orbital (FAMO) analysis


## Documents

+ [Installation](doc/INSTALLATION.md)
+ [Wavefunction processing (Gaussian-based)](doc/WFN.md)
+ [Wavefunction processing (Plane-wave-based)](doc/PW.md)
+ [Natural atomic orbital](doc/NAO.md)
+ [Principal interacting orbital analysis](doc/PIO.md)
+ [Natural fragment bond orbital analysis](doc/NFBO.md)
+ [Fragment-aligned molecular orbital analysis](doc/FAMO.md)
+ [Orbital localization](doc/LMO.md)


## Citation
Zhang, Y. Orbaplaw: Orbital alignment analysis for plane wave basis sets. https://github.com/FreemanTheMaverick/Orbaplaw, 2024.

If `Orbaplaw` benefits your research, please cite this program and the whole method toolchain you use in your manuscript. For example, if you have used the NFBO method, you had better cite (1) the program `Orbaplaw`, (2) the original paper on NFBO and (3) the original paper on NAO (because NAO is a prerequisite for NFBO and thus part of the toolchain).

![](doc/please_cite.gif)
