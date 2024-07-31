# Fragment aligned molecular orbital

Fragment aligned molecular orbitals is transformed from canonical molecular orbitals of a molecule so that some orbitals resemble the CMOs of the fragment(s) to the largest degree.
It has an auspicious abbreviation, FAMO.
The first syllable "FA" goes like /'fa:/, sounding like "making a fortune" in Chinese, as pointed out by Professor Lin.

## Contents
+ [Theory](#theory)
+ [Usage](#usage)
+ [Gallery](#gallery)


## Theory
Sometimes we want to study the electronic structure of a particular part of a molecule, only to find that the canonical molecular orbitals (CMOs) are too delicalized that the orbitals of the targeted part are contaminated by the neighbours.
How can we eliminate the contamination and focus on the electronic structure of the interesting part?

For this problem, fragment aligned molecular orbital (FAMO) analysis adopts the idea of computing the neighbouring fragment CMOs and projecting them out from the wholistic CMOs.
Then ideally, only the orbitals semilocalized on the interesting subsystems are left.
By doing so, the wholist CMOs are divided into two spaces, the "matched space" matched with the neighbouring fragment CMOs and the "mismatched space" corresponding to the orbitals on the rest of the molecule, whose electronic structure is to be investigated.

Let's formulate the problem mathematically.
The **occupied** CMOs of the whole molecule $\left\{|\psi_i\rangle\right\}$ and the fragments $\left\{|\phi_i\rangle\right\}$ are represented in two basis sets $\left\{|\chi_\mu\rangle\right\}$ and $\left\{|\omega_\nu\rangle\right\}$.
The two basis sets can be the same.

$$
|\psi_i\rangle=\sum_\mu^m A_{\mu i}|\chi_\mu\rangle\\
|\phi_j\rangle=\sum_\nu^p B_{\nu j}|\omega_\nu\rangle
$$

$\mathbf{A}$ and $\mathbf{B}$ are the two coefficient matrices.
$m$ and $p$ are the numbers of basis functions used in *ab inito* computation of the whole molecule and the fragments, respectively.
We want to find a unitary transformation, $\mathbf{X}$, of $|\psi_i\rangle$

$$
\begin{aligned}
|\tilde{\psi}_k\rangle
&=\sum_i^n X_{ik}|\psi_i\rangle\\
&=\sum_i^n\sum_\mu^m A_{\mu i} X_{ik}|\chi_\mu\rangle
\end{aligned}
$$

where $n$ is the number of occupied CMOs in the whole molecule, so that the resulted orbitals $\left\{|\tilde{\psi}_k \rangle\right\}$ have the best matched and mismatched spaces.
Finding the best matched and mismatched spaces boils down to the optimization problem, where the wholistic CMOs are transform so that (1) part of the orbitals have as large overlap with the fragment CMOs as possible to construct the matched space and (2) the other part have as small overlap with the fragment CMOs as possible to construct the mismatched space.
That is to find the minimum of the loss function

$$
\begin{aligned}
    L
    &=\sum_j^q\sum_k^n\left(\langle \phi_j |\tilde{\psi}_k\rangle -\delta_{jk}\right)^2\\
    &=\sum_j^q\sum_k^n\left(\sum_\nu^p B_{\nu j}^* \sum_i^n\sum_\mu^m A_{\mu i} X_{ik} \langle \omega_\nu | \chi_\mu \rangle -\delta_{jk}\right)^2\\
    &=\sum_j^q\sum_k^n\left[\left(\mathbf{B^\dagger S^\dagger AX-I}\right)_{jk}\right]^2\\
    &=\left\| \mathbf{\mathbf{B^\dagger S^\dagger AX-I}} \right\|_\text{F}^2
\end{aligned}
$$

with the unitary constraint on $\mathbf{X}$, where $q$ is the number of occupied CMOs in the fragments, $\mathbf{S}$ is the overlap between the two basis sets and $\mathbf{I}$ is the rectangular identity matrix.

$$
\begin{aligned}
    S_{\mu\nu}&=\langle \chi_\mu | \omega_\nu \rangle\\
    I_{jk}&=\delta_{jk}
\end{aligned}
$$

It is a typical Procrustes problem,

$$
\min_\mathbf{X}\Bigl\{\Vert \mathbf{\mathbf{B^\dagger S^\dagger AX-I}} \Vert_\text{F}^2,\ \mathbf{X}^{-1}=\mathbf{X}^\dagger\Bigr\}
$$

whose solution is

$$
\mathbf{X=UV^\dagger}
$$

by singular value decomposition.

$$
\mathbf{A^\dagger SB I=U\Sigma V^\dagger}
$$

The resulted orbitals are called fragment aligned molecular orbitals (FAMOs).
Their corresponding singular values $\left\{\sigma_f\right\}$ are the overlap between the matched space and the fragment CMOs, ranging in $\left(0,1\right)$.
A set of singular values close to 1 indicates that the fragments are well-chosen.
However, the FAMOs in the mismatched space do not have singular values, because they are meant to have zero overlap with the fragment CMOs.
They mix each other due to degeneracy and have no chemical meaning.
To lift the degeneracy, we project the Fock matrix to the FAMO basis set and diagonalized the mismatched-mismatched block to obtain the final mismatched orbitals.
The final mismatched FAMOs have the CMO features in accordance with general chemical intuition.

The readers may take a look at the original paper[^famo] for more information.

## Usage
+ *Ab initio* computation.

*Ab initio* computation needs to be done on both the molecule and the fragments.
Some quantum chemistry packages, such as `Gaussian`, rotates the input geometry to a standard orientation by default.
In FAMO analysis, users must disable this action by using extra keywords, such as `nosym` in `Gaussian`.

The input file for the whole molecule:
```
%nprocshared=40
%mem=60GB
%chk=mol.chk
# b3lyp 6-31g(d) nosym

Title Card Required

0 1
[Geometry]
```
The geometry of the two fragments are extracted from the whole molecule and written in separate input files.
```
%nprocshared=40
%mem=60GB
%chk=frag1.chk
# b3lyp 6-31g(d) nosym

Title Card Required

0 1
[Geometry]
```
```
%nprocshared=40
%mem=60GB
%chk=frag2.chk
# b3lyp 6-31g(d) nosym

Title Card Required

0 1
[Geometry]
```
Note that the fragments should not cover the whole molecule, because the mismatched space is meaningful and had better exist.

+ Loading necessary packages.
```
from Orbaplaw import WaveFunction as wfn
from Orbaplaw import OrbitalAlignment as oa
```

+ Loading the wavefunction file.
```
mol_mwfn=wfn.MultiWaveFunction("mol.mwfn")
frag1_mwfn=wfn.MultiWaveFunction("frag1.mwfn")
frag2_mwfn=wfn.MultiWaveFunction("frag2.mwfn")
```

+ Conducting FAMO analysis.
```
famo_mwfn=oa.FragmentAlignment(phcl_mwfn,[frag1_mwfn,frag2_mwfn])
```

+ Exporting the FAMOs to `job_famo.mwfn`.
```
famo_mwfn.Export("job_famo.mwfn")
```

## Gallery


[^famo]: Sheong, F. K.; Zhang, J.-X.; Lin, Z. Fragment aligned molecular orbital analysis: An innovative tool for analyzing atypical chemical bonding. *J. Chem. Theory Comput.* **2024**, . [link](https://doi.org/10.1021/acs.jctc.4c00456)