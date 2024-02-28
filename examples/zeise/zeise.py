import numpy as np
from Orbaplaw import WaveFunction as wfn
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo=wfn.MultiWaveFunction("zeise.mwfn")
nao=nbo.NaturalAtomicOrbital(mo)
nao.Export("zeise_nao.mwfn")
nho,nbo=nbo.NaturalBondOrbital(nao,frags=[[0],[1],[2],[3],[4,5,6,7,8,9]],threshold=0.05)
nho.Export("zeise_nfho.mwfn")
nbo.Export("zeise_nfbo.mwfn")
