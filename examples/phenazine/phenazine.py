import numpy as np
from Orbaplaw import WaveFunction as wfn
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo=wfn.MultiWaveFunction("phenazine.mwfn")
nao=nbo.NaturalAtomicOrbital(mo)
nao.Export("phenazine_nao.mwfn")
nho,nbo=nbo.NaturalBondOrbital(nao,frags=[[i for i in range(10)],[i for i in range(10,20)],[20],[21]],threshold=0.05)
nho.Export("phenazine_nfho.mwfn")
nbo.Export("phenazine_nfbo.mwfn")
