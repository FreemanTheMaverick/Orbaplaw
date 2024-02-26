import numpy as np
from Orbaplaw import WaveFunction as wfn
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo=wfn.MultiWaveFunction("da.mwfn")
nao=nbo.NaturalAtomicOrbital(mo)
nao.Export("da_nao.mwfn")
nho,nbo=nbo.NaturalBondOrbital(nao,frags=[[i for i in range(26)],[26,27,28,29,30,31]])
nho.Export("da_nfho.mwfn")
nbo.Export("da_nfbo.mwfn")
