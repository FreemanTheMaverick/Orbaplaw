from Orbaplaw import WaveFunction as wfn
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo=wfn.MultiWaveFunction("boron12.mwfn")
nao=nbo.NaturalAtomicOrbital(mo)
nao.Export("boron12_nao.mwfn")
nho,nbo=nbo.NaturalBondOrbital(nao,frags=[[i for i in range(11)],[i for i in range(11,22)]],threshold=0.05)
nho.Export("boron12_nfho.mwfn")
nbo.Export("boron12_nfbo.mwfn")

