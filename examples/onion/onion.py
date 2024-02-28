from Orbaplaw import WaveFunction as wfn
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo=wfn.MultiWaveFunction("onion.mwfn")
nao=nbo.NaturalAtomicOrbital(mo)
nho,nbo=nbo.NaturalBondOrbital(nao,frags=[[0],[i for i in range(1,13)],[i for i in range(13,33)]],threshold=0.025)
nho.Export("onion_nfho.mwfn")
nbo.Export("onion_nfbo.mwfn")

