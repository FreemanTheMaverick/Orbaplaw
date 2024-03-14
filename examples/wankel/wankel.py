from Orbaplaw import WaveFunction as wfn
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo=wfn.MultiWaveFunction("wankel.mwfn")
nao=nbo.NaturalAtomicOrbital(mo)
nao.Export("wankel_nao.mwfn")
nho,nbo=nbo.NaturalBondOrbital(nao,frags=[[0],[i for i in range(1,6)],[i for i in range(6,19)]],threshold=0.1)
nho.Export("wankel_nfho.mwfn")
nbo.Export("wankel_nfbo.mwfn")

