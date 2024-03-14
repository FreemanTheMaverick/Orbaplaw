from Orbaplaw import WaveFunction as wfn
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo=wfn.MultiWaveFunction("thorium.mwfn")
nao=nbo.NaturalAtomicOrbital(mo)
nao.Export("thorium_nao.mwfn")
nho,nbo=nbo.NaturalBondOrbital(nao,frags=[[0,2],[1,4],[3,5],[i for i in range(6,23)],[i for i in range(23,40)],[i for i in range(40,57)]],threshold=0.1)
nho.Export("thorium_nfho.mwfn")
nbo.Export("thorium_nfbo.mwfn")
