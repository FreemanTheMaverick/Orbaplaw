from Orbaplaw import WaveFunction as wfn
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo=wfn.MultiWaveFunction("naphthalene.mwfn")
nao=nbo.NaturalAtomicOrbital(mo)
nao.Export("naphthalene_nao.mwfn")
nho,nbo=nbo.NaturalBondOrbital(nao,frags=[[i for i in range(10)],[i for i in range(10,20)],[20,21],[22,23]],threshold=0.05)
nho.Export("naphthalene_nfho.mwfn")
nbo.Export("naphthalene_nfbo.mwfn")
