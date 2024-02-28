from Orbaplaw import WaveFunction as wfn
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo=wfn.MultiWaveFunction("carbon_gold.mwfn")
nao=nbo.NaturalAtomicOrbital(mo)
nao.Export("carbon_gold_nao.mwfn")
nho,nbo=nbo.NaturalBondOrbital(nao,frags=[[0],[i for i in range(1,31)]],threshold=0.025)
nho.Export("carbon_gold_nfho.mwfn")
nbo.Export("carbon_gold_nfbo.mwfn")

