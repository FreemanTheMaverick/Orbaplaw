from Orbaplaw import WaveFunction as wfn
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo=wfn.MultiWaveFunction("bigmac.mwfn")
nao=nbo.NaturalAtomicOrbital(mo)
nao.Export("bigmac_nao.mwfn")
nho,nbo=nbo.NaturalBondOrbital(nao,frags=[[0,1,2,11,12,13],[3,10],[4,5,6,7,8,9]],threshold=0.05)
nho.Export("bigmac_nfho.mwfn")
nbo.Export("bigmac_nfbo.mwfn")

