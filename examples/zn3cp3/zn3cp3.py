from Orbaplaw import WaveFunction as wfn
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo = wfn.MultiWaveFunction("zn3cp3.mwfn")
nao = nbo.NaturalAtomicOrbital(mo)
nfho, nfbo = nbo.NaturalBondOrbital(nao, frags = [
	[i for i in range(11)],
	[i for i in range(11, 22)],
	[i for i in range(22, 33)]
])
nfho.Export("zn3cp3_nfho.mwfn")
nfbo.Export("zn3cp3_nfbo.mwfn")
