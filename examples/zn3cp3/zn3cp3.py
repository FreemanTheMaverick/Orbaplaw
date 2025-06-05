import libmwfn as lm
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo = lm.Mwfn("zn3cp3.mwfn")
nao, nao_info = nbo.NaturalAtomicOrbital(mo)
nfho, nfbo = nbo.NaturalBondOrbital(nao, nao_info, frags = [
	[i for i in range(11)],
	[i for i in range(11, 22)],
	[i for i in range(22, 33)]
])
nfho.Export("zn3cp3_nfho.mwfn")
nfbo.Export("zn3cp3_nfbo.mwfn")
