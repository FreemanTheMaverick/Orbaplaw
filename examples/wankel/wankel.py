import libmwfn as lm
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo = lm.Mwfn("wankel.mwfn")
nao, nao_info = nbo.NaturalAtomicOrbital(mo)
nao.Export("wankel_nao.mwfn")
nho, nbo = nbo.NaturalBondOrbital(nao, nao_info, frags = [
	[0],
	[i for i in range(1, 6)],
	[i for i in range(6, 19)]
], occ_thres = 0.9)
nho.Export("wankel_nfho.mwfn")
nbo.Export("wankel_nfbo.mwfn")
