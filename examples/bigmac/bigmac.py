import libmwfn as lm
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo = lm.Mwfn("bigmac.mwfn")
nao, nao_info = nbo.NaturalAtomicOrbital(mo)
nao.Export("bigmac_nao.mwfn")
nho, nbo = nbo.NaturalBondOrbital(nao, nao_info, frags = [
	[0, 1, 2, 11, 12, 13],
	[3, 10],
	[4, 5, 6, 7, 8, 9]
], occ_thres = 0.95)
nho.Export("bigmac_nfho.mwfn")
nbo.Export("bigmac_nfbo.mwfn")

