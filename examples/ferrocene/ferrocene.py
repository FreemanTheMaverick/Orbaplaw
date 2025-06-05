import libmwfn as lm
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo = lm.Mwfn("ferrocene.mwfn")
nao, nao_info = nbo.NaturalAtomicOrbital(mo)
nao.Export("ferrocene_nao.mwfn")
nfho, nfbo = nbo.NaturalBondOrbital(nao, nao_info, frags = [
			[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
			[10],
			[11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
], occ_thres = 0.9)
nfbo.Export("ferrocene_nfbo.mwfn")
nfho.Export("ferrocene_nfho.mwfn")
