import libmwfn as lm
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo = lm.Mwfn("thorium.mwfn")
nao, nao_info = nbo.NaturalAtomicOrbital(mo)
nao.Export("thorium_nao.mwfn")
nfho, nfbo = nbo.NaturalBondOrbital(nao, nao_info, frags = [
		[0, 2], [1, 4], [3, 5],
		[i for i in range(6, 23)],
		[i for i in range(23, 40)],
		[i for i in range(40, 57)]
], occ_thres = 0.9)
nfho.Export("thorium_nfho.mwfn")
nfbo.Export("thorium_nfbo.mwfn")
