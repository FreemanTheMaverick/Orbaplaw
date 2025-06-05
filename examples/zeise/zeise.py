import libmwfn as lm
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo = lm.Mwfn("zeise.mwfn")
print("The three-center-four-electron scheme:")
nao, nao_info = nbo.NaturalAtomicOrbital(mo)
nao.Export("zeise_3c4e_nao.mwfn")
nfho, nfbo = nbo.NaturalBondOrbital(nao, nao_info, frags = [
	[0], [1], [2], [3],
	[4, 5, 6, 7, 8, 9]
], multi_thres = 0.90)
nfho.Export("zeise_3c4e_nfho.mwfn")
nfbo.Export("zeise_3c4e_nfbo.mwfn")
print("The hybridized-orbital scheme:")
nao, nao_info = nbo.NaturalAtomicOrbital(mo, modify_minimal_shells = [(0, [9, 9, 9, 9, 9, 9, 9])])
nao.Export("zeise_hybridized_nao.mwfn")
nfho, nfbo = nbo.NaturalBondOrbital(nao, nao_info, frags = [
	[0], [1], [2], [3],
	[4, 5, 6, 7, 8, 9]
], multi_thres = 0.95, deg_thres = 0.05)
nfho.Export("zeise_hybridized_nfho.mwfn")
nfbo.Export("zeise_hybridized_nfbo.mwfn")
