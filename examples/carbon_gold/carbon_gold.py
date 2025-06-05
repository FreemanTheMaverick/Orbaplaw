import libmwfn as lm
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo = lm.Mwfn("carbon_gold.mwfn")
nao, nao_info = nbo.NaturalAtomicOrbital(mo)
nao.Export("carbon_gold_nao.mwfn")
nfho, nfbo = nbo.NaturalBondOrbital(nao, nao_info, frags = [
	[0],
	[i for i in range(1, 31)]
])
nfho.Export("carbon_gold_nfho.mwfn")
nfbo.Export("carbon_gold_nfbo.mwfn")
pio, pimo = nbo.PrincipalInteractingOrbital(nao, nao_info, [
	[0],
	[i for i in range(1, 31)]
])
pio.Export("carbon_gold_pio.mwfn")
pimo.Export("carbon_gold_pimo.mwfn")
