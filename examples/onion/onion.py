import libmwfn as lm
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo = lm.Mwfn("onion.mwfn")
nao, nao_info = nbo.NaturalAtomicOrbital(mo)
nfho, nfbo = nbo.NaturalBondOrbital(nao, nao_info, frags = [
	[0],
	[i for i in range(1, 33)],
])
nfho.Export("onion_nfho.mwfn")
nfbo.Export("onion_nfbo.mwfn")
pio, pimo = nbo.PrincipalInteractingOrbital(nao, nao_info, [
	[0],
	[i for i in range(1, 33)],
])
pio.Export("onion_pio.mwfn")
pimo.Export("onion_pimo.mwfn")
