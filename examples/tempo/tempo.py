import libmwfn as lm
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo = lm.Mwfn("tempo.mwfn")
nao, nao_info = nbo.NaturalAtomicOrbital(mo)
nao.Export("tempo_nao.mwfn")
nfho, nfbo = nbo.NaturalBondOrbital(nao, nao_info, frags = [
	[i for i in range(16)], [16]
])
nfho.Export("tempo_nfho.mwfn")
nfbo.Export("tempo_nfbo.mwfn")
pio, pimo = nbo.PrincipalInteractingOrbital(nao, nao_info, [
	[i for i in range(16)], [16]
])
pio.Export("tempo_pio.mwfn")
pimo.Export("tempo_pimo.mwfn")
