import libmwfn as lm
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo = lm.Mwfn("boron12.mwfn")
nao, nao_info = nbo.NaturalAtomicOrbital(mo)
nao.Export("boron12_nao.mwfn")
nfho, nfbo = nbo.NaturalBondOrbital(nao, nao_info, frags = [
	[i for i in range(11)],
	[i for i in range(11, 22)]
])
nfho.Export("boron12_nfho.mwfn")
nfbo.Export("boron12_nfbo.mwfn")
pio, pimo = nbo.PrincipalInteractingOrbital(nao, nao_info, [
	[i for i in range(11)],
	[i for i in range(11, 22)]
])
pio.Export("boron12_pio.mwfn")
pimo.Export("boron12_pimo.mwfn")
