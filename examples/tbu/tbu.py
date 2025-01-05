from Orbaplaw import WaveFunction as wfn
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo = wfn.MultiWaveFunction("tbu.mwfn")
nao = nbo.NaturalAtomicOrbital(mo)
nao.Export("tbu_nao.mwfn")
nfho, nfbo = nbo.NaturalBondOrbital(nao, frags = [
	[1, 2, 3, 4],
	[0] + [i for i in range(5, 13)]
])
nfho.Export("tbu_nfho.mwfn")
nfbo.Export("tbu_nfbo.mwfn")
pio, pimo = nbo.PrincipalInteractingOrbital(nao, [
	[1, 2, 3, 4],
	[0] + [i for i in range(5, 13)]
])
pio.Export("tbu_pio.mwfn")
pimo.Export("tbu_pimo.mwfn")
