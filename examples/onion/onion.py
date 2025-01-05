from Orbaplaw import WaveFunction as wfn
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo = wfn.MultiWaveFunction("onion.mwfn")
nao = nbo.NaturalAtomicOrbital(mo)
nfho, nfbo = nbo.NaturalBondOrbital(nao, frags = [
	[0],
	[i for i in range(1, 33)],
])
nfho.Export("onion_nfho.mwfn")
nfbo.Export("onion_nfbo.mwfn")
pio, pimo = nbo.PrincipalInteractingOrbital(nao, [
	[0],
	[i for i in range(1, 33)],
])
pio.Export("onion_pio.mwfn")
pimo.Export("onion_pimo.mwfn")
