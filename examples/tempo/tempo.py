from Orbaplaw import WaveFunction as wfn
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo = wfn.MultiWaveFunction("tempo.mwfn")
nao = nbo.NaturalAtomicOrbital(mo)
nao.Export("tempo_nao.mwfn")
nfho, nfbo = nbo.NaturalBondOrbital(nao, frags = [[i for i in range(16)], [16]])
nfho.Export("tempo_nfho.mwfn")
nfbo.Export("tempo_nfbo.mwfn")
pio, pimo = nbo.PrincipalInteractingOrbital(nao, [[i for i in range(16)], [16]])
pio.Export("tempo_pio.mwfn")
pimo.Export("tempo_pimo.mwfn")
