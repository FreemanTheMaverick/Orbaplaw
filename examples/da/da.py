from Orbaplaw import WaveFunction as wfn
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo = wfn.MultiWaveFunction("da.mwfn")
nao = nbo.NaturalAtomicOrbital(mo)
nao.Export("da_nao.mwfn")
nfho, nfbo = nbo.NaturalBondOrbital(nao, frags = [[i for i in range(26)],[26,27,28,29,30,31]], occ_thres = 0.9)
nfho.Export("da_nfho.mwfn")
nfbo.Export("da_nfbo.mwfn")
pio, pimo = nbo.PrincipalInteractingOrbital(nao, [[i for i in range(26)],[26,27,28,29,30,31]])
pio.Export("da_pio.mwfn")
pimo.Export("da_pimo.mwfn")
