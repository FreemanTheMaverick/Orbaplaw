from Orbaplaw import WaveFunction as wfn
from Orbaplaw import NaturalBondOrbitalMethods as nbo
from Orbaplaw import Localization as loc
from Orbaplaw import Population as pop

da_mo = wfn.MultiWaveFunction("da.mwfn")
pop.PopulationAnalyzer(da_mo)

da_loc = loc.Localizer(da_mo, space="occ")
da_loc = loc.Localizer(da_loc, space="vir")
da_loc.Export("da_loc.mwfn")

nao = nbo.NaturalAtomicOrbital(da_mo)
nao.Export("da_nao.mwfn")

nfho, nfbo = nbo.NaturalBondOrbital(nao, frags = [
	[i for i in range(26)],
	[i for i in range(26, 32)]
], occ_thres = 0.9)
nfho.Export("da_nfho.mwfn")
nfbo.Export("da_nfbo.mwfn")

pio, pimo = nbo.PrincipalInteractingOrbital(nao, [
	[i for i in range(26)],
	[i for i in range(26, 32)]
])
pio.Export("da_pio.mwfn")
pimo.Export("da_pimo.mwfn")
