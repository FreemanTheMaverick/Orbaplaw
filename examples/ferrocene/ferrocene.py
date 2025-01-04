from Orbaplaw import WaveFunction as wfn
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo = wfn.MultiWaveFunction("ferrocene.mwfn")
nao = nbo.NaturalAtomicOrbital(mo)
nao.Export("ferrocene_nao.mwfn")
nfho, nfbo = nbo.NaturalBondOrbital(nao,
		frags = [
			[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
			[10],
			[11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
], occ_thres = 0.9)
nfbo.Export("ferrocene_nfbo.mwfn")
nfho.Export("ferrocene_nfho.mwfn")
