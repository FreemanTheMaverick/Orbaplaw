from Orbaplaw import WaveFunction as wfn
from Orbaplaw import NaturalBondOrbitalMethods as nbo

mo = wfn.MultiWaveFunction("thorium.mwfn")
nao = nbo.NaturalAtomicOrbital(mo)
nao.Export("thorium_nao.mwfn")
nfho, nfbo = nbo.NaturalBondOrbital(nao, frags = [
		[0, 2], [1, 4], [3, 5],
		[i for i in range(6, 23)],
		[i for i in range(23, 40)],
		[i for i in range(40, 57)]
], occ_thres = 0.9)
nfho.Export("thorium_nfho.mwfn")
nfbo.Export("thorium_nfbo.mwfn")
