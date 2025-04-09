from Orbaplaw import WaveFunction as wfn
from Orbaplaw import Localization as loc

polyene_mo = wfn.MultiWaveFunction("polyene.mwfn")

polyene_loc = loc.Localizer(polyene_mo, space="mix", method="ol")
polyene_loc.Export("polyene_loc.mwfn")
