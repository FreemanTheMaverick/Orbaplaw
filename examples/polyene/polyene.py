import libmwfn as lm
from Orbaplaw import Localization as loc

polyene_mo = lm.Mwfn("polyene.mwfn")

polyene_loc = loc.Localizer(polyene_mo, space="mix", method="ol")
polyene_loc.Export("polyene_loc.mwfn")
