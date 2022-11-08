"""Define all available cuts"""

# columnar analysis
import awkward as ak
# local
from analysis.definitions.objects import obj_defs


obj_cut_defs = {
    "pvs" : {
        "ndof > 4" : lambda objs: objs["pvs"].ndof > 4,
        "|z| < 24 cm" : lambda objs: abs(objs["pvs"].z) < 24,
        "|rho| < 0.2 mm" : lambda objs: abs(objs["pvs"].rho) < 0.2, # fixme: double check mm
    },
    "ljs" : {
        "pT > 30 GeV" : lambda objs: objs["ljs"].p4.pt > 30,
        "|eta| < 2.4" : lambda objs: abs(objs["ljs"].p4.eta) < 2.4,
        # fixme: figure out way to implement ljs = ljs[:, :2] and similar
    },
}

evt_cut_defs = {
    "PV filter" : lambda objs: ak.num(objs["pvs"]) >= 1,
    "Cosmic veto" : lambda objs: objs["cosmicveto"].result,
    ">=2 LJs" : lambda objs: ak.num(objs["ljs"]) >= 2,
    "4mu" : lambda objs: ak.num(obj_defs["mu_ljs"](objs)) == 2,
    "2mu2e" : lambda objs: (
        (ak.num(obj_defs["mu_ljs"](objs)) == 1) & (ak.num(obj_defs["egm_ljs"](objs)) == 1)
    ),
}
