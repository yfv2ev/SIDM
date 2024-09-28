"""Define all commonly used objects"""

import awkward as ak
from sidm.tools.utilities import matched

# define helper functions
def pid(part, val):
    return part[part.pdgId == val]

def toPid(part, val):
    return part[ak.all(abs(part.children.pdgId) == val, axis=-1)]

def yesMu(lj):
    return lj[lj.muon_n > 0]

def noMu(lj):
    return lj[lj.muon_n == 0]

def noDsa(lj):
    return lj[lj.dsamu_n == 0]

def noPf(lj):
    return lj[lj.pfmu_n == 0]

def noE(lj):
    return lj[lj.electron_n == 0]

def noPhoton(lj):
    return lj[lj.photon_n == 0]

def nE(lj):
    return lj[lj.electron_n == n]

def nPhoton(lj):
    return lj[lj.photon_n == n]

# define objects whose definitions don't depend on LJs
obj_defs = {}
obj_defs["pvs"]        = lambda evts: evts.PV
obj_defs["bs"]         = lambda evts: evts.BS
obj_defs["met"]        = lambda evts: evts.MET
obj_defs["electrons"]  = lambda evts: evts.Electron
obj_defs["photons"]    = lambda evts: evts.Photon
obj_defs["muons"]      = lambda evts: evts.Muon
obj_defs["dsaMuons"]   = lambda evts: evts.DSAMuon
obj_defs["weight"]     = lambda evts: evts.genWeight
obj_defs["gens"]       = lambda evts: evts.GenPart
obj_defs["genMus"]     = lambda evts: pid(obj_defs["gens"](evts), 13)
obj_defs["genEs"]      = lambda evts: pid(obj_defs["gens"](evts), 11)
obj_defs["genAs"]      = lambda evts: pid(obj_defs["gens"](evts), 32)
obj_defs["genAs_toMu"] = lambda evts: toPid(obj_defs["genAs"](evts), 13)
obj_defs["genAs_toE"]  = lambda evts: toPid(obj_defs["genAs"](evts), 11)

# define objects whose definitions rely on LJs
# note that objs["ljs"] will be created in sidm_processor
lj_objs = {}
lj_objs["mu_ljs"]         = lambda objs: yesMu(objs["ljs"])
lj_objs["egm_ljs"]        = lambda objs: noMu(objs["ljs"])
lj_objs["pfmu_ljs"]       = lambda objs: noDsa(lj_objs["mu_ljs"](objs))
lj_objs["dsamu_ljs"]      = lambda objs: noPf(lj_objs["mu_ljs"](objs))
lj_objs["electron_ljs"]   = lambda objs: noPhoton(lj_objs["egm_ljs"](objs))
lj_objs["photon_ljs"]     = lambda objs: noE(lj_objs["egm_ljs"](objs))
lj_objs["n_electron_ljs"] = lambda objs, n: nE(lj_objs["electron_ljs"](objs), n)
lj_objs["n_photon_ljs"]   = lambda objs, n: nPhoton(lj_objs["photon_ljs"](objs), n)
lj_objs["genAs_matched_lj"]        = lambda objs, r: matched(objs["genAs"], objs["ljs"], r)
lj_objs["genAs_toMu_matched_lj"]   = lambda objs, r: matched(objs["genAs_toMu"], objs["ljs"], r)
lj_objs["genAs_toE_matched_lj"]    = lambda objs, r: matched(objs["genAs_toE"], objs["ljs"], r)
lj_objs["genAs_matched_muLj"]      = lambda objs, r: matched(objs["genAs"], lj_objs["mu_ljs"](objs), r)
lj_objs["genAs_toMu_matched_muLj"] = lambda objs, r: matched(objs["genAs_toMu"], lj_objs["mu_ljs"](objs), r)
lj_objs["genAs_matched_egmLj"]     = lambda objs, r: matched(objs["genAs"], lj_objs["egm_ljs"](objs), r)
lj_objs["genAs_toE_matched_egmLj"] = lambda objs, r: matched(objs["genAs_toE"], lj_objs["egm_ljs"](objs), r)
