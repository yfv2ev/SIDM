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
preLj_objs = {}
preLj_objs["pvs"]        = lambda evts: evts.PV
preLj_objs["bs"]         = lambda evts: evts.BS
preLj_objs["met"]        = lambda evts: evts.MET
preLj_objs["electrons"]  = lambda evts: evts.Electron
preLj_objs["photons"]    = lambda evts: evts.Photon
preLj_objs["muons"]      = lambda evts: evts.Muon
preLj_objs["dsaMuons"]   = lambda evts: evts.DSAMuon
preLj_objs["weight"]     = lambda evts: evts.genWeight
preLj_objs["gens"]       = lambda evts: evts.GenPart
preLj_objs["genMus"]     = lambda evts: pid(preLj_objs["gens"](evts), 13)
preLj_objs["genEs"]      = lambda evts: pid(preLj_objs["gens"](evts), 11)
preLj_objs["genAs"]      = lambda evts: pid(preLj_objs["gens"](evts), 32)
preLj_objs["genAs_toMu"] = lambda evts: toPid(preLj_objs["genAs"](evts), 13)
preLj_objs["genAs_toE"]  = lambda evts: toPid(preLj_objs["genAs"](evts), 11)

# define objects whose definitions rely on LJs
# note that objs["ljs"] will be created in sidm_processor
postLj_objs = {}
postLj_objs["mu_ljs"]         = lambda objs: yesMu(objs["ljs"])
postLj_objs["egm_ljs"]        = lambda objs: noMu(objs["ljs"])
postLj_objs["pfmu_ljs"]       = lambda objs: noDsa(postLj_objs["mu_ljs"](objs))
postLj_objs["dsamu_ljs"]      = lambda objs: noPf(postLj_objs["mu_ljs"](objs))
postLj_objs["electron_ljs"]   = lambda objs: noPhoton(postLj_objs["egm_ljs"](objs))
postLj_objs["photon_ljs"]     = lambda objs: noE(postLj_objs["egm_ljs"](objs))
postLj_objs["n_electron_ljs"] = lambda objs, n: nE(postLj_objs["electron_ljs"](objs), n)
postLj_objs["n_photon_ljs"]   = lambda objs, n: nPhoton(postLj_objs["photon_ljs"](objs), n)
postLj_objs["genAs_matched_lj"]        = lambda objs, r: matched(objs["genAs"], objs["ljs"], r)
postLj_objs["genAs_toMu_matched_lj"]   = lambda objs, r: matched(objs["genAs_toMu"], objs["ljs"], r)
postLj_objs["genAs_toE_matched_lj"]    = lambda objs, r: matched(objs["genAs_toE"], objs["ljs"], r)
postLj_objs["genAs_matched_muLj"]      = lambda objs, r: matched(objs["genAs"], postLj_objs["mu_ljs"](objs), r)
postLj_objs["genAs_toMu_matched_muLj"] = lambda objs, r: matched(objs["genAs_toMu"], postLj_objs["mu_ljs"](objs), r)
postLj_objs["genAs_matched_egmLj"]     = lambda objs, r: matched(objs["genAs"], postLj_objs["egm_ljs"](objs), r)
postLj_objs["genAs_toE_matched_egmLj"] = lambda objs, r: matched(objs["genAs_toE"], postLj_objs["egm_ljs"](objs), r)
