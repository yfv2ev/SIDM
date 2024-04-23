"""Define all commonly used objects"""

import awkward as ak
from sidm.tools.utilities import dR, matched

# define objects whose definitions depend only on the event record
primary_objs = {
    "cosmicveto": lambda evts: evts.cosmicveto,
    "pvs": lambda evts: evts.pv,
    "electrons": lambda evts: evts.electron,
    "photons": lambda evts: evts.pfphoton, # fixme: understand differences between photon and pfphoton in v4 ntuples
    "muons": lambda evts: evts.muon,
    "dsaMuons": lambda evts: evts.dsamuon,
    "ntuple_ljs": lambda evts: evts.pfjet,
    "ljsources": lambda evts: evts.ljsource,
    "gens": lambda evts: evts.gen,
    "genEs": lambda evts: evts.gen[abs(evts.gen.pid) == 11],
    "genMus": lambda evts: evts.gen[abs(evts.gen.pid) == 13],
    "genAs": lambda evts: evts.gen[abs(evts.gen.pid) == 32],
    "genAs_toMu": lambda evts: evts.gen[(abs(evts.gen.pid)== 32) & (abs(evts.gen.daupid) == 13)],
    "genAs_toE": lambda evts: evts.gen[(abs(evts.gen.pid)== 32) & (abs(evts.gen.daupid) == 11)],
    "weight" : lambda evts: evts.weightProduct,
}

llpNanoAod_objs = {
    "pvs": lambda evts: evts.PV,
    "electrons": lambda evts: evts.Electron,
    "photons": lambda evts: evts.Photon,
    "muons" : lambda evts: evts.Muon,
    "dsaMuons" : lambda evts: evts.DSAMuon,
    "gens": lambda evts: evts.GenPart,
    "genEs": lambda evts: evts.GenPart[abs(evts.GenPart.pdgId) == 11],
    "genMus": lambda evts: evts.GenPart[abs(evts.GenPart.pdgId) == 13],
    "genAs": lambda evts: evts.GenPart[abs(evts.GenPart.pdgId) == 32],
    "genAs_toMu": lambda evts: evts.GenPart[(abs(evts.GenPart.pdgId)== 32) & ak.all(abs(evts.GenPart.children.pdgId) == 13, axis=-1)],
    "genAs_toE": lambda evts: evts.GenPart[(abs(evts.GenPart.pdgId)== 32) & ak.all(abs(evts.GenPart.children.pdgId) == 11, axis=-1)],
    "weight" : lambda evts: evts.genWeight,
}

# define objects whose definitions depend on analysis choices
derived_objs = {
    "mu_ljs": lambda objs: objs["ljs"][(objs["ljs"].muon_n >= 2)],
    "egm_ljs": lambda objs: objs["ljs"][(objs["ljs"].muon_n == 0)],
    "electron_ljs": lambda objs, n: objs["ljs"][(objs["ljs"].muon_n == 0) & (objs["ljs"].photon_n == 0) & (objs["ljs"].electron_n == n)],
    "photon_ljs": lambda objs, n: objs["ljs"][(objs["ljs"].muon_n == 0) & (objs["ljs"].photon_n == n) & (objs["ljs"].electron_n == 0)],
    "genAs_matched_lj": lambda objs, r: matched(objs["genAs"], objs["ljs"], r),
    "genAs_toMu_matched_lj": lambda objs, r: matched(objs["genAs_toMu"], objs["ljs"], r),
    "genAs_toE_matched_lj": lambda objs, r: matched(objs["genAs_toE"], objs["ljs"], r),
    "genAs_matched_muLj": lambda objs, r: matched(objs["genAs"], objs["ljs"][(objs["ljs"].muon_n >= 2)], r),
    "genAs_toMu_matched_muLj": lambda objs, r: matched(objs["genAs_toMu"], objs["ljs"][(objs["ljs"].muon_n >= 2)], r),
    "genAs_matched_egmLj": lambda objs, r: matched(objs["genAs"], objs["ljs"][(objs["ljs"].muon_n == 0)], r),
    "genAs_toE_matched_egmLj": lambda objs, r: matched(objs["genAs_toE"], objs["ljs"][(objs["ljs"].muon_n == 0)], r),
}

