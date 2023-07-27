"""Define all commonly used objects"""

from sidm.tools.utilities import dR

# define objects whose definitions depend only on the event record
primary_objs = {
    "cosmicveto": lambda evts: evts.cosmicveto,
    "pvs": lambda evts: evts.pv,
    "electrons": lambda evts: evts.electron,
    "photons": lambda evts: evts.pfphoton,
    "muons": lambda evts: evts.muon,
    "dsaMuons": lambda evts: evts.dsamuon,
    "ljs": lambda evts: evts.pfjet,
    "ljsources": lambda evts: evts.ljsource,
    "gens": lambda evts: evts.gen,
    "genEs": lambda evts: evts.gen[abs(evts.gen.pid) == 11],
    "genMus": lambda evts: evts.gen[abs(evts.gen.pid) == 13],
    "genAs": lambda evts: evts.gen[abs(evts.gen.pid) == 32],
    "genAs_toMu": lambda evts : evts.gen[(abs(evts.gen.pid)== 32) & (abs(evts.gen.daupid) == 13)],
    "genAs_toE": lambda evts : evts.gen[(abs(evts.gen.pid)== 32) & (abs(evts.gen.daupid) == 11)],
}

# define objects whose definitions depend on analysis choices
derived_objs = {
    "mu_ljs": lambda objs: objs["ljs"][(objs["ljs"].muon_n >= 2)],
    "egm_ljs": lambda objs: objs["ljs"][(objs["ljs"].muon_n == 0)],
    "matched_genAs": lambda objs, r: objs["genAs"][dR(objs["genAs"], objs["ljs"]) < r]
    "matched_genAs_mu": lambda objs, r: objs["genAs"][dR(objs["genAs"], objs["ljs"][(objs["ljs"].muon_n >= 2)]) < r],
}
