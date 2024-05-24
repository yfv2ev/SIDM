"""Define all available cuts"""

# columnar analysis
import awkward as ak
# local
from sidm.definitions.objects import derived_objs
from sidm.tools.utilities import dR, dR_outer, lxy, check_bit, as_int


obj_cut_defs = {
    "pvs": {
        "ndof > 4": lambda objs: objs["pvs"].ndof > 4,
        "|z| < 24 cm": lambda objs: abs(objs["pvs"].z) < 24,
        "|rho| < 0.2 mm": lambda objs: abs(objs["pvs"].pos.rho) < 0.2,
    },
    "ljs": {
        "pT > 30 GeV": lambda objs: objs["ljs"].pt > 30,
        "|eta| < 2.4": lambda objs: abs(objs["ljs"].eta) < 2.4,
        "dR(LJ, A) < 0.2": lambda objs: dR(objs["ljs"], objs["genAs"]) < 0.2,
        "eLj": lambda objs: objs["ljs"].electron_n > 0,
        "gLj": lambda objs: (objs["ljs"].photon_n > 0 ) & (objs["ljs"].electron_n == 0),
        "muLj": lambda objs: objs["ljs"][(objs["ljs"].muon_n >= 2)],
        "dsaMuLj": lambda objs: (objs["ljs"].muon_n == 2) & (ak.any(objs["ljs"].pfcand['type'] == 8, axis=-1)),
        "2dsaMuLj": lambda objs: (objs["ljs"].muon_n == 2) & (ak.all(objs["ljs"].pfcand['type'] == 8, axis=-1)),
        "1dsa1pfMuLj": lambda objs: ((objs["ljs"].muon_n == 2)
                                     & ((ak.any(objs["ljs"].pfcand['type'] == 8, axis=-1))
                                     & (ak.any(objs["ljs"].pfcand['type'] == 3, axis=-1)))),
        "pfMuLj": lambda objs: (objs["ljs"].muon_n == 2) & (ak.all(objs["ljs"].pfcand['type'] != 8, axis=-1)),
    },
    "genAs": {
        "dR(A, LJ) < 0.2": lambda objs: dR(objs["genAs"], objs["ntuple_ljs"]) < 0.2,
        "dR(A, LJ) < 0.4": lambda objs: dR(objs["genAs"], objs["ntuple_ljs"]) < 0.4,
        "lxy < 10 cm": lambda objs: lxy(objs["genAs"]) < 10,
        "lxy < 40 cm": lambda objs: lxy(objs["genAs"]) < 40,
        "10 cm <= lxy < 100 cm": lambda objs: ((lxy(objs["genAs"]) >= 10)
                                               & (lxy(objs["genAs"]) < 100)),
        "100 cm <= lxy < 135 cm": lambda objs: ((lxy(objs["genAs"]) >= 100)
                                               & (lxy(objs["genAs"]) < 135)),
        "lxy >= 100 cm": lambda objs: lxy(objs["genAs"]) >= 100,
        "lxy <= 100 cm": lambda objs: lxy(objs["genAs"]) <= 100,
        "lxy <= 150 cm": lambda objs: lxy(objs["genAs"]) <= 150,
        "lxy <= 250 cm": lambda objs: lxy(objs["genAs"]) <= 250,
        "pT > 30 GeV": lambda objs: objs["genAs"].pt > 30,
        "pT < 300 GeV": lambda objs: objs["genAs"].pt < 300,
    },
    "genAs_toMu": {
        "dR(A, LJ) < 0.2": lambda objs: dR(objs["genAs_toMu"], objs["ntuple_ljs"]) < 0.2,
        "dR(A, LJ) < 0.4": lambda objs: dR(objs["genAs_toMu"], objs["ntuple_ljs"]) < 0.4,
        "lxy < 10 cm": lambda objs: lxy(objs["genAs_toMu"]) < 10,
        "10 cm <= lxy < 100 cm": lambda objs: (lxy(objs["genAs_toMu"]) >= 10) & (lxy(objs["genAs_toMu"]) < 100),
        "lxy >= 100 cm": lambda objs: lxy(objs["genAs_toMu"]) >= 100,
        "lxy <= 100 cm": lambda objs: lxy(objs["genAs_toMu"]) <= 100,
        "lxy <= 150 cm": lambda objs: lxy(objs["genAs_toMu"]) <= 150,
        "lxy <= 250 cm": lambda objs: lxy(objs["genAs_toMu"]) <= 250,
        "lxy <= 400 cm": lambda objs: lxy(objs["genAs_toMu"]) <= 400,
        "pT > 30 GeV": lambda objs: objs["genAs_toMu"].pt > 30,
        "pT < 300 GeV": lambda objs: objs["genAs_toMu"].pt < 300,
    },
    "genAs_toE": {
        "dR(A, LJ) < 0.2": lambda objs: dR(objs["genAs_toE"], objs["ntuple_ljs"]) < 0.2,
        "dR(A, LJ) < 0.4": lambda objs: dR(objs["genAs_toE"], objs["ntuple_ljs"]) < 0.4,
        "lxy <= 5 cm": lambda objs: lxy(objs["genAs_toE"]) <= 5,
        "lxy <= 2.5 cm": lambda objs: lxy(objs["genAs_toE"]) <= 2.5,
        "lxy < 10 cm": lambda objs: lxy(objs["genAs_toE"]) < 10,
        "lxy < 40 cm": lambda objs: lxy(objs["genAs_toE"]) < 40,
        "10 cm <= lxy < 100 cm": lambda objs: (lxy(objs["genAs_toE"]) >= 10) & (lxy(objs["genAs_toE"]) < 100),
        "40 cm <= lxy < 77 cm": lambda objs: (lxy(objs["genAs_toE"]) >= 40) & (lxy(objs["genAs_toE"]) < 77),
        "100 cm <= lxy < 135 cm": lambda objs: (lxy(objs["genAs_toE"]) >= 100) & (lxy(objs["genAs_toE"]) < 135),
        "125 cm <= lxy <= 135 cm": lambda objs: (lxy(objs["genAs_toE"]) >= 125) & (lxy(objs["genAs_toE"]) < 135),
        "lxy >= 100 cm": lambda objs: lxy(objs["genAs_toE"]) >= 100,
        "lxy <= 100 cm": lambda objs: lxy(objs["genAs_toE"]) <= 100,
        "lxy <= 150 cm": lambda objs: lxy(objs["genAs_toE"]) <= 250,
        "lxy <= 250 cm": lambda objs: lxy(objs["genAs_toE"]) <= 250,
        "pT > 30 GeV": lambda objs: objs["genAs_toE"].pt > 30,
        "pT < 300 GeV": lambda objs: objs["genAs_toE"].pt < 300,
    },
    "electrons": {
        "pT > 10 GeV": lambda objs: objs["electrons"].pt > 10,
        "|eta| < 1.479": lambda objs: abs(objs["electrons"].eta) < 1.479,
        "1.479 < |eta| < 2.4": lambda objs: ((abs(objs["electrons"].eta) > 1.479)
                                             & (abs(objs["electrons"].eta) < 2.4)),
        "|eta| < 2.4": lambda objs: abs(objs["electrons"].eta) < 2.4,
        "dR(e, A) < 0.5": lambda objs: dR(objs["electrons"], objs["genAs_toE"]) < 0.5,
        #Loose ID = bit 1
        "looseID": lambda objs: objs["electrons"].cutBased == 2, # (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
        "barrel SigmaIEtaIEtaCut": lambda objs: objs["electrons"].sieie < 0.0112,
        "barrel DEtaInSeedCut": lambda objs: abs(objs["electrons"].deltaEtaSC) < 0.00377,
        "barrel DPhiInCut": lambda objs: abs(objs["electrons"].GsfEleDPhiInCut_0) < 0.0884, # fixme: can't find dPhiIn in nanoAOD
        "barrel InverseCut": lambda objs: abs(objs["electrons"].eInvMinusPInv) < 0.193,
        "barrel Iso": lambda objs: (objs["electrons"].pfRelIso03_all
                                    < (0.112 + 0.506/(objs["electrons"].pt))),
        "barrel ConversionVeto": lambda objs: objs["electrons"].convVeto,
        "barrel H/E": lambda objs: objs["electrons"].hoe < 0.05,
        "barrel MissingHits": lambda objs: (abs(objs["electrons"].lostHits) < 1),
    },
    "muons": {
        "looseID": lambda objs: objs["muons"].looseId,
        "pT > 5 GeV": lambda objs: objs["muons"].pt > 5,
        "|eta| < 2.4": lambda objs: abs(objs["muons"].eta) < 2.4,
    },
    "photons":{
        "pT > 20 GeV": lambda objs: objs["photons"].pt > 20,
        "|eta| < 2.5": lambda objs: abs(objs["photons"].eta) < 2.5, # fixme: do we want eta or scEta
        #Loose ID = bit 0
        "looseID": lambda objs: objs["photons"].cutBased == 2,
    },
    "dsaMuons": {
        "pT > 10 GeV": lambda objs: objs["dsaMuons"].pt > 10,
        "|eta| < 2.4": lambda objs: abs(objs["dsaMuons"].eta) < 2.4,
        "DT + CSC hits > 12": lambda objs: (objs["dsaMuons"].trkNumDTHits
                                            + objs["dsaMuons"].trkNumCSCHits) > 12,
        "ifcsczero": lambda objs: ak.where(((objs["dsaMuons"].trkNumCSCHits == 0) 
                                           & (objs["dsaMuons"].trkNumDTHits <= 18)), False, True),
        "normChi2 < 2.5": lambda objs: objs["dsaMuons"].normChi2 < 2.5,
        "ptErrorOverPT < 1": lambda objs: (objs["dsaMuons"].ptErr / objs["dsaMuons"].pt) < 1.0,
        # I can't find all the necessary info to properly match DSA to PF; for now only consider
        # PF with most matched segments and require >=2 matched segments OR dR_outer < 0.1
        "no PF match" : lambda objs: ((objs["dsaMuons"].muonMatch1 < 2) & (
            dR_outer(objs["dsaMuons"], objs["muons"][as_int(objs["dsaMuons"].muonMatch1idx)]) > 0.1)),
    }
}

evt_cut_defs = {
    # This following will be True for every event. There's probably a more intuitive way to do this
    "Keep all evts": lambda objs: objs["pvs"].npvs >= 0,
    ">=1 muon": lambda objs: ak.num(objs["muons"]) >= 1,
    "PV filter": lambda objs: objs["pvs"].npvs >= 1,
    #"Cosmic veto": lambda objs: objs["cosmicveto"].result,
    ">=2 LJs": lambda objs: ak.num(objs["ljs"]) >= 2,
    ">=2 matched As": lambda objs: ak.num(derived_objs["genAs_matched_lj"](objs, 0.2)) >= 2,
    # 4mu: leading two LJs are both mu-type
    "4mu": lambda objs: ak.count_nonzero(objs["ljs"][:, :2].muon_n >= 2, axis=-1) == 2,
    # 2mu2e: leading two LJs contain exactly 1 mu-type and exactly 1 egm-type
    "2mu2e": lambda objs: ((ak.count_nonzero(objs["ljs"][:, :2].muon_n >= 2, axis=-1) == 1)
                           & (ak.count_nonzero(objs["ljs"][:, :2].muon_n == 0, axis=-1) == 1)),
    "genAs_toE_matched_egmLj": lambda objs: ak.num(derived_objs["genAs_toE_matched_egmLj"](objs, 0.4)) >= 1,
    "genAs_toMu_matched_muLj": lambda objs: ak.num(derived_objs["genAs_toMu_matched_muLj"](objs,0.4)) >= 1,
    "genAs_toE": lambda objs: ak.num(objs["genAs_toE"]) >= 1,
    "genAs_toMu": lambda objs: ak.num(objs["genAs_toMu"]) >= 1,           
    "ljs": lambda objs: ak.num(objs["ntuple_ljs"]) >= 1,           
}
