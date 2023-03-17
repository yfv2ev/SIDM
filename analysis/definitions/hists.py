"""Define all available histograms

All hists are defined as Histogram objects whose axes are given as a list of Axis objects, which
bundle a hist.axis with a function that defines how the axis will be filled. The underlying
hist.Hists storage is weight unless otherwise specified.
"""

# python
import math
import importlib
# columnar analysis
import hist
import awkward as ak
# local
from analysis.tools import histogram as h
from analysis.tools.utilities import dR
from analysis.definitions.objects import obj_defs
# always reload local modules to pick up changes during development
importlib.reload(h)


common_axes = {
    "lj_pt": hist.axis.Regular(100, 0, 100, name="lj_pt", label="Lepton jet pT [GeV]")
}

hist_defs = {
    # pv
    "pv_n": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 100, name="pv_n"),
                   lambda objs: ak.num(objs["pvs"])),
        ],
    ),
    "pv_ndof": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0, 100, name="pv_ndof"),
                   lambda objs: objs["pvs"].ndof),
        ],
    ),
    "pv_z": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, -50, 50, name="pv_z"),
                   lambda objs: objs["pvs"].z),
        ],
    ),
    "pv_rho": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, -0.5, 0.5, name="pv_rho"),
                   lambda objs: objs["pvs"].rho),
        ],
    ),
    # pfelectron
    "electron_n": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="electron_n"),
                   lambda objs: ak.num(objs["electrons"])),
        ],
    ),
    "electron_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="electron_pt"),
                   lambda objs: objs["electrons"].p4.pt),
        ],
    ),
    "electron_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="electron_eta"),
                   lambda objs: objs["electrons"].p4.eta),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="electron_phi"),
                   lambda objs: objs["electrons"].p4.phi),
        ],
    ),
    # pfphoton
    "photon_n": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="photon_n"),
                   lambda objs: ak.num(objs["photons"])),
        ],
    ),
    "photon_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="photon_pt"),
                   lambda objs: objs["photons"].p4.pt),
        ],
    ),
    "photon_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="photon_eta"),
                   lambda objs: objs["photons"].p4.eta),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="photon_phi"),
                   lambda objs: objs["photons"].p4.phi),
        ],
    ),
    # pfmuon
    "muon_n": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="muon_n"),
                   lambda objs: ak.num(objs["muons"])),
        ],
    ),
    "muon_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="muon_pt"),
                   lambda objs: objs["muons"].p4.pt),
        ],
    ),
    "muon_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="muon_eta"),
                   lambda objs: objs["muons"].p4.eta),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="muon_phi"),
                   lambda objs: objs["muons"].p4.phi),
        ],
    ),
    # dsamuon
    "dsaMuon_n": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="dsaMuon_n"),
                   lambda objs: ak.num(objs["dsaMuons"])),
        ],
    ),
    "dsaMuon_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="dsaMuon_pt"),
                   lambda objs: objs["dsaMuons"].p4.pt),
        ],
    ),
    "dsaMuon_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="dsaMuon_eta"),
                   lambda objs: objs["dsaMuons"].p4.eta),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="dsaMuon_phi"),
                   lambda objs: objs["dsaMuons"].p4.phi),
        ],
    ),
    # lj
    "lj_n": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="lj_n"),
                   lambda objs: ak.num(objs["ljs"])),
        ],
    ),
    "lj_pt": h.Histogram(
        [
            h.Axis(common_axes["lj_pt"],
                   lambda objs: objs["ljs"].p4.pt),
        ],
    ),
    "lj0_pt": h.Histogram(
        [
            h.Axis(common_axes["lj_pt"],
                   lambda objs: objs["ljs"][ak.num(objs["ljs"]) > 0, 0].p4.pt),
        ],
    ),
    "lj1_pt": h.Histogram(
        [
            h.Axis(common_axes["lj_pt"],
                   lambda objs: objs["ljs"][ak.num(objs["ljs"]) > 1, 1].p4.pt),
        ],
    ),
    "lj_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="lj_eta"),
                   lambda objs: objs["ljs"].p4.eta),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="lj_phi"),
                   lambda objs: objs["ljs"].p4.phi),
        ],
    ),
    "egm_lj_pt": h.Histogram(
        [
            h.Axis(common_axes["lj_pt"],
                   lambda objs: obj_defs["egm_ljs"](objs).p4.pt),
        ],
    ),
    "mu_lj_pt": h.Histogram(
        [
            h.Axis(common_axes["lj_pt"],
                   lambda objs: obj_defs["mu_ljs"](objs).p4.pt),
        ],
    ),
    "lj_electronN": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="lj_electronN"),
                   lambda objs: objs["ljs"].electron_n),
        ],
    ),
    "lj_photonN": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="lj_photonN"),
                   lambda objs: objs["ljs"].photon_n),
        ],
    ),
    "lj_electronPhotonN": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="lj_electronPhotonN"),
                   lambda objs: objs["ljs"].electron_n + objs["ljs"].photon_n),
        ],
    ),
    "lj_muonN": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="lj_muonN"),
                   lambda objs: objs["ljs"].muon_n),
        ],
    ),
    # ljsource
    "ljsource_n": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="ljsource_n"),
                   lambda objs: ak.num(objs["ljsources"])),
        ],
    ),
    "ljsource_pt": h.Histogram(
        [
            h.Axis(common_axes["lj_pt"],
                   lambda objs: objs["ljsources"].p4.pt),
        ],
    ),
    "ljsource_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="ljsource_eta"),
                   lambda objs: objs["ljsources"].p4.eta),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="ljsource_phi"),
                   lambda objs: objs["ljsources"].p4.phi),
        ],
    ),
    "ljsource_charge": h.Histogram(
        [
            h.Axis(hist.axis.Integer(-1, 1, name="ljsource_charge"),
                   lambda objs: objs["ljsources"].charge),
        ],
    ),
    "ljsource_type": h.Histogram(
        [
            h.Axis(hist.axis.IntCategory([2, 3, 4, 8], name="lj_type"),
                   lambda objs: objs["ljsources"]["type"]), # avoid ak.Array.type
        ],
    ),
    # pfelectron-lj
    "electron_lj_dR": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="electron_lj_dR"),
                   lambda objs: dR(objs["electrons"].p4, objs["ljs"].p4))
        ],
    ),
    "electron_lj_dR_lowRange": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="electron_lj_dR_lowRange"),
                   lambda objs: dR(objs["electrons"].p4, objs["ljs"].p4))
        ],
    ),
    # pfphoton-lj
    "photon_lj_dR": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="photon_lj_dR"),
                   lambda objs: dR(objs["photons"].p4, objs["ljs"].p4))
        ],
    ),
    "photon_lj_dR_lowRange": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="photon_lj_dR_lowRange"),
                   lambda objs: dR(objs["photons"].p4, objs["ljs"].p4))
        ],
    ),
    "photon_lj_dR_reallyLowRange": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 0.1, name="photon_lj_dR_reallyLowRange"),
                   lambda objs: dR(objs["photons"].p4, objs["ljs"].p4))
        ],
    ),
    # pfmuon-lj
    "muon_lj_dR": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="muon_lj_dR"),
                   lambda objs: dR(objs["muons"].p4, objs["ljs"].p4))
        ],
    ),
    "muon_lj_dR_lowRange": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="muon_lj_dR_lowRange"),
                   lambda objs: dR(objs["muons"].p4, objs["ljs"].p4))
        ],
    ),
    # dsamuon-lj
    "dsaMuon_lj_dR": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="dsaMuon_lj_dR"),
                   lambda objs: dR(objs["dsaMuons"].p4, objs["ljs"].p4))
        ],
    ),
    "dsaMuon_lj_dR_lowRange": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="dsaMuon_lj_dR_lowRange"),
                   lambda objs: dR(objs["dsaMuons"].p4, objs["ljs"].p4))
        ],
    ),
    # lj-lj
    # fixme: these assume exactly two LJ per event
    "lj_lj_absdphi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="ljlj_absdphi"),
                   lambda objs: abs(objs["ljs"][:, 1].p4.phi - objs["ljs"][:, 0].p4.phi)),
        ],
    ),
    "lj_lj_invmass": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 2000, name="ljlj_mass"),
                   lambda objs: objs["ljs"].p4.sum().mass),
        ],
    ),
    "lj_lj_invmass_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name="ljlj_mass"),
                   lambda objs: objs["ljs"].p4.sum().mass),
        ],
    ),
    # gen
    "gen_abspid": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 40, name="gen_abspid"),
                   lambda objs: abs(objs["gens"].pid)),
        ],
    ),
    # genelectron
    "genE_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="genE_pt"),
                   lambda objs: abs(objs["genEs"].p4.pt)),
        ],
    ),
    # genelectron-genelectron
    "genE_genE_dR": h.Histogram(
        [
            # dR(subleading gen E, leading gen E) # fixme: assumes two gen electrons
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="genE_genE_dR"),
                   lambda objs: objs["genEs"][:, 1].p4.delta_r(objs["genEs"][:, 0].p4)),
        ],
    ),
    "genE_genE_pt": h.Histogram(
        [
            # fixme: assumes two gen electrons
            h.Axis(hist.axis.Regular(100, 0, 200, name="genE_genE_pt"),
                   lambda objs: objs["genEs"][:, 1].p4.add(objs["genEs"][:, 0].p4).pt),
        ],
    ),
    # genmuon
    "genMu_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="genMu_pt"),
                   lambda objs: abs(objs["genMus"].p4.pt)),
        ],
    ),
    # genmuon-genmuon
    "genMu_genMu_dR": h.Histogram(
        [
            # dR(subleading gen Mu, leading gen Mu) # fixme: assumes two gen muons
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="genMu_genMu_dR"),
                   lambda objs: objs["genMus"][:, 1].p4.delta_r(objs["genMus"][:, 0].p4)),
        ],
    ),
    "genMu_genMu_pt": h.Histogram(
        [
            # fixme: assumes two gen muons
            h.Axis(hist.axis.Regular(100, 0, 200, name="genMu_genMu_pt"),
                   lambda objs: objs["genMus"][:, 1].p4.add(objs["genMus"][:, 0].p4).pt),
        ],
    ),
    # gen dark photons (A)
    "genA_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="genA_pt"),
                   lambda objs: abs(objs["genAs"].p4.pt)),
        ],
    ),
    "genA_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="genA_eta"),
                   lambda objs: objs["genAs"].p4.eta),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="genA_phi"),
                   lambda objs: objs["genAs"].p4.phi),
        ],
    ),
    # genA-genA
    "genA_genA_dphi": h.Histogram( # fixme: assumes exactly two genA per event
        [
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="genA_genA_dphi"),
                   lambda objs: abs(objs["genAs"][:, 1].p4.phi - objs["genAs"][:, 0].p4.phi)),
        ],
    ),
    # genA-LJ
    "genA_lj_dR": h.Histogram(
        [
            # dR(A, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="genA_lj_dR"),
                   lambda objs: dR(objs["genAs"].p4, objs["ljs"].p4))
        ],
    ),
    "genA_lj_dR_lowRange": h.Histogram(
        [
            # dR(A, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="genA_lj_dR_lowRange"),
                   lambda objs: dR(objs["genAs"].p4, objs["ljs"].p4))
        ],
    ),
    "lj_genA_ptRatio": h.Histogram(
        [
            # (LJ pT)/(nearest A pT)
            h.Axis(hist.axis.Regular(50, 0, 2.0, name="lj_genA_ptRatio"),
                   lambda objs: objs["ljs"].p4.pt
                       / objs["ljs"].p4.nearest(objs["genAs"].p4).pt),
        ],
    ),
    "egm_lj_genA_ptRatio": h.Histogram(
        [
            # (LJ pT)/(nearest A pT)
            h.Axis(hist.axis.Regular(50, 0, 2.0, name="egm_lj_genA_ptRatio"),
                   lambda objs: obj_defs["egm_ljs"](objs).p4.pt
                       / obj_defs["egm_ljs"](objs).p4.nearest(objs["genAs"].p4).pt),
        ],
    ),
    "mu_lj_genA_ptRatio": h.Histogram(
        [
            # (LJ pT)/(nearest A pT)
            h.Axis(hist.axis.Regular(50, 0, 2.0, name="mu_lj_genA_ptRatio"),
                   lambda objs: obj_defs["mu_ljs"](objs).p4.pt
                       / obj_defs["mu_ljs"](objs).p4.nearest(objs["genAs"].p4).pt),
        ],
    ),
}
