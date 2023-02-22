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
        weight_key="evt",
    ),
    "pv_ndof": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0, 100, name="pv_ndof"),
                   lambda objs: ak.flatten(objs["pvs"].ndof)),
        ],
        weight_key="pv",
    ),
    "pv_z": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, -50, 50, name="pv_z"),
                   lambda objs: ak.flatten(objs["pvs"].z)),
        ],
        weight_key="pv",
    ),
    "pv_rho": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, -0.5, 0.5, name="pv_rho"),
                   lambda objs: ak.flatten(objs["pvs"].rho)),
        ],
        weight_key="pv",
    ),
    # pfelectron
    "electron_n": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="electron_n"),
                   lambda objs: ak.num(objs["electrons"])),
        ],
        weight_key="evt",
    ),
    "electron_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="electron_pt"),
                   lambda objs: ak.flatten(objs["electrons"].p4.pt)),
        ],
        weight_key="electron",
    ),
    "electron_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="electron_eta"),
                   lambda objs: ak.flatten(objs["electrons"].p4.eta)),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="electron_phi"),
                   lambda objs: ak.flatten(objs["electrons"].p4.phi)),
        ],
        weight_key="electron",
    ),
    # pfphoton
    "photon_n": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="photon_n"),
                   lambda objs: ak.num(objs["photons"])),
        ],
        weight_key="evt",
    ),
    "photon_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="photon_pt"),
                   lambda objs: ak.flatten(objs["photons"].p4.pt)),
        ],
        weight_key="photon",
    ),
    "photon_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="photon_eta"),
                   lambda objs: ak.flatten(objs["photons"].p4.eta)),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="photon_phi"),
                   lambda objs: ak.flatten(objs["photons"].p4.phi)),
        ],
        weight_key="photon",
    ),
    # pfmuon
    "muon_n": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="muon_n"),
                   lambda objs: ak.num(objs["muons"])),
        ],
        weight_key="evt",
    ),
    "muon_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="muon_pt"),
                   lambda objs: ak.flatten(objs["muons"].p4.pt)),
        ],
        weight_key="muon",
    ),
    "muon_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="muon_eta"),
                   lambda objs: ak.flatten(objs["muons"].p4.eta)),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="muon_phi"),
                   lambda objs: ak.flatten(objs["muons"].p4.phi)),
        ],
        weight_key="muon",
    ),
    # dsamuon
    "dsaMuon_n": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="dsaMuon_n"),
                   lambda objs: ak.num(objs["dsaMuons"])),
        ],
        weight_key="evt",
    ),
    "dsaMuon_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="dsaMuon_pt"),
                   lambda objs: ak.flatten(objs["dsaMuons"].p4.pt)),
        ],
        weight_key="dsaMuon",
    ),
    "dsaMuon_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="dsaMuon_eta"),
                   lambda objs: ak.flatten(objs["dsaMuons"].p4.eta)),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="dsaMuon_phi"),
                   lambda objs: ak.flatten(objs["dsaMuons"].p4.phi)),
        ],
        weight_key="dsaMuon",
    ),
    # pfjet
    "pfjet_n": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="pfjet_n"),
                   lambda objs: ak.num(objs["pfjets"])),
        ],
        weight_key="evt",
    ),
    "pfjet_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="pfjet_pt"),
                   lambda objs: ak.flatten(objs["pfjets"].p4.pt)),
        ],
        weight_key="pfjet",
    ),
    "pfjet_electronN": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="pfjet_electronN"),
                   lambda objs: ak.flatten(objs["pfjets"].electron_n)),
        ],
        weight_key="pfjet",
    ),
    "pfjet_photonN": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="pfjet_photonN"),
                   lambda objs: ak.flatten(objs["pfjets"].photon_n)),
        ],
        weight_key="pfjet",
    ),
    "pfjet_electronPhotonN": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="pfjet_electronPhotonN"),
                   lambda objs: ak.flatten(objs["pfjets"].electron_n + objs["pfjets"].photon_n)),
        ],
        weight_key="pfjet",
    ),
    "pfjet_muonN": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="pfjet_muonN"),
                   lambda objs: ak.flatten(objs["pfjets"].muon_n)),
        ],
        weight_key="pfjet",
    ),
    # lj
    "lj_n": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="lj_n"),
                   lambda objs: ak.num(objs["ljs"])),
        ],
        weight_key="evt",
    ),
    "lj_charge": h.Histogram(
        [
            h.Axis(hist.axis.Integer(-5, 5, name="lj_charge"),
                   lambda objs: ak.flatten(objs["ljs"].charge)),
        ],
        weight_key="lj",
    ),
    "lj_pt_type": h.Histogram(
        [
            h.Axis(common_axes["lj_pt"],
                   lambda objs: ak.flatten(objs["ljs"].p4.pt)),
            h.Axis(hist.axis.IntCategory([2, 3, 4, 8], name="lj_type"),
                   lambda objs: ak.flatten(objs["ljs"]["type"])), # avoid ak.Array.type
        ],
        weight_key="lj",
    ),
    "lj0_pt": h.Histogram(
        [
            h.Axis(common_axes["lj_pt"],
                   lambda objs: objs["ljs"][ak.num(objs["ljs"]) > 0, 0].p4.pt),
        ],
        weight_key="evt",
    ),
    "lj1_pt": h.Histogram(
        [
            h.Axis(common_axes["lj_pt"],
                   lambda objs: objs["ljs"][ak.num(objs["ljs"]) > 1, 1].p4.pt),
        ],
        weight_key="evt",
    ),
    "lj_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="lj_eta"),
                   lambda objs: ak.flatten(objs["ljs"].p4.eta)),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="lj_phi"),
                   lambda objs: ak.flatten(objs["ljs"].p4.phi)),
        ],
        weight_key="lj",
    ),
    # pfelectron-lj
    "electron_lj_dR": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="electron_lj_dR"),
                   lambda objs: ak.flatten(dR(objs["electrons"].p4, objs["ljs"].p4)))
        ],
        weight_key="electron"
    ),
    "electron_lj_dR_lowRange": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="electron_lj_dR_lowRange"),
                   lambda objs: ak.flatten(dR(objs["electrons"].p4, objs["ljs"].p4)))
        ],
        weight_key="electron"
    ),
    # pfphoton-lj
    "photon_lj_dR": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="photon_lj_dR"),
                   lambda objs: ak.flatten(dR(objs["photons"].p4, objs["ljs"].p4)))
        ],
        weight_key="photon"
    ),
    "photon_lj_dR_lowRange": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="photon_lj_dR_lowRange"),
                   lambda objs: ak.flatten(dR(objs["photons"].p4, objs["ljs"].p4)))
        ],
        weight_key="photon"
    ),
    "photon_lj_dR_reallyLowRange": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 0.1, name="photon_lj_dR_reallyLowRange"),
                   lambda objs: ak.flatten(dR(objs["photons"].p4, objs["ljs"].p4)))
        ],
        weight_key="photon"
    ),
    # pfmuon-lj
    "muon_lj_dR": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="muon_lj_dR"),
                   lambda objs: ak.flatten(dR(objs["muons"].p4, objs["ljs"].p4)))
        ],
        weight_key="muon"
    ),
    "muon_lj_dR_lowRange": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="muon_lj_dR_lowRange"),
                   lambda objs: ak.flatten(dR(objs["muons"].p4, objs["ljs"].p4)))
        ],
        weight_key="muon"
    ),
    # dsamuon-lj
    "dsaMuon_lj_dR": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="dsaMuon_lj_dR"),
                   lambda objs: ak.flatten(dR(objs["dsaMuons"].p4, objs["ljs"].p4)))
        ],
        weight_key="dsaMuon"
    ),
    "dsaMuon_lj_dR_lowRange": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="dsaMuon_lj_dR_lowRange"),
                   lambda objs: ak.flatten(dR(objs["dsaMuons"].p4, objs["ljs"].p4)))
        ],
        weight_key="dsaMuon"
    ),
    # lj-lj
    # fixme: these assume exactly two LJ per event
    "lj_lj_absdphi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="ljlj_absdphi"),
                   lambda objs: abs(objs["ljs"][:, 1].p4.phi - objs["ljs"][:, 0].p4.phi)),
        ],
        weight_key="evt"
    ),
    "lj_lj_invmass": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 2000, name="ljlj_mass"),
                   lambda objs: objs["ljs"].p4.sum().mass),
        ],
        weight_key="evt"
    ),
    "lj_lj_invmass_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name="ljlj_mass"),
                   lambda objs: objs["ljs"].p4.sum().mass),
        ],
        weight_key="evt"
    ),
    # gen
    "gen_abspid": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 40, name="gen_abspid"),
                   lambda objs: ak.flatten(abs(objs["gens"].pid))),
        ],
        weight_key="gen"
    ),
    # genelectron
    "genE_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="genE_pt"),
                   lambda objs: ak.flatten(abs(objs["genEs"].p4.pt))),
        ],
        weight_key="genE"
    ),
    # genelectron-genelectron
    "genE_genE_dR": h.Histogram(
        [
            # dR(subleading gen E, leading gen E) # fixme: assumes two gen electrons
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="genE_genE_dR"),
                   lambda objs: objs["genEs"][:, 1].p4.delta_r(objs["genEs"][:, 0].p4)),
        ],
        weight_key="evt"
    ),
    "genE_genE_pt": h.Histogram(
        [
            # fixme: assumes two gen electrons
            h.Axis(hist.axis.Regular(100, 0, 200, name="genE_genE_pt"),
                   lambda objs: objs["genEs"][:, 1].p4.add(objs["genEs"][:, 0].p4).pt),
        ],
        weight_key="evt"
    ),
    # genmuon
    "genMu_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="genMu_pt"),
                   lambda objs: ak.flatten(abs(objs["genMus"].p4.pt))),
        ],
        weight_key="genMu"
    ),
    # genmuon-genmuon
    "genMu_genMu_dR": h.Histogram(
        [
            # dR(subleading gen Mu, leading gen Mu) # fixme: assumes two gen muons
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="genMu_genMu_dR"),
                   lambda objs: objs["genMus"][:, 1].p4.delta_r(objs["genMus"][:, 0].p4)),
        ],
        weight_key="evt"
    ),
    "genMu_genMu_pt": h.Histogram(
        [
            # fixme: assumes two gen muons
            h.Axis(hist.axis.Regular(100, 0, 200, name="genMu_genMu_pt"),
                   lambda objs: objs["genMus"][:, 1].p4.add(objs["genMus"][:, 0].p4).pt),
        ],
        weight_key="evt"
    ),
    # gen dark photons (A)
    "genA_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="genA_pt"),
                   lambda objs: ak.flatten(abs(objs["genAs"].p4.pt))),
        ],
        weight_key="genA"
    ),
    "genA_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="genA_eta"),
                   lambda objs: ak.flatten(objs["genAs"].p4.eta)),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="genA_phi"),
                   lambda objs: ak.flatten(objs["genAs"].p4.phi)),
        ],
        weight_key="genA"
    ),
    # genA-genA
    "genA_genA_dphi": h.Histogram( # fixme: assumes exactly two genA per event
        [
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="genA_genA_dphi"),
                   lambda objs: abs(objs["genAs"][:, 1].p4.phi - objs["genAs"][:, 0].p4.phi)),
        ],
        weight_key="evt"
    ),
    # genA-LJ
    "genA_lj_dR": h.Histogram(
        [
            # dR(A, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="genA_lj_dR"),
                   lambda objs: ak.flatten(dR(objs["genAs"].p4, objs["ljs"].p4)))
        ],
        weight_key="genA"
    ),
    "genA_lj_dR_lowRange": h.Histogram(
        [
            # dR(A, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="genA_lj_dR_lowRange"),
                   lambda objs: ak.flatten(dR(objs["genAs"].p4, objs["ljs"].p4)))
        ],
        weight_key="genA"
    ),
    "lj_genA_ptRatio_lj_type": h.Histogram(
        [
            # (LJ pT)/(nearest A pT)
            h.Axis(hist.axis.Regular(50, 0, 2.0, name="lj_genA_ptRatio"),
                   lambda objs: ak.flatten((objs["ljs"].p4.pt
                       / objs["ljs"].p4.nearest(objs["genAs"].p4).pt))),
            h.Axis(hist.axis.IntCategory([2, 3, 4, 8], name="lj_type"),
                   lambda objs: ak.flatten(objs["ljs"]["type"])), # avoid ak.Array.type
        ],
        weight_key="lj"
    ),
}
