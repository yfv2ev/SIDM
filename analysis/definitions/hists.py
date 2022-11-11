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
# always reload local modules to pick up changes during development
importlib.reload(h)


common_axes = {
    "lj_pt" : hist.axis.Regular(100, 0, 100, name="lj_pt", label="Lepton jet pT [GeV]")
}

hist_defs = {
    # pv
    "pv_n" : h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 100, name="pv_n"),
                   lambda objs: ak.num(objs["pvs"])),
        ],
        weight_key="evt",
    ),
    "pv_ndof" : h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0, 100, name="pv_ndof"),
                   lambda objs: ak.flatten(objs["pvs"].ndof)),
        ],
        weight_key="pv",
    ),
    "pv_z" : h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, -50, 50, name="pv_z"),
                   lambda objs: ak.flatten(objs["pvs"].z)),
        ],
        weight_key="pv",
    ),
    "pv_rho" : h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, -0.5, 0.5, name="pv_rho"),
                   lambda objs: ak.flatten(objs["pvs"].rho)),
        ],
        weight_key="pv",
    ),
    # lj
    "lj_n" : h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="lj_n"),
                   lambda objs: ak.num(objs["ljs"])),
        ],
        weight_key="evt",
    ),
    "lj_charge" : h.Histogram(
        [
            h.Axis(hist.axis.Integer(-5, 5, name="lj_charge"),
                   lambda objs: ak.flatten(objs["ljs"].charge)),
        ],
        weight_key="lj",
    ),
    "lj_pt_type" : h.Histogram(
        [
            h.Axis(common_axes["lj_pt"],
                   lambda objs: ak.flatten(objs["ljs"].p4.pt)),
            h.Axis(hist.axis.IntCategory([2, 3, 4, 8], name="lj_type"),
                   lambda objs: ak.flatten(objs["ljs"]["type"])), # avoid ak.Array.type
        ],
        weight_key="lj",
    ),
    "lj0_pt" : h.Histogram(
        [
            h.Axis(common_axes["lj_pt"],
                   lambda objs: objs["ljs"][ak.num(objs["ljs"]) > 0, 0].p4.pt),
        ],
        weight_key="evt",
    ),
    "lj1_pt" : h.Histogram(
        [
            h.Axis(common_axes["lj_pt"],
                   lambda objs: objs["ljs"][ak.num(objs["ljs"]) > 1, 1].p4.pt),
        ],
        weight_key="evt",
    ),
    "lj_eta_phi" : h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="lj_eta"),
                   lambda objs: ak.flatten(objs["ljs"].p4.eta)),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="lj_phi"),
                   lambda objs: ak.flatten(objs["ljs"].p4.phi)),
        ],
        weight_key="lj",
    ),
    # lj-lj
    # fixme: these assume exactly two LJ per event
    "lj_lj_absdphi" : h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="ljlj_absdphi"),
                   lambda objs: abs(objs["ljs"][:, 1].p4.phi - objs["ljs"][:, 0].p4.phi)),
        ],
        weight_key="evt"
    ),
    "lj_lj_invmass" : h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 2000, name="ljlj_mass"),
                   lambda objs: objs["ljs"].p4.sum().mass),
        ],
        weight_key="evt"
    ),
    "lj_lj_invmass_lowRange" : h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name="ljlj_mass"),
                   lambda objs: objs["ljs"].p4.sum().mass),
        ],
        weight_key="evt"
    ),
    # gen
    "gen_abspid" : h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 40, name="gen_abspid"),
                   lambda objs: ak.flatten(abs(objs["gens"].pid))),
        ],
        weight_key="gen"
    ),
    "genA_genA_dphi" : h.Histogram( # fixme: assumes exactly two genA per event
        [
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="genA_genA_dphi"),
                   lambda objs: abs(objs["genAs"][:, 1].p4.phi - objs["genAs"][:, 0].p4.phi)),
        ],
        weight_key="evt"
    ),
    # gen-LJ
    "genA_lj_dR" : h.Histogram(
        [
            # dR(A, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="gen_genA_lj_dR"),
                   lambda objs: ak.flatten(
                      objs["genAs"].p4.nearest(objs["ljs"].p4, return_metric=True)[1])),
        ],
        weight_key="genA"
    ),
    "genA_lj_dR_lowRange" : h.Histogram(
        [
            # dR(A, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="gen_genA_lj_dR_lowRange"),
                   lambda objs: ak.flatten(
                      objs["genAs"].p4.nearest(objs["ljs"].p4, return_metric=True)[1])),
        ],
        weight_key="genA"
    ),
}
