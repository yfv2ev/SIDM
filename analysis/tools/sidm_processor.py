"""Module to define the base SIDM processor"""

# python
import copy
# columnar analysis
from coffea import processor
from coffea.analysis_tools import PackedSelection
import hist
import awkward as ak
#local
from analysis.tools import cutflow
# always reload local modules to pick up changes during development
import importlib
importlib.reload(cutflow)


class SidmProcessor(processor.ProcessorABC):
    """Class to perform first processor tests.

    I expect the contents to evolve along the following lines:
        1. do a few basic tests
        2. implement a basic SIDM preselection
        3. make a few test histograms
        4. then convert this into the base SIDM processor from which all other SIDM
        processors will inherit (spinning the selection and histogram details off to other files
        in the process)
    """

    def __init__(self):
        """No need to do anything here -- everything happens in process"""

    def process(self, events):
        """Apply selections, make histograms and cutflow"""

        # handle metadata
        sample = events.metadata["sample"]

        # define PV selection
        pvs = events.pv
        pvs = pvs[
            (pvs.ndof > 4)
            & (abs(pvs.z) < 24) # assume cm
            & (abs(pvs.rho) < 0.2) # assume cm (fixme: weinan disagrees: https://github.com/phylsix/Firefighter/blob/master/recoStuff/python/ffPrimaryVertexFilter_cfi.py)
        ]

        # define LJ selection
        ljs = events.ljsource
        ljs = ljs[ak.argsort(ljs.p4.pt, ascending=False)] # pt ordering
        ljs = ljs[
            (ljs.p4.pt > 30) # GeV, not applied in ntuples
            & (abs(ljs.p4.eta) < 2.4) # already applied in ntuples
            # fixme: add muon number and charge constraints if not already applied in ntuples
        ]
        ljs = ljs[:, :2] # only consider leading two LJs; fixme: would be nice to not do this
        
        # define muon- and egamma-type LJs
        mu_ljs = ljs[(ljs["type"] == 3) | (ljs["type"] == 8)] # 3=pfmuon, 8=dsamuon
        egm_ljs = ljs[(ljs["type"] == 2) | (ljs["type"] == 4)] # 2=pfelectron, 4=pfphoton

        # define base selection
        all_cuts = PackedSelection()
        all_cuts.add("PV filter", ak.num(pvs) >= 1)
        all_cuts.add("Cosmic veto", events.cosmicveto.result)
        all_cuts.add(">=2 LJs", ak.num(ljs) >= 2)
        base_selection = all_cuts.names[:]
        # define analysis channels (assuming exactly 2 LJs in each event)
        all_cuts.add("4mu", ak.num(mu_ljs) == 2)
        all_cuts.add("2mu2e", (ak.num(mu_ljs) == 1) & (ak.num(egm_ljs) == 1))
        channels = {
            "4mu" : base_selection + ["4mu"],
            "2mu2e" : base_selection + ["2mu2e"],
        }

        # define hists
        channel_axis = hist.axis.StrCategory(channels.keys(), name="channel")
        lj_pt_axis = hist.axis.Regular(100, 0, 100, name="lj_pt", label="Lepton jet pT [GeV]")
        lj_type_axis = hist.axis.IntCategory([2, 3, 4, 8], name="lj_type")
        hists = {
            # pv
            "pv_n" : hist.Hist(channel_axis, hist.axis.Integer(0, 100, name="pv_n")),
            "pv_ndof" : hist.Hist(channel_axis, hist.axis.Integer(0, 20, name="pv_ndof")),
            "pv_z" : hist.Hist(channel_axis, hist.axis.Regular(100, -50, 50, name="pv_z")),
            "pv_rho" : hist.Hist(channel_axis, hist.axis.Regular(100, -0.5, 0.5, name="pv_rho")),
            # lj
            "lj_n" : hist.Hist(channel_axis, hist.axis.Integer(0, 10, name="lj_n")),
            "lj_charge" : hist.Hist(channel_axis, hist.axis.Integer(-5, 5, name="lj_charge")),
            "lj_pt_type" : hist.Hist(channel_axis, lj_pt_axis, lj_type_axis),
            "lj_0_pt" : hist.Hist(channel_axis, lj_pt_axis),
            "lj_1_pt" : hist.Hist(channel_axis, lj_pt_axis),
            "lj_eta" : hist.Hist(channel_axis, hist.axis.Regular(100, -3, 3, name="lj_eta")),
            "lj_phi" : hist.Hist(channel_axis, hist.axis.Regular(100, -3.14, 3.14, name="lj_phi")),
        }

        cutflows = {}
        for channel, selection in channels.items():
            # apply full selection
            cutflows[channel] = cutflow.Cutflow(all_cuts, selection)
            # update object collections to only include selected events
            sel_pvs = pvs[all_cuts.all(*selection)]
            sel_ljs = ljs[all_cuts.all(*selection)]

            # fill hists
            # pv
            hists["pv_n"].fill(channel=channel, pv_n=ak.num(pvs))
            hists["pv_ndof"].fill(channel=channel, pv_ndof=ak.flatten(pvs.ndof))
            hists["pv_z"].fill(channel=channel, pv_z=ak.flatten(pvs.z))
            hists["pv_rho"].fill(channel=channel, pv_rho=ak.flatten(pvs.rho))
            # lj
            hists["lj_n"].fill(channel=channel, lj_n=ak.num(sel_ljs))
            hists["lj_charge"].fill(channel=channel, lj_charge=ak.flatten(sel_ljs.charge))
            hists["lj_pt_type"].fill(
                channel=channel, lj_pt=ak.flatten(sel_ljs.p4.pt), lj_type=ak.flatten(sel_ljs['type'])
            )
            hists["lj_0_pt"].fill(channel=channel, lj_pt=sel_ljs[ak.num(sel_ljs) > 0, 0].p4.pt)
            hists["lj_1_pt"].fill(channel=channel, lj_pt=sel_ljs[ak.num(sel_ljs) > 1, 1].p4.pt)
            hists["lj_eta"].fill(channel=channel, lj_eta=ak.flatten(sel_ljs.p4.eta))
            hists["lj_phi"].fill(channel=channel, lj_phi=ak.flatten(sel_ljs.p4.phi))

        out = {
            "cutflow" : cutflows,
            "hists" : hists,
        }
        return {sample : out}

    def postprocess(self, accumulator):
        raise NotImplementedError
