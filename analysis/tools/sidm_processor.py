"""Module to define the base SIDM processor"""

# python
import math
import yaml
# columnar analysis
from coffea import processor
from coffea.analysis_tools import PackedSelection
import hist
import awkward as ak
#local
from analysis.tools import selection, cutflow, utilities
# always reload local modules to pick up changes during development
import importlib
importlib.reload(selection)
importlib.reload(cutflow)
importlib.reload(utilities)


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

    selections_cfg = "../configs/selections.yaml" # fixme: there must be a better way

    def __init__(self):
        """No need to do anything here -- everything happens in process"""

    def process(self, events, channel_names):
        """Apply selections, make histograms and cutflow"""

        # handle metadata
        sample = events.metadata["sample"]

        # pt order objects
        events.ljsource = events.ljsource[ak.argsort(events.ljsource.p4.pt, ascending=False)]

        channels = self.build_analysis_channels(events, channel_names)

        # define hists
        channel_axis = hist.axis.StrCategory([c.name for c in channels], name="channel")
        lj_pt_axis = hist.axis.Regular(100, 0, 100, name="lj_pt", label="Lepton jet pT [GeV]")
        lj_type_axis = hist.axis.IntCategory([2, 3, 4, 8], name="lj_type")
        hists = {
            # pv
            "pv_n" : hist.Hist(
                channel_axis,
                hist.axis.Integer(0, 100, name="pv_n"),
                storage="weight",
            ),
            "pv_ndof" : hist.Hist(
                channel_axis,
                hist.axis.Integer(0, 20, name="pv_ndof"),
                storage="weight",
            ),
            "pv_z" : hist.Hist(
                channel_axis,
                hist.axis.Regular(100, -50, 50, name="pv_z"),
                storage="weight",
            ),
            "pv_rho" : hist.Hist(
                channel_axis,
                hist.axis.Regular(100, -0.5, 0.5, name="pv_rho"),
                storage="weight",
            ),
            # lj
            "lj_n" : hist.Hist(
                channel_axis,
                hist.axis.Integer(0, 10, name="lj_n"),
                storage="weight",
            ),
            "lj_charge" : hist.Hist(
                channel_axis,
                hist.axis.Integer(-5, 5, name="lj_charge"),
                storage="weight",
            ),
            "lj_pt_type" : hist.Hist(
                channel_axis,
                lj_pt_axis,
                lj_type_axis,
                storage="weight",
            ),
            "lj_0_pt" : hist.Hist(
                channel_axis,
                lj_pt_axis,
                storage="weight",
            ),
            "lj_1_pt" : hist.Hist(
                channel_axis,
                lj_pt_axis,
                storage="weight",
            ),
            "lj_eta" : hist.Hist(
                channel_axis,
                hist.axis.Regular(50, -3, 3, name="lj_eta"),
                storage="weight",
            ),
            "lj_phi" : hist.Hist(
                channel_axis,
                hist.axis.Regular(50, -1*math.pi, math.pi, name="lj_phi"),
                storage="weight",
            ),
            # di-lj
            "lj_lj_absdphi" : hist.Hist(
                channel_axis,
                hist.axis.Regular(50, 0, 2*math.pi, name="ljlj_absdphi"),
                storage="weight",
            ),
            "lj_lj_invmass" : hist.Hist(
                channel_axis,
                hist.axis.Regular(200, 0, 2000, name="ljlj_mass"),
                storage="weight",
            ),
        }

        cutflows = {}
        for channel in channels:
            # apply full selection
            cutflows[channel.name] = cutflow.Cutflow(channel.all_evt_cuts, channel.evt_cuts, events.weightProduct)
            # update object collections to include only selected objects
            sel_pvs = events.pv[channel.obj_masks["pv"]]
            sel_ljs = events.ljsource[channel.obj_masks["ljsource"]]
            sel_ljs = sel_ljs[:, :2] # fixme: temporary hacky solution to only keep leading 2 LJs
            sel_pvs = sel_pvs[channel.all_evt_cuts.all(*channel.evt_cuts)]
            sel_ljs = sel_ljs[channel.all_evt_cuts.all(*channel.evt_cuts)]

            # get arrays of event weights to apply to objects when filling object-level hists
            evt_weights = events.weightProduct[channel.all_evt_cuts.all(*channel.evt_cuts)]
            pv_weights = evt_weights*ak.ones_like(sel_pvs.z)
            lj_weights = evt_weights*ak.ones_like(sel_ljs.p4.pt)

            # fill hists
            # pv
            hists["pv_n"].fill(
                channel=channel.name,
                pv_n=ak.num(sel_pvs),
                weight=evt_weights,
            )
            hists["pv_ndof"].fill(
                channel=channel.name,
                pv_ndof=ak.flatten(sel_pvs.ndof),
                weight=ak.flatten(pv_weights),
            )
            hists["pv_z"].fill(
                channel=channel.name,
                pv_z=ak.flatten(sel_pvs.z),
                weight=ak.flatten(pv_weights),
            )
            hists["pv_rho"].fill(
                channel=channel.name,
                pv_rho=ak.flatten(sel_pvs.rho),
                weight=ak.flatten(pv_weights),
            )
            # lj
            hists["lj_n"].fill(
                channel=channel.name,
                lj_n=ak.num(sel_ljs),
                weight=evt_weights,
            )
            hists["lj_charge"].fill(
                channel=channel.name,
                lj_charge=ak.flatten(sel_ljs.charge),
                weight=ak.flatten(lj_weights),
            )
            hists["lj_pt_type"].fill(
                channel=channel.name,
                lj_pt=ak.flatten(sel_ljs.p4.pt),
                lj_type=ak.flatten(sel_ljs['type']),
                weight=ak.flatten(lj_weights),
            )
            hists["lj_0_pt"].fill(
                channel=channel.name,
                lj_pt=sel_ljs[ak.num(sel_ljs) > 0, 0].p4.pt,
                weight=evt_weights,
            )
            hists["lj_1_pt"].fill(
                channel=channel.name,
                lj_pt=sel_ljs[ak.num(sel_ljs) > 1, 1].p4.pt,
                weight=evt_weights,
            )
            hists["lj_eta"].fill(
                channel=channel.name,
                lj_eta=ak.flatten(sel_ljs.p4.eta),
                weight=ak.flatten(lj_weights),
            )
            hists["lj_phi"].fill(
                channel=channel.name,
                lj_phi=ak.flatten(sel_ljs.p4.phi),
                weight=ak.flatten(lj_weights),
            )
            # di-lj
            hists["lj_lj_absdphi"].fill(
                channel=channel.name,
                ljlj_absdphi=abs(sel_ljs[:, 1].p4.phi - sel_ljs[:, 0].p4.phi),
                weight=evt_weights,
            )
            hists["lj_lj_invmass"].fill(
                channel=channel.name,
                ljlj_mass=sel_ljs.p4.sum().mass,
                weight=evt_weights,
            )

        out = {
            "cutflow" : cutflows,
            "hists" : hists,
        }
        return {sample : out}

    def build_analysis_channels(self, events, channel_names):
        """Create list of Selection objects that define analysis channels"""
        with open(type(self).selections_cfg) as sel_cfg:
            selection_menu = yaml.safe_load(sel_cfg)

        channels = []
        for name in channel_names:
            cuts = selection_menu[name]
            # flatten object and event cut lists
            for obj_cuts in cuts["obj_cuts"].items():
                obj_cuts = utilities.flatten(obj_cuts)
            cuts["evt_cuts"] = utilities.flatten(cuts["evt_cuts"])

            channels.append(selection.Selection(name, cuts, events))
        return channels

    def postprocess(self, accumulator):
        raise NotImplementedError
