"""Module to define the base SIDM processor"""

# python
import math
import yaml
import copy
# columnar analysis
from coffea import processor
import awkward as ak
#local
from analysis.tools import selection, cutflow, histogram, utilities
from analysis.definitions.hists import hist_definitions
# always reload local modules to pick up changes during development
import importlib
importlib.reload(selection)
importlib.reload(cutflow)
importlib.reload(histogram)
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

    def __init__(
        self,
        channel_names,
        hist_collection_names,
        selections_cfg="../configs/selections.yaml",
        histograms_cfg="../configs/hist_collections.yaml"
    ):
        self.channel_names = channel_names
        self.hist_collection_names = hist_collection_names
        self.selections_cfg = selections_cfg
        self.histograms_cfg = histograms_cfg

    def process(self, events):
        """Apply selections, make histograms and cutflow"""

        # handle metadata
        sample = events.metadata["sample"]

        # pt order objects
        events.ljsource = events.ljsource[ak.argsort(events.ljsource.p4.pt, ascending=False)]

        # evaluate selections for all analysis channels
        channels = self.build_analysis_channels(events)

        # define histograms
        hists = self.build_histograms()

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
            objs = {
                "ch" : channel,
                "pvs" : sel_pvs,
                "ljs" : sel_ljs,
            }

            # get arrays of event weights to apply to objects when filling object-level hists
            evt_weights = events.weightProduct[channel.all_evt_cuts.all(*channel.evt_cuts)]
            pv_weights = evt_weights*ak.ones_like(sel_pvs.z)
            lj_weights = evt_weights*ak.ones_like(sel_ljs.p4.pt)
            wgts = {
                "evt" : evt_weights,
                "pv" : pv_weights,
                "lj" : lj_weights,
            }

            # fill all hists
            for h in hists.values():
                h.fill(objs, wgts)

        out = {
            "cutflow" : cutflows,
            "hists" : {n : h.hist for n, h in hists.items()}, # output hist.Hists, not Histograms
        }
        return {sample : out}

    def build_analysis_channels(self, events):
        """Create list of Selection objects that define analysis channels"""
        with open(self.selections_cfg) as sel_cfg:
            selection_menu = yaml.safe_load(sel_cfg)

        channels = []
        for name in self.channel_names:
            cuts = selection_menu[name]
            # flatten object and event cut lists
            for obj_cuts in cuts["obj_cuts"].items():
                obj_cuts = utilities.flatten(obj_cuts)
            cuts["evt_cuts"] = utilities.flatten(cuts["evt_cuts"])

            channels.append(selection.Selection(name, cuts, events))
        return channels

    def build_histograms(self):
        """Create dictionary of Histogram objects"""
        with open(self.histograms_cfg) as hist_cfg:
            hist_menu = yaml.safe_load(hist_cfg)

        # build dictionary and create hist.Hist objects
        hists = {}
        for collection in self.hist_collection_names:
            hist_list = utilities.flatten(hist_menu[collection])
            for hist_name in hist_list:
                hists[hist_name] = copy.deepcopy(hist_definitions[hist_name])
                hists[hist_name].make_hist(self.channel_names)
        return hists

    def postprocess(self, accumulator):
        raise NotImplementedError
