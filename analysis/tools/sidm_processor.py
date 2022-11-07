"""Module to define the base SIDM processor"""

# python
import yaml
import copy
# columnar analysis
from coffea import processor
import awkward as ak
#local
from analysis.tools import selection, cutflow, histogram, utilities
from analysis.definitions.hists import hist_defs
# always reload local modules to pick up changes during development
import importlib
importlib.reload(selection)
importlib.reload(cutflow)
importlib.reload(histogram)
importlib.reload(utilities)


class SidmProcessor(processor.ProcessorABC):
    """Class to apply selections, make histograms, and make cutflows

    Accepts NanoEvents records that are assumed to have been produced by FFSchema. Selections are
    chosen by supplying a list of selection names (as defined in selections.yaml), and histograms
    are chosen by providing a list of histogram collection names (as definined in
    hist_collections.yaml).
    """

    def __init__(
        self,
        channel_names,
        hist_collection_names,
        selections_cfg="../configs/selections.yaml", # fixme: relative path could be bad idea
        histograms_cfg="../configs/hist_collections.yaml" # fixme: relative path could be bad idea
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
        # fixme: do this for all objects with a p4.pt attribute
        events.ljsource = events.ljsource[ak.argsort(events.ljsource.p4.pt, ascending=False)]

        # define objects
        objs = {
            "cosmicveto" : events.cosmicveto,
            "pvs" : events.pv,
            "ljs" : events.ljsource,
        }

        # evaluate object selections for all analysis channels
        channels = self.build_analysis_channels(objs)

        # define histograms
        hists = self.build_histograms()

        cutflows = {}
        for channel in channels:
            # apply object selection
            sel_objs = channel.apply_obj_masks(objs)
            sel_objs["ljs"] = sel_objs["ljs"][:, :2] # fixme: hacky way to only keep leading 2 LJs
            # apply event selection
            sel_objs = channel.apply_evt_cuts(sel_objs)

            # get arrays of event weights to apply to objects when filling object-level hists
            evt_weights = events.weightProduct[channel.all_evt_cuts.all(*channel.evt_cuts)]
            pv_weights = evt_weights*ak.ones_like(sel_objs["pvs"].z)
            lj_weights = evt_weights*ak.ones_like(sel_objs["ljs"].p4.pt)
            wgts = {
                "evt" : evt_weights,
                "pv" : pv_weights,
                "lj" : lj_weights,
            }

            # fill all hists
            sel_objs["ch"] = channel.name
            for h in hists.values():
                h.fill(sel_objs, wgts)

            # make cutflow
            cutflows[channel.name] = cutflow.Cutflow(channel.all_evt_cuts, channel.evt_cuts,
                                                     events.weightProduct)

        out = {
            "cutflow" : cutflows,
            "hists" : {n : h.hist for n, h in hists.items()}, # output hist.Hists, not Histograms
        }
        return {sample : out}

    def build_analysis_channels(self, objs):
        """Create list of Selection objects that define analysis channels"""
        with open(self.selections_cfg, encoding="utf8") as sel_cfg:
            selection_menu = yaml.safe_load(sel_cfg)

        channels = []
        for name in self.channel_names:
            cuts = selection_menu[name]
            # flatten object and event cut lists
            for obj_cuts in cuts["obj_cuts"].items():
                obj_cuts = utilities.flatten(obj_cuts)
            cuts["evt_cuts"] = utilities.flatten(cuts["evt_cuts"])

            # build Selection objects
            channels.append(selection.Selection(name, cuts, objs))
        return channels

    def build_histograms(self):
        """Create dictionary of Histogram objects"""
        with open(self.histograms_cfg, encoding="utf8") as hist_cfg:
            hist_menu = yaml.safe_load(hist_cfg)

        # build dictionary and create hist.Hist objects
        hists = {}
        for collection in self.hist_collection_names:
            collection = utilities.flatten(hist_menu[collection])
            for hist_name in collection:
                hists[hist_name] = copy.deepcopy(hist_defs[hist_name])
                hists[hist_name].make_hist(self.channel_names)
        return hists

    def postprocess(self, accumulator):
        raise NotImplementedError
