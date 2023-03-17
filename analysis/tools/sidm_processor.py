"""Module to define the base SIDM processor"""

# python
import copy
import importlib
import yaml
# columnar analysis
from coffea import processor
import awkward as ak
#local
from analysis.tools import selection, cutflow, histogram, utilities
from analysis.definitions.hists import hist_defs
# always reload local modules to pick up changes during development
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

        # define objects
        # fixme: this should be defined elsewhere
        objs = {
            "cosmicveto": events.cosmicveto,
            "pvs": events.pv,
            "electrons": events.electron,
            "photons": events.pfphoton,
            "muons": events.muon,
            "dsaMuons": events.dsamuon,
            "ljs": events.pfjet,
            "ljsources": events.ljsource,
            "gens": events.gen,
            "genEs": events.gen[abs(events.gen.pid) == 11],
            "genMus": events.gen[abs(events.gen.pid) == 13],
            "genAs": events.gen[events.gen.pid == 32],
        }

        # pt order objects
        for obj in objs:
            if hasattr(objs[obj], "p4"):
                objs[obj] = objs[obj][ak.argsort(objs[obj].p4.pt, ascending=False)]

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

            # fill all hists
            sel_objs["ch"] = channel.name
            evt_weights = events.weightProduct[channel.all_evt_cuts.all(*channel.evt_cuts)]
            for h in hists.values():
                h.fill(sel_objs, evt_weights)

            # make cutflow
            cutflows[channel.name] = cutflow.Cutflow(channel.all_evt_cuts, channel.evt_cuts,
                                                     events.weightProduct)

        out = {
            "cutflow": cutflows,
            "hists": {n: h.hist for n, h in hists.items()}, # output hist.Hists, not Histograms
        }
        return {events.metadata["dataset"]: out}

    def build_analysis_channels(self, objs):
        """Create list of Selection objects that define analysis channels"""
        with open(self.selections_cfg, encoding="utf8") as sel_cfg:
            selection_menu = yaml.safe_load(sel_cfg)

        channels = []
        for name in self.channel_names:
            cuts = selection_menu[name]
            # flatten object and event cut lists
            for obj, obj_cuts in cuts["obj_cuts"].items():
                cuts["obj_cuts"][obj] = utilities.flatten(obj_cuts)
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
        pass
