"""Module to define the base SIDM processor"""

# python
import copy
import importlib
# columnar analysis
from coffea import processor
from coffea.nanoevents.schemas.base import zip_forms
from coffea.nanoevents.methods import vector as cvec

import awkward as ak
import fastjet
import vector
#local
from sidm.tools import selection, cutflow, histogram, utilities
from sidm.definitions.hists import hist_defs
from sidm.definitions.objects import primary_objs
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
        lj_reco_choices=[0],
        selections_cfg="../configs/selections.yaml",
        histograms_cfg="../configs/hist_collections.yaml"
    ):
        self.channel_names = channel_names
        self.hist_collection_names = hist_collection_names
        self.lj_reco_choices = lj_reco_choices
        self.selections_cfg = selections_cfg
        self.histograms_cfg = histograms_cfg

    def process(self, events):
        """Apply selections, make histograms and cutflow"""

        # create object collections
        objs = {}
        for obj_name, obj_def in primary_objs.items():
            objs[obj_name] = obj_def(events)
            # pt order
            objs[obj_name] = self.order(objs[obj_name])

        cutflows = {}        
            
        # define histograms
        hists = self.build_histograms()

        # evaluate object selections for all analysis channels
        channels = self.build_analysis_channels(objs)

        # loop through lj reco choices and selections, treating each lj+selection pair as a unique
        # analysis channel
        for channel in channels:
            
            for lj_reco in self.lj_reco_choices:
                
                #>>>>>>>>>>>>> Would be better if this (applying object selection) was only done once per channel! 
                # But I ran into issues because adding "ch" and "lj_reco" to sel_objs made it have the 
                # wrong shape for apply_evt_cuts <<<<<<<<<<<<<<<<<<<<<<<
                # apply object selection
                
                sel_objs = channel.apply_obj_masks(objs)

                # reconstruct lepton jets
                sel_objs["ljs"] = self.build_lepton_jets(sel_objs, int(lj_reco))

                #apply lj selection here <<<<<<<<<<<<<<<<<<<<<<<
                
                # apply event selection
                sel_objs = channel.apply_evt_cuts(sel_objs)

                # fill all hists
                sel_objs["ch"] = channel.name
                sel_objs["lj_reco"] = lj_reco
                evt_weights = events.weightProduct[channel.all_evt_cuts.all(*channel.evt_cuts)]

                for h in hists.values():
                    h.fill(sel_objs, evt_weights)
                
                 # make cutflow
                if not lj_reco in cutflows:
                    cutflows[lj_reco] = {}
                cutflows[lj_reco][channel.name] = cutflow.Cutflow(
                    channel.all_evt_cuts, channel.evt_cuts, events.weightProduct)
                
        # lose lj_reco dimension to cutflows if only one reco was run
        if len(self.lj_reco_choices) == 1:
            cutflows = cutflows[self.lj_reco_choices[0]]
            
        out = {
            "cutflow": cutflows,
            "hists": {n: h.hist for n, h in hists.items()}, # output hist.Hists, not Histograms
        }
        return {events.metadata["dataset"]: out}
    
    def build_lepton_jets(self, objs, lj_reco):
        """Reconstruct lepton jets according to defintion given by lj_reco"""
        # fixme: define LJ reco choices in separate config instead of if, elif structure
        # fixme: what other LJ info do we want to store other than 4-vector? What's possible?
        
        # take lepton jets directly from ntuples
        if lj_reco == 0:
            ljs = objs["ljs"]

        # reconstruct lepton jets from scratch
        elif lj_reco == -1 or lj_reco > 0:
            
            if lj_reco == -1: #Use ljsource collection
                jet_def = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4) 
                lj_inputs = vector.zip(
                    {
                        "type": objs["ljsources"]["type"],
                        "charge": objs["ljsources"].charge,
                        "px": objs["ljsources"].px,
                        "py": objs["ljsources"].py,
                        "pz": objs["ljsources"].pz,
                        "E":  objs["ljsources"].energy,
                    },
                )
                
            else: #Use electron/muon/photon/dsamuon collections with a custom distance parameter           
                
                distance_param = lj_reco/10.0
                jet_def = fastjet.JetDefinition(fastjet.antikt_algorithm,distance_param) 
                muon_inputs = vector.zip(
                    {
                        "type": 3,
                        "charge": objs["muons"].charge,
                        "px": objs["muons"].px,
                        "py": objs["muons"].py,
                        "pz": objs["muons"].pz,
                        "E":  objs["muons"].energy,
                    },
                )

                dsamuon_inputs = vector.zip(
                    {
                        "type": 8,
                        "charge": objs["dsaMuons"].charge,
                        "px": objs["dsaMuons"].px,
                        "py": objs["dsaMuons"].py,
                        "pz": objs["dsaMuons"].pz,
                        "E":  objs["dsaMuons"].energy,
                    },
                )

                ele_inputs = vector.zip(
                    {
                        "type": 2,
                        "charge": objs["electrons"].charge,
                        "px": objs["electrons"].px,
                        "py": objs["electrons"].py,
                        "pz": objs["electrons"].pz,
                        "E":  objs["electrons"].energy,
                    },
                )

                photon_inputs = vector.zip(
                    {
                        "type": 4,
                        "charge": ak.without_parameters(ak.zeros_like(objs["photons"].px), behavior={}),
                        "px": objs["photons"].px,
                        "py": objs["photons"].py,
                        "pz": objs["photons"].pz,
                        "E":  objs["photons"].energy,
                    },
                )

                lj_inputs = ak.concatenate([muon_inputs, dsamuon_inputs, ele_inputs, photon_inputs],axis=-1)

            cluster = fastjet.ClusterSequence(lj_inputs, jet_def)
            jets = cluster.inclusive_jets()

            # turn lepton jets back into LorentzVectors that match existing structures
            ljs = ak.zip(
                {"x": jets.x,
                 "y": jets.y,
                 "z": jets.z,
                  "t": jets.t},
                with_name="LorentzVector",
                behavior=cvec.behavior
            )
    
            # add jet constituent info
            ljs["constituents"] = cluster.constituents()
            ### >>>>>>>>> Dummy placeholders! Replace with an actual value at some point.
            ljs["dRSpread"]= ak.num(ljs.constituents[ljs.constituents["type"] == 2])
            ljs["pfIsolation05"]= ak.num(ljs.constituents[ljs.constituents["type"] == 2])
            ljs["pfIsolationPtNoPU05"]= ak.num(ljs.constituents[ljs.constituents["type"] == 2])
            #### <<<<<<<<<
            ljs["muon_n"] = ak.num(ljs.constituents[(ljs.constituents["type"] == 3)
                                                    | (ljs.constituents["type"] == 8)],axis=-1)
            ljs["electron_n"] = ak.num(ljs.constituents[ljs.constituents["type"] == 2],axis=-1)
            ljs["photon_n"] = ak.num(ljs.constituents[ljs.constituents["type"] == 4],axis=-1)
 
            # apply cuts to match cuts applied in ntuples
            # fixme: implement num(mu) % 2 == 0 and other cuts
            ljs = ljs[ljs.pt>30]
            #ljs = ljs[ljs.p4.pt>30]
            # fixme: calculate muon_n and other fields
            
        else:
            raise NotImplementedError(f"{lj_reco} is not a recognized LJ reconstruction choice")
        # pt order and return
        return self.order(ljs)

    def build_analysis_channels(self, objs):
        """Create list of Selection objects that define analysis channels"""
        selection_menu = utilities.load_yaml(self.selections_cfg)

        channels = []
        evaluated_obj_cuts = {}
        for name in self.channel_names:
            cuts = selection_menu[name]
            # flatten object and event cut lists
            for obj, obj_cuts in cuts["obj_cuts"].items():
                cuts["obj_cuts"][obj] = utilities.flatten(obj_cuts)
            cuts["evt_cuts"] = utilities.flatten(cuts["evt_cuts"])

            # build Selection objects
            channel = selection.Selection(name, cuts)
            evaluated_obj_cuts = channel.evaluate_obj_cuts(objs, evaluated_obj_cuts)
            channel.make_obj_masks(evaluated_obj_cuts)
            channels.append(channel)

        return channels

    def build_histograms(self):
        """Create dictionary of Histogram objects"""
        hist_menu = utilities.load_yaml(self.histograms_cfg)

        # build dictionary and create hist.Hist objects
        hists = {}
        for collection in self.hist_collection_names:
            collection = utilities.flatten(hist_menu[collection])
            for hist_name in collection:
                hists[hist_name] = copy.deepcopy(hist_defs[hist_name])
                # Add lj_reco axis only when more than one reco is run
                lj_reco_names = self.lj_reco_choices if len(self.lj_reco_choices) > 1 else None
                hists[hist_name].make_hist(self.channel_names, lj_reco_names)
        return hists

    def order(self, obj):
        """Explicitly order objects"""
        # pt order objects with a pt attribute
        if hasattr(obj, "p4"):
            obj = obj[ak.argsort(obj.pt, ascending=False)]
        # fixme: would be good to explicitly order other objects as well
        return obj

    def postprocess(self, accumulator):
        pass
