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
from sidm.tools import selection, jagged_selection, cutflow, histogram, utilities
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
        
        ### Make list of all object-level cuts; define object-level, lj-level, and event-level cuts per channel
        selection_menu = utilities.load_yaml(self.selections_cfg)
        
        all_obj_cuts = {}
        channel_cuts = {}
        
        #fixme: could make this the new build_analysis_channels function
        for channel in self.channel_names:
            channel_cuts[channel] = {}
            channel_cuts[channel]["obj"] = {}
            channel_cuts[channel]["evt"] = {}
            channel_cuts[channel]["lj"] = {}
            channel_cuts[channel]["lj"]["ljs"] = {}
     
            #Merge all object level cuts into a single list to be evaluated once
            cuts = selection_menu[channel]
            for obj, obj_cuts in cuts["obj_cuts"].items():
                
                if obj not in all_obj_cuts:
                    all_obj_cuts[obj] = []
                all_obj_cuts[obj] = utilities.add_unique_and_flatten(all_obj_cuts[obj],obj_cuts)
                
                if obj not in channel_cuts[channel]["obj"]:
                    channel_cuts[channel]["obj"][obj] = []
                channel_cuts[channel]["obj"][obj] = utilities.flatten(obj_cuts)
                
            channel_cuts[channel]["evt"] = utilities.flatten(cuts["evt_cuts"])
            channel_cuts[channel]["lj"]["ljs"] = utilities.flatten(cuts["lj_cuts"])
                
        #Evaluate all object-level cuts
        obj_selection = jagged_selection.JaggedSelection(all_obj_cuts)
        obj_selection.evaluate_obj_cuts(objs)
        
        # loop through lj reco choices and channels, treating each lj+channel pair as a unique Selection
        for channel in self.channel_names:
            
            # apply object selection
            channel_objs = obj_selection.make_and_apply_obj_masks(objs, channel_cuts[channel]["obj"])
              
            for lj_reco in self.lj_reco_choices: 
                
                sel_objs = channel_objs

                # reconstruct lepton jets
                sel_objs["ljs"] = self.build_lepton_jets(channel_objs, int(lj_reco))

                #apply lj selection
                lj_selection = jagged_selection.JaggedSelection(channel_cuts[channel]["lj"])
                lj_selection.evaluate_obj_cuts(sel_objs)
                sel_objs = lj_selection.make_and_apply_obj_masks(sel_objs,channel_cuts[channel]["lj"])
                                   
                # build Selection objects and apply event selection
                evt_selection = selection.Selection(channel_cuts[channel]["evt"])
                sel_objs = evt_selection.apply_evt_cuts(sel_objs)

                # fill all hists
                sel_objs["ch"] = channel
                sel_objs["lj_reco"] = lj_reco
            ## fixme: add option to make evt_weights all 1
               # evt_weights =  ak.ones_like(events.weightProduct[channel.all_evt_cuts.all(*channel.evt_cuts)])
                evt_weights = events.weightProduct[evt_selection.all_evt_cuts.all(*evt_selection.evt_cuts)]

                for h in hists.values():
                    h.fill(sel_objs, evt_weights)
                
                 # make cutflow
                if not lj_reco in cutflows:
                    cutflows[str(lj_reco)] = {}
                cutflows[str(lj_reco)][channel] = cutflow.Cutflow(evt_selection.all_evt_cuts, evt_selection.evt_cuts, events.weightProduct)
                
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
        # fixme: can define other LJ variables 
        
        # take lepton jets directly from ntuples
        if lj_reco == 0:
            ljs = objs["ntuple_ljs"]

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
            ### >>>>>>>>> Fixme! Dummy placeholders! Replace with an actual value at some point.
            ljs["dRSpread"]= ak.num(ljs.constituents[ljs.constituents["type"] == 2])
            ljs["pfIsolation05"]= ak.num(ljs.constituents[ljs.constituents["type"] == 2])
            ljs["pfIsolationPtNoPU05"]= ak.num(ljs.constituents[ljs.constituents["type"] == 2])
            #### <<<<<<<<<
            ljs["muon_n"] = ak.num(ljs.constituents[(ljs.constituents["type"] == 3)
                                                    | (ljs.constituents["type"] == 8)],axis=-1)
            ljs["electron_n"] = ak.num(ljs.constituents[ljs.constituents["type"] == 2],axis=-1)
            ljs["photon_n"] = ak.num(ljs.constituents[ljs.constituents["type"] == 4],axis=-1)
 
            # Todo: to apply cuts to match cuts applied in ntuples, use the normal selections framework 
            # and add cuts to cuts.py and selections.yaml
           
            
        else:
            raise NotImplementedError(f"{lj_reco} is not a recognized LJ reconstruction choice")
        # pt order and return
        return self.order(ljs)

    def build_analysis_channels(self, objs):
        #fixme: not actually used anymore
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
