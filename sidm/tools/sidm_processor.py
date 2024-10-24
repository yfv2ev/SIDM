"""Module to define the base SIDM processor"""

# python
import copy
import numpy as np
# columnar analysis
from coffea import processor
from coffea.nanoevents.methods import nanoaod
from coffea.nanoevents.methods import vector as cvec
import awkward as ak
import fastjet
import vector
#local
from sidm import BASE_DIR
from sidm.tools import selection, cutflow, utilities
from sidm.definitions.hists import hist_defs, counter_defs
from sidm.definitions.objects import preLj_objs, postLj_objs


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
        lj_reco_choices=["0.4"],
        selections_cfg=f"{BASE_DIR}/configs/selections.yaml",
        histograms_cfg=f"{BASE_DIR}/configs/hist_collections.yaml",
        unweighted_hist=False,
        verbose=False,
    ):
        self.channel_names = channel_names
        self.hist_collection_names = hist_collection_names
        self.lj_reco_choices = lj_reco_choices
        self.selections_cfg = selections_cfg
        self.histograms_cfg = histograms_cfg
        self.unweighted_hist = unweighted_hist
        self.obj_defs = preLj_objs
        self.verbose = verbose

    def process(self, events):
        """Apply selections, make histograms and cutflow"""

        # create object collections
        objs = {}
        for obj_name, obj_def in self.obj_defs.items():
            try:
                obj = obj_def(events)
            except AttributeError:
                print(f"Warning: {obj_name} not found in this sample. Skipping.")
                continue
            objs[obj_name] = obj

            # pt order
            objs[obj_name] = self.order(objs[obj_name])

            # use nanoevents.Muon behaviors for dsa muons
            if obj_name == "dsaMuons":
                forms = {f : objs[obj_name][f] for f in objs[obj_name].fields}
                objs[obj_name] = ak.zip(forms, with_name="Muon", behavior=nanoaod.behavior)

            # add lxy attribute to particles with children
            if hasattr(obj, "children"):
                objs[obj_name]["lxy"] = utilities.lxy(objs[obj_name])

            # add dimension to one-per-event objects to allow independent obj and evt cuts
            # skip objects with no fields
            if objs[obj_name].ndim == 1 and "x" in obj.fields:
                counts = ak.ones_like(objs[obj_name].x, dtype=np.int32)
                objs[obj_name] = ak.unflatten(objs[obj_name], counts)

        cutflows = {}
        counters = {}

        # define histograms
        hists = self.build_histograms()

        ### Make list of all object-level cuts; define object-level, post-lj-level, and event-level cuts per channel
        all_obj_cuts, ch_cuts = self.build_cuts()

        # evaluate all object-level cuts
        obj_selection = selection.JaggedSelection(all_obj_cuts, self.verbose)
        obj_selection.evaluate_obj_cuts(objs)

        # loop through lj reco choices and channels, treating each lj+channel pair as a unique Selection
        for channel in self.channel_names:

            # apply object selection
            channel_objs = obj_selection.make_and_apply_obj_masks(objs, ch_cuts[channel]["obj"])

            for lj_reco in self.lj_reco_choices:

                sel_objs = channel_objs

                # reconstruct lepton jets
                sel_objs["ljs"] = self.build_lepton_jets(channel_objs, float(lj_reco))

                # apply obj selection to ljs
                lj_selection = selection.JaggedSelection(ch_cuts[channel]["lj"], self.verbose)
                lj_selection.evaluate_obj_cuts(sel_objs)
                sel_objs = lj_selection.make_and_apply_obj_masks(sel_objs, ch_cuts[channel]["lj"])

                # add post-lj objects to sel_objs
                for obj in postLj_objs:
                    sel_objs[obj] = postLj_objs[obj](sel_objs)

                # apply post-lj obj selection
                postLj_selection = selection.JaggedSelection(ch_cuts[channel]["postLj_obj"], self.verbose)
                postLj_selection.evaluate_obj_cuts(sel_objs)
                sel_objs = postLj_selection.make_and_apply_obj_masks(sel_objs, ch_cuts[channel]["postLj_obj"])

                # build Selection objects and apply event selection
                evt_selection = selection.Selection(ch_cuts[channel]["evt"], self.verbose)
                sel_objs = evt_selection.apply_evt_cuts(sel_objs)

                # fill all hists
                sel_objs["ch"] = channel
                sel_objs["lj_reco"] = lj_reco

                # define event weights
                if self.unweighted_hist:
                    evt_weights =  ak.ones_like(self.obj_defs["weight"](events)[evt_selection.all_evt_cuts.all(*evt_selection.evt_cuts)])
                else:
                    evt_weights = self.obj_defs["weight"](events)[evt_selection.all_evt_cuts.all(*evt_selection.evt_cuts)]

                # fill histograms for this channel+lj_reco pair
                for h in hists.values():
                    h.fill(sel_objs, evt_weights)

                # make cutflow
                if lj_reco not in cutflows:
                    cutflows[str(lj_reco)] = {}
                cutflows[str(lj_reco)][channel] = cutflow.Cutflow(evt_selection.all_evt_cuts, evt_selection.evt_cuts, self.obj_defs["weight"](events))

                # Fill counters
                if lj_reco not in counters:
                    counters[lj_reco] = {}
                counters[lj_reco][channel] = {}

                for name, counter in counter_defs.items():
                    try:
                        counters[lj_reco][channel][name] = counter(sel_objs)
                    except (KeyError, AttributeError) as e:
                        print(f"Warning: cannot fill counter {name}. Skipping.")

        # lose lj_reco dimension to cutflows if only one reco was run
        if len(self.lj_reco_choices) == 1:
            cutflows = cutflows[self.lj_reco_choices[0]]

        out = {
            "cutflow": cutflows,
            "hists": {n: h.hist for n, h in hists.items()}, # output hist.Hists, not Histograms
            "counters": counters
        }

        return {events.metadata["dataset"]: out}

    def make_vector(self, objs, collection, type_id=None, mass=None, charge=None):
        shape = ak.ones_like(objs[collection].pt)
        return vector.zip(
            {
                "part_type": objs[collection]["type"] if type_id is None else type_id*shape,
                "charge": objs[collection].charge if charge is None else charge*shape,
                "pt": objs[collection].pt,
                "eta": objs[collection].eta,
                "phi": objs[collection].phi,
                "mass": objs[collection].mass if mass is None else mass*shape,
            }
        )

    def build_lepton_jets(self, objs, lj_reco):
        """Reconstruct lepton jets according to defintion given by lj_reco"""
        # fixme: can define other LJ variables

        # take lepton jets directly from ntuples
        if lj_reco == 0:
            ljs = objs["ntuple_ljs"]

        # reconstruct lepton jets from scratch
        else:
            if lj_reco < 0: # Use ljsource collection
                lj_inputs = self.make_vector(objs, "ljsources")

            else: #Use electron/muon/photon/dsamuon collections with a custom distance parameter
                muon_inputs = self.make_vector(objs, "muons", type_id=3)
                dsa_inputs = self.make_vector(objs, "dsaMuons", type_id=8, mass=0.106)
                ele_inputs = self.make_vector(objs, "electrons", type_id=2)
                photon_inputs = self.make_vector(objs, "photons", type_id=4, charge=0)
                lj_inputs = ak.concatenate([muon_inputs, dsa_inputs, ele_inputs, photon_inputs],
                                           axis=-1)

            distance_param = abs(lj_reco)
            jet_def = fastjet.JetDefinition(fastjet.antikt_algorithm, distance_param)

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
            const_vec = ak.zip(
                {"x": cluster.constituents().x,
                 "y": cluster.constituents().y,
                 "z": cluster.constituents().z,
                 "t": cluster.constituents().t,
                 "charge": cluster.constituents().charge,
                  "part_type":cluster.constituents().part_type},
                 with_name="LorentzVector",
                 behavior=cvec.behavior)

            ljs["constituents"] = const_vec
            ljs["pfMuons"] = ljs.constituents[ljs.constituents.part_type == 3]
            ljs["dsaMuons"] = ljs.constituents[ljs.constituents.part_type == 8]
            ljs["muons"] = ljs.constituents[(ljs.constituents.part_type == 3)
                                            | (ljs.constituents["part_type"] == 8)]
            ljs["electrons"] = ljs.constituents[ljs.constituents.part_type == 2]
            ljs["photons"] = ljs.constituents[ljs.constituents.part_type == 4]

            #Confusing to read, but to calculate dRSpread (the maximum dR betwen any pair of constituents in each lepton jet):
            #a) for each constituent, find the dR between it and all other constituents in the same LJ
            #b) flatten that into a list of dRs per LJ
            #c) and then take the maximum dR per LJ, leaving us with a single value per LJ
            ljs["dRSpread"]= ak.max( ak.flatten(const_vec.metric_table(const_vec, axis = 2) , axis = -1) ,  axis = -1)

            ljs["pfMu_n"] = ak.num(ljs.constituents[ljs.constituents.part_type == 3], axis=-1)
            ljs["dsaMu_n"] = ak.num(ljs.constituents[ljs.constituents.part_type == 8], axis=-1)
            ljs["muon_n"] = ak.num(ljs.constituents[(ljs.constituents["part_type"] == 3)
                                                    | (ljs.constituents["part_type"] == 8)],axis=-1)
            ljs["electron_n"] = ak.num(ljs.constituents[ljs.constituents["part_type"] == 2],axis=-1)
            ljs["photon_n"] = ak.num(ljs.constituents[ljs.constituents["part_type"] == 4],axis=-1)
            ljs["pfMu_n"] = ak.num(ljs.constituents[ljs.constituents.part_type == 3], axis=-1)
            ljs["dsaMu_n"] = ak.num(ljs.constituents[ljs.constituents.part_type == 8], axis=-1)

            # pt order the new LJs
            ljs = self.order(ljs)

        # return the new LJ collection
        return ljs

    def build_cuts(self):
        """ Make list of pre-lj object, lj, post-lj obj, and event cuts per channel"""

        selection_menu = utilities.load_yaml(self.selections_cfg)

        all_obj_cuts = {}
        ch_cuts = {}

        for channel in self.channel_names:
            ch_cuts[channel] = {}
            ch_cuts[channel]["obj"] = {}
            ch_cuts[channel]["lj"] = {}
            ch_cuts[channel]["postLj_obj"] = {}
            ch_cuts[channel]["evt"] = {}

            #Merge all object level cuts into a single list to be evaluated once
            cuts = selection_menu[channel]
            for obj, obj_cuts in cuts["obj_cuts"].items():
                if obj not in all_obj_cuts:
                    all_obj_cuts[obj] = []
                all_obj_cuts[obj] = utilities.add_unique_and_flatten(all_obj_cuts[obj],obj_cuts)

                if obj not in ch_cuts[channel]["obj"]:
                    ch_cuts[channel]["obj"][obj] = []
                ch_cuts[channel]["obj"][obj] = utilities.flatten(obj_cuts)

            if "postLj_obj_cuts" in cuts:
                for obj, obj_cuts in cuts["postLj_obj_cuts"].items():
                    if obj == "ljs":
                        ch_cuts[channel]["lj"][obj] = utilities.flatten(obj_cuts)
                    else:
                        ch_cuts[channel]["postLj_obj"][obj] = utilities.flatten(obj_cuts)

            if "evt_cuts" in cuts:
                ch_cuts[channel]["evt"] = utilities.flatten(cuts["evt_cuts"])

        return all_obj_cuts, ch_cuts

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
                hists[hist_name].make_hist(hist_name, self.channel_names, lj_reco_names)
        return hists

    def order(self, obj):
        """Explicitly order objects"""
        # pt order objects with a pt attribute
        if hasattr(obj, "pt"):
            obj = obj[ak.argsort(obj.pt, ascending=False)]
        # fixme: would be good to explicitly order other objects as well
        return obj

    def postprocess(self, accumulator):
        pass
