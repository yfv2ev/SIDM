"""Module to define the Selection class"""

# columnar analysis
from coffea.analysis_tools import PackedSelection
# local
import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '../.')) # This is definitely a bit hacky, but I don't mind
from definitions.cuts import obj_cut_defs, evt_cut_defs


class Selection:
    """Class to represent the collection of cuts that define a selection

    A selection consists of object-level and event-level cuts. Object-level cuts slim object
    collections, and event-level cuts reject whole events. Object-level cuts are stored as a
    dictionary of masks, and event-level cuts are stored as a PackedSelection.

    All available cuts are defined in sidm.definitions.cuts. The specific cuts that define each
    selection are accepted by Selection() as lists of strings.
    """

    def __init__(self, name, cuts):
        self.name = name
        self.obj_cuts = cuts["obj_cuts"] # dictionary of names of cuts to be applied
        self.evt_cuts = cuts["evt_cuts"] # list of names of cuts to be applied
        self.all_evt_cuts = PackedSelection() # will be filled later when cuts are evaluated
        self.obj_masks = None

    def evaluate_obj_cuts(self, objs, evaluated_obj_cuts):
        """Evaluate all relevant object-level cuts that have not already been evaluated"""
        for obj, cuts in self.obj_cuts.items():
            if obj not in evaluated_obj_cuts:
                evaluated_obj_cuts[obj] = {}
            for cut in cuts:
                if cut not in evaluated_obj_cuts[obj]:
                    evaluated_obj_cuts[obj][cut] = obj_cut_defs[obj][cut](objs)
        return evaluated_obj_cuts

    def make_obj_masks(self, evaluated_obj_cuts):
        """Create one mask per object for all cuts in obj_cuts"""
        obj_masks = {}
        for obj, cuts in self.obj_cuts.items():
            obj_masks[obj] = evaluated_obj_cuts[obj][cuts[0]]
            for cut in cuts[1:]:
                obj_masks[obj] = obj_masks[obj] & evaluated_obj_cuts[obj][cut]
        self.obj_masks = obj_masks

    def apply_obj_masks(self, objs):
        """Filter object collections based on object masks """
        sel_objs = {}
        for name, obj in objs.items():
            # filter objects if mask exists, return collection unfiltered if mask does not exist
            sel_objs[name] = obj[self.obj_masks[name]] if name in self.obj_masks else obj
        return sel_objs

    def apply_evt_cuts(self, objs):
        """Evaluate all event cuts and apply results to object collections"""
        # evaluate all selected cuts
        for cut in self.evt_cuts:
            self.all_evt_cuts.add(cut, evt_cut_defs[cut](objs))

        # apply event cuts to object collections
        sel_objs = {}
        for name, obj in objs.items():
            sel_objs[name] = obj[self.all_evt_cuts.all(*self.evt_cuts)]
        return sel_objs
