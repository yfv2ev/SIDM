"""Module to define the Selection class"""

# columnar analysis
from coffea.analysis_tools import PackedSelection
# local
from sidm.definitions.cuts import evt_cut_defs


class Selection:
    """Class to represent the collection of cuts that define a selection

    A selection consists of event-level cuts which reject whole events. 
    Cuts are stored as a PackedSelection.

    All available cuts are defined in sidm.definitions.cuts. The specific cuts that define each
    selection are accepted by Selection() as lists of strings.
    """

    def __init__(self, cuts):
        self.evt_cuts = cuts # list of names of cuts to be applied
        self.all_evt_cuts = PackedSelection() # will be filled later when cuts are evaluated

    def apply_evt_cuts(self, objs,verbose=False):
        """Evaluate all event cuts and apply results to object collections"""
        # evaluate all selected cuts
        for cut in self.evt_cuts:
            if verbose:
                print("Adding ",cut)
            self.all_evt_cuts.add(cut, evt_cut_defs[cut](objs))

        # apply event cuts to object collections
        sel_objs = {}
        for name, obj in objs.items():
            sel_objs[name] = obj[self.all_evt_cuts.all(*self.evt_cuts)]
        return sel_objs
