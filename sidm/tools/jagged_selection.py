"""Module to define the JaggedSelection class"""

# local
from sidm.definitions.cuts import obj_cut_defs


class JaggedSelection:
    """Class to represent the collection of cuts that define a selection

    A selection consists of object-level cuts (including lepton-jet-level cuts) and event-level cuts. 
    Object-level cuts slim object collections, and event-level cuts reject whole events. Object-level cuts are stored as a
    dictionary of masks, and event-level cuts are stored as a PackedSelection.

    All available cuts are defined in sidm.definitions.cuts. The specific cuts that define each
    selection are accepted by Selection() as lists of strings.
    """

    def __init__(self, cuts):
        self.obj_cuts = cuts # dict of names of cuts to be applied
        self.evaluated_obj_cuts = {}        

    def evaluate_obj_cuts(self, objs,verbose=False):
        """Evaluate all relevant object-level cuts that have not already been evaluated"""
        for obj, cuts in self.obj_cuts.items():            
            if obj not in self.evaluated_obj_cuts:
                self.evaluated_obj_cuts[obj] = {}
            for cut in cuts:
                if cut not in self.evaluated_obj_cuts[obj]:
                    if verbose:
                        print("Evaluating ", obj," ",cut)
                    self.evaluated_obj_cuts[obj][cut] = obj_cut_defs[obj][cut](objs)

    def make_obj_masks(self, channel_cut_list,verbose=False):
        """Create one mask per object, using the subset of cuts specified in channel_cut_list"""
        obj_masks = {}
        for obj, cuts in channel_cut_list.items():
            for cut in cuts:
                if cut not in self.evaluated_obj_cuts[obj]:
                    print("Uh oh, haven't evaluated this cut yet! Make sure it was included in the list of cuts you used to initialize this JaggedSelection.  ", obj, ": ",cut)
                else:      
                    if verbose:
                        print("Adding the following cut on ",obj, "to the mask: ",cut)
                    if obj not in obj_masks:
                        obj_masks[obj] = self.evaluated_obj_cuts[obj][cut]
                    else:
                        obj_masks[obj] = obj_masks[obj] & self.evaluated_obj_cuts[obj][cut]
        return obj_masks

    def apply_obj_masks(self, objs, obj_masks,verbose=False):
        """Filter object collections based on object masks """
        
        sel_objs = {}
        for name, obj in objs.items():
            # filter objects if mask exists, return collection unfiltered if mask does not exist
            sel_objs[name] = obj[obj_masks[name]] if name in obj_masks else obj
            
            if verbose:
                if name in obj_masks: 
                    print("Applying mask to collection: ", name)
                else:
                    print("No mask available for collection; returning unfiltered: ",name)
            
        return sel_objs
    
    
    def make_and_apply_obj_masks(self, objs, channel_cut_list, verbose=False):
        return self.apply_obj_masks(objs, self.make_obj_masks(channel_cut_list,verbose), verbose)

