"""Module to define the Selection class"""

# columnar analysis
import awkward as ak
from coffea.analysis_tools import PackedSelection


class Selection:
    """Class to represent the collection of cuts that define a selection

    A selection consists of object-level and event-level cuts. Object-level cuts slim object
    collections, and event-level cuts reject whole events. Object-level cuts are stored as a
    dictionary of masks, and event-level cuts are stored as a PackedSelection.

    All available cuts are defined in evaluate_all_obj_cuts() and evaluate_all_evt_cuts(). The
    specific cuts that define each selection are accepted by Selection() as lists of strings.
    """

    def __init__(self, name, cuts, events):
        self.name = name
        self.obj_cuts = cuts["obj_cuts"] # dictionary of names of cuts to be applied
        self.evt_cuts = cuts["evt_cuts"] # list of names of cuts to be applied
        self.events = events

        # evaluate all available object cuts if not previously evaluated
        # result is independent of given selection, so store as class variable
        if not hasattr(type(self), "all_obj_cuts"):
            type(self).all_obj_cuts = self.evaluate_all_obj_cuts()

        # get object mask for given selection
        self.obj_masks = self.make_obj_masks()

        # evaluate all available event cuts
        # result depends on chosen object-level cuts, so store as instance variable
        self.all_evt_cuts = self.evaluate_all_evt_cuts()

    def evaluate_all_obj_cuts(self):
        """Evaluate all available object-level cuts"""
        # define objects to be used in cuts
        pvs = self.events.pv
        ljs = self.events.ljsource

        # evaluate cuts
        all_obj_cuts = {
            "pv" : {
                "ndof > 4" : pvs.ndof > 4,
                "|z| < 24 cm" : abs(pvs.z) < 24,
                "|rho| < 0.2 mm" : abs(pvs.rho) < 0.2, # fixme: double check mm
            },
            "ljsource" : {
                "pT > 30 GeV" : ljs.p4.pt > 30,
                "|eta| < 2.4" : abs(ljs.p4.eta) < 2.4,
                # fixme: figure out way to implement ljs = ljs[:, :2] and similar
            }
        }
        return all_obj_cuts

    def evaluate_all_evt_cuts(self):
        """Evaluate all available event-level cuts"""
        # define objects to be used in cuts
        evts = self.events
        pvs = self.events.pv[self.obj_masks["pv"]]
        ljs = self.events.ljsource[self.obj_masks["ljsource"]]
        ljs = ljs[:, :2] # fixme: hacky temporary solution to only keep leading 2 LJs
        mu_ljs = ljs[(ljs["type"] == 3) | (ljs["type"] == 8)] # 3=pfmuon, 8=dsamuon
        egm_ljs = ljs[(ljs["type"] == 2) | (ljs["type"] == 4)] # 2=pfelectron, 4=pfphoton

        # evaluate cuts and store results as PackedSelection
        all_evt_cuts = PackedSelection()
        all_evt_cuts.add("PV filter", ak.num(pvs) >= 1)
        all_evt_cuts.add("Cosmic veto", evts.cosmicveto.result)
        all_evt_cuts.add(">=2 LJs", ak.num(ljs) >= 2)
        all_evt_cuts.add("4mu", ak.num(mu_ljs) == 2)
        all_evt_cuts.add("2mu2e", (ak.num(mu_ljs) == 1) & (ak.num(egm_ljs) == 1))

        return all_evt_cuts

    def make_obj_masks(self):
        """Create one mask per object for all cuts in obj_cuts"""
        obj_masks = {}
        for obj, cuts in self.obj_cuts.items():
            obj_masks[obj] = self.all_obj_cuts[obj][cuts[0]]
            for cut in cuts[1:]:
                obj_masks[obj] = obj_masks[obj] & self.all_obj_cuts[obj][cut]
        return obj_masks
