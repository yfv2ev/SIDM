"""Module to define the Selection class"""

# columnar analysis
import awkward as ak
from coffea.analysis_tools import PackedSelection


class Selection:
    """Class to represent the collection of cuts that define a selection

    A selection consists of object-level and event-level cuts. Object-level cuts slim object
    collections, and event-level cuts reject whole events. 
    """

    def __init__(self, name, events, obj_cuts, evt_cuts):
        self.name = name
        self.obj_cuts = obj_cuts
        self.evt_cuts = evt_cuts
        self.events = events

        # evaluate all available object cuts if not previously evaluated
        # result is independent of given selection, so store as class variable
        if not hasattr(type(self), "all_obj_cuts"):
            type(self).all_obj_cuts = self.evaluate_all_obj_cuts()

        # get object mask for given selection
        self.obj_mask = self.get_obj_mask()

        # evaluate all available event cuts
        self.all_evt_cuts = self.evaluate_all_evt_cuts()

    def evaluate_all_obj_cuts(self):
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
        """Define all available event-level cuts"""
        # define objects to be used in cuts
        evts = self.events
        pvs = self.events.pv[self.obj_mask["pv"]]
        ljs = self.events.ljsource[self.obj_mask["ljsource"]]
        ljs = ljs[:, :2] # fixme: hacky temporary solution to only keep leading 2 LJs
        mu_ljs = ljs[(ljs["type"] == 3) | (ljs["type"] == 8)] # 3=pfmuon, 8=dsamuon
        egm_ljs = ljs[(ljs["type"] == 2) | (ljs["type"] == 4)] # 2=pfelectron, 4=pfphoton

        # evaluate cuts
        all_evt_cuts = PackedSelection()
        all_evt_cuts.add("PV filter", ak.num(pvs) >= 1)
        all_evt_cuts.add("Cosmic veto", evts.cosmicveto.result)
        all_evt_cuts.add(">=2 LJs", ak.num(ljs) >= 2)
        all_evt_cuts.add("4mu", ak.num(mu_ljs) == 2)
        all_evt_cuts.add("2mu2e", (ak.num(mu_ljs) == 1) & (ak.num(egm_ljs) == 1))

        return all_evt_cuts

    def get_obj_mask(self):
        obj_mask = {}
        for obj, cuts in self.obj_cuts.items():
            obj_mask[obj] = self.all_obj_cuts[obj][cuts[0]]
            for cut in cuts[1:]:
                obj_mask[obj] = obj_mask[obj] & self.all_obj_cuts[obj][cut]
        return obj_mask

    def apply_evt_cuts(self):
        raise NotImplementedError
