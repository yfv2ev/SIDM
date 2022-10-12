"""Module to define the base SIDM processor"""

# columnar analysis
from coffea import processor
from coffea.analysis_tools import PackedSelection
import hist
import awkward as ak
#local
from analysis.tools import cutflow
# always reload local modules to pick up changes during development
import importlib
importlib.reload(cutflow)


class SidmProcessor(processor.ProcessorABC):
    """Class to perform first processor tests.

    I expect the contents to evolve along the following lines:
        1. do a few basic tests
        2. implement a basic SIDM preselection
        3. make a few test histograms
        4. then convert this into the base SIDM processor from which all other SIDM
        processors will inherit (spinning the selection and histogram details off to other files
        in the process)
    """

    def __init__(self):
        """No need to do anything here -- everything happens in process"""

    def process(self, events):
        """Apply selections, make histograms and cutflow"""

        # handle metadata
        sample = events.metadata["sample"]

        # define PV selection (just a test; PV filter is already applied in FF nTuples)
        pvs = events.pv
        pvs = pvs[
            (pvs.ndof > 4)
            & (abs(pvs.z) < 24) # assume cm
            & (abs(pvs.rho) < 0.2) # assume cm (fixme: weinan disagrees: https://github.com/phylsix/Firefighter/blob/master/recoStuff/python/ffPrimaryVertexFilter_cfi.py)
        ]

        # define LJ selection
        ljs = events.ljsource
        ljs = ljs[ak.argsort(ljs.p4.pt, ascending=False)] # pt ordering
        ljs = ljs[
            (ljs.p4.pt > 30) # GeV, not applied in ntuples
            & (abs(ljs.p4.eta) < 2.4) # already applied in ntuples
            # fixme: add muon number and charge constraints if not already applied in ntuples
        ]

        # define event selection
        selection = PackedSelection()
        selection.add("PV filter", ak.num(pvs) >= 1)
        selection.add("Cosmic veto", events.cosmicveto.result)
        selection.add(">=2 LJs", ak.num(ljs) >= 2)

        # define hists
        lj_pt_axis = hist.axis.Regular(100, 0, 100, name="lj_pt", label="Lepton jet pT [GeV]")
        lj_type_axis = hist.axis.IntCategory([2, 3, 4, 8], name="lj_type")
        hists = {
            # pv
            "pv_n" : hist.Hist.new.Regular(100, 0, 100, name="pv_n").Int64(),
            "pv_ndof" : hist.Hist.new.Regular(20, 0, 20, name="pv_ndof").Int64(),
            "pv_z" : hist.Hist.new.Regular(100, -50, 50, name="pv_z").Double(),
            "pv_rho" : hist.Hist.new.Regular(100, -0.5, 0.5, name="pv_rho").Double(),
            # lj
            "lj_n" : hist.Hist.new.Regular(10, 0, 10, name="lj_n").Int64(),
            "lj_charge" : hist.Hist.new.Regular(10, -5, 5, name="lj_charge").Int64(),
            "lj_pt_type" : hist.Hist(lj_pt_axis, lj_type_axis),
            "lj_0_pt" : hist.Hist(lj_pt_axis),
            "lj_1_pt" : hist.Hist(lj_pt_axis),
            "lj_2_pt" : hist.Hist(lj_pt_axis),
            "lj_eta" : hist.Hist.new.Regular(100, -3, 3, name="lj_eta").Double(),
            "lj_phi" : hist.Hist.new.Regular(100, -3.14, 3.14, name="lj_phi").Double(),
        }

        # apply full selection
        cf = cutflow.Cutflow(selection)
        events = events[selection.all(*selection.names)]
        # update object collections to only include selected events
        # fixme: could just apply object cuts directly to events.object instead
        pvs = pvs[selection.all(*selection.names)]
        ljs = ljs[selection.all(*selection.names)]

        # fill hists
        # pv
        hists["pv_n"].fill(pv_n=ak.num(pvs))
        hists["pv_ndof"].fill(pv_ndof=ak.flatten(pvs.ndof))
        hists["pv_z"].fill(pv_z=ak.flatten(pvs.z))
        hists["pv_rho"].fill(pv_rho=ak.flatten(pvs.rho))
        # lj
        hists["lj_n"].fill(lj_n=ak.num(ljs))
        hists["lj_charge"].fill(lj_charge=ak.flatten(ljs.charge))
        hists["lj_pt_type"].fill(lj_pt=ak.flatten(ljs.p4.pt), lj_type=ak.flatten(ljs['type']))
        hists["lj_0_pt"].fill(lj_pt=ljs[ak.num(ljs) > 0, 0].p4.pt)
        hists["lj_1_pt"].fill(lj_pt=ljs[ak.num(ljs) > 1, 1].p4.pt)
        hists["lj_2_pt"].fill(lj_pt=ljs[ak.num(ljs) > 2, 2].p4.pt)
        hists["lj_eta"].fill(lj_eta=ak.flatten(ljs.p4.eta))
        hists["lj_phi"].fill(lj_phi=ak.flatten(ljs.p4.phi))

        out = {
            "cutflow" : cf,
            "hists" : hists,
        }
        return {sample : out}

    def postprocess(self, accumulator):
        raise NotImplementedError
