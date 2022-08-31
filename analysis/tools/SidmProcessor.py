#python
import importlib
# columnar analysis
from coffea import processor
from coffea.analysis_tools import PackedSelection
import hist
import awkward as ak
#local
from analysis.tools import Cutflow
importlib.reload(Cutflow)

class SidmProcessor(processor.ProcessorABC):
    """
    Class to perform first processor tests. I expect the contents to evolve along the following
    lines:
        1. do a few basic tests
        2. implement a basic SIDM preselection
        3. make a few test histograms 
        4. then convert this into the base SIDM processor from which all other SIDM
        processors will inherit (spinning the selection and histogram details off to other files
        in the process)
    """
    
    def __init__(self):
        pass

    def process(self, events):
        sample = events.metadata["sample"]
        
        # define PV selection (just a test; PV filter is already applied in FF nTuples)
        pvs = events.pv
        pvs = pvs[
            (pvs.ndof > 4)
            & (abs(pvs.z) < 24) # assume cm
            & (abs(pvs.rho) < 0.2) # assume cm (fixme: weinan disagrees: https://github.com/phylsix/Firefighter/blob/master/recoStuff/python/ffPrimaryVertexFilter_cfi.py)
        ]

        # define event selection
        selection = PackedSelection()
        selection.add("PV filter", ak.num(pvs) >= 1)
        selection.add("Cosmic veto", events.cosmicveto.result)

        # define hists
        hists = {
            "pv_z" : hist.Hist.new.Regular(100, -50, 50, name="pv_z").Double(),
            "pv_rho" : hist.Hist.new.Regular(100, -0.5, 0.5, name="pv_rho").Double(),
        }

        # apply full selection
        cutflow = Cutflow.Cutflow(events, selection)
        events = events[selection.all(*selection.names)]

        # fill hists
        hists["pv_z"].fill(pv_z=ak.flatten(pvs.z))
        hists["pv_rho"].fill(pv_rho=ak.flatten(pvs.rho))

        out = {
            "cutflow" : cutflow,
            "hists" : hists,
        }
        return {sample : out}

    def postprocess(self, accumulator):
        raise NotImplementedError