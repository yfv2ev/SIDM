from coffea import processor

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

        out = {"sumw" : len(events)}
        return {sample : out}

    def postprocess(self, accumulator):
        raise NotImplementedError