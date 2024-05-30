"""Module to define the Cutflow and CutflowElement classes"""

# python
from tabulate import tabulate
# columnar analysis
from coffea import processor
from coffea.analysis_tools import PackedSelection
import awkward as ak
import numpy as np

class Cutflow(processor.AccumulatorABC):
    """Class to represent the number of events that pass each cut in a selection

    Cutflow can print tables of the following values:
    - n_ind: number of events that pass each cut individually
    - n_all: number of events that pass the logical AND of the current and all preceding cuts
    - f_ind: fraction of events that pass each cut individually
    - f_mar: fraction of events passing all preceding cuts that pass the current cut
    - f_all: fraction of events that pass the logical AND of the current and all preceding cuts
    """

    def __init__(self, all_cuts, selection, weights):
        """Make Cutflow, starting with 'No selection' row"""
        self.selection = selection # list of cut names to apply
        # make behavior-free array with weights set to zero for making additive identity Cutflows
        self.zero_weights = ak.without_parameters(ak.zeros_like(weights), behavior={})

        # make all cutflow rows
        self.flow = [CutflowElement("No selection", all_cuts, self, weights, is_first_element=True)]
        for cut in selection:
            self.flow.append(CutflowElement(cut, all_cuts, self, weights))

        # make all unweighted cutflow rows
        one_weights = ak.without_parameters(ak.ones_like(weights), behavior={})
        self.unweighted_flow = [CutflowElement("No selection", all_cuts, self, one_weights,
                                               is_first_element=True)]
        for cut in selection:
            self.unweighted_flow.append(CutflowElement(cut, all_cuts, self, one_weights))

    def identity(self):
        """Create additive identity Cutflow to allow accumlator behavior"""
        all_cuts = PackedSelection()
        for cut in self.selection:
            all_cuts.add(cut, ak.values_astype(self.zero_weights, bool))
        return Cutflow(all_cuts, self.selection, self.zero_weights)

    def add(self, other):
        """Add two cutflows"""
        for i, _ in enumerate(self.flow):
            self.flow[i] = self.flow[i] + other.flow[i]
            self.unweighted_flow[i] = self.unweighted_flow[i] + other.unweighted_flow[i]

            
    def efficiency(self):
        """Outputs the fraction of events passing the cutflow as a fraction of 1"""
        return float(list(enumerate(self.flow))[-1][1].n_all / list(enumerate(self.flow))[-1][1].n_evts)
    
    def cut_breakdown(self, fraction=False, unweighted=False, giveCuts=False):
        """Outputs a list of the number of events passing each cut. Effectively isolates the cumulative column of the cut table"""
        """The giveCuts argument decides whether the function returns the column of cut names, useful for plotting / making a table"""
        flow = self.unweighted_flow if unweighted else self.flow
        data = []
        if giveCuts:
            data = [e.cut for e in flow]
        else:
            for i in range(len(list(enumerate(flow)))):
                data.append(list(enumerate(flow))[i][1].n_all)
            if fraction == True:
                data = [100.0 * x / list(enumerate(flow))[-1][1].n_evts for x in data]
        return data
            
    def print_table(self, fraction=False, unweighted=False):
        """Print simple cutflow table to stdout"""
        flow = self.unweighted_flow if unweighted else self.flow
        if fraction:
            data = []
            for i, e in enumerate(flow):
                previous_element = flow[i - 1] if i > 0 else None
                e.calculate_fractions(previous_element)
                data.append([e.cut, 100*e.f_ind, 100*e.f_mar, 100*e.f_all])
            headers = [
                "cut name",
                "individual %",
                "marginal %",
                "cumulative %",
            ]
        else:
            data = [[e.cut, e.n_ind, e.n_all] for e in flow]
            headers = [
                "cut name",
                "individual cut N",
                "all cut N",
            ]
        print(tabulate(data, headers, floatfmt=".1f"))
        

 
    def print_multi_table(self, cutflows, headers, fraction=False, unweighted=False, title=""):
        """Prints a table with multiple cutflows listed, one in each column. Total number of cuts on each sample are listed. 
        It would be better to make this its own function independent of the cutflow class. That would allow a complete list of cutflows to be passed
        rather than just removing one and calling the function from it."""
        data = np.array([self.cut_breakdown(fraction, unweighted, giveCuts=True), self.cut_breakdown(fraction, unweighted)])
        for cutflow in cutflows:
            data = np.append(data, [cutflow.cut_breakdown(fraction, unweighted)], axis=0)
        data = data.transpose()
        headerline = ["cut name",]
        for header in headers:
            if fraction == False:
                headerline.append("Total cuts: \n" + header)
            else:
                headerline.append("% cuts: \n" + header)
        if title != "":
            print(title)
            for header in headerline:
                for i in range(len(header)):
                    print("-", end='')
                print("----", end='')
        print('\n' + tabulate(data, headerline, floatfmt=".2f") + '\n')
        
    def n_input_evts(self, unweighted=False):
        """Return number of events in sample before applying any cuts"""
        flow = self.unweighted_flow if unweighted else self.flow
        return flow[0].n_evts

class CutflowElement(processor.AccumulatorABC):
    """Class to represent individual rows of a cutflow table"""

    def __init__(self, cut, all_cuts, cutflow, weights, is_first_element=False):
        """Create each cutflow table row"""
        self.cut = cut
        self.cutflow = cutflow
        self.n_evts = ak.sum(weights)
        self.is_first_element = is_first_element
        self.f_ind = None
        self.f_all = None
        self.f_mar = None

        if is_first_element or self.n_evts == 0:
            self.n_ind = self.n_evts
            self.n_all = self.n_evts
        else:
            cumulative_cuts = self.cutflow.selection[:self.cutflow.selection.index(cut) + 1]
            self.n_ind = ak.sum(weights[all_cuts.all(cut)])
            self.n_all = ak.sum(weights[all_cuts.all(*cumulative_cuts)])

    def identity(self):
        """Create additive identity CutflowElement"""
        return CutflowElement(self.cut, PackedSelection(), self.cutflow, self.cutflow.zero_weights,
                              self.is_first_element)

    def add(self, other):
        """Add two CutflowElements"""
        self.n_evts = self.n_evts + other.n_evts
        self.n_ind = self.n_ind + other.n_ind
        self.n_all = self.n_all + other.n_all

    def calculate_fractions(self, previous_element):
        """Calculate individual, cumulative, and marginal fractional cutflow values"""
        # only calculate if fractions have not already been calculated
        if self.is_first_element:
            self.f_ind = 1.0
            self.f_all = 1.0
            self.f_mar = 1.0
        else:
            self.f_ind = self.n_ind / self.n_evts
            self.f_all = self.n_all / self.n_evts
            try:
                # note that this produces runtime warnings when using coffea.processor.Runner
                self.f_mar = self.n_all / previous_element.n_all
            except ZeroDivisionError:
                self.f_mar = 0.0
