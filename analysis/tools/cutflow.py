"""Module to define the Cutflow and CutflowElement classes"""

# python
from tabulate import tabulate
# columnar analysis
from coffea import processor
from coffea.analysis_tools import PackedSelection
import awkward as ak


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

    def identity(self):
        """Create additive identity Cutflow to allow accumlator behavior"""
        return Cutflow(PackedSelection(), self.selection, self.zero_weights)

    def add(self, other):
        """Add two cutflows"""
        for i, _ in enumerate(self.flow):
            self.flow[i] = self.flow[i] + other.flow[i]

    def print_table(self, fraction=False):
        """Print simple cutflow table to stdout"""
        if fraction:
            data = []
            for e in self.flow:
                e.calculate_fractions()
                data.append([e.cut, 100*e.f_ind, 100*e.f_mar, 100*e.f_all])
            headers = [
                "cut name",
                "individual %",
                "marginal %",
                "cumulative %",
            ]
        else:
            data = [[e.cut, e.n_ind, e.n_all] for e in self.flow]
            headers = [
                "cut name",
                "individual cut N",
                "all cut N",
            ]
        print(tabulate(data, headers, floatfmt=".1f"))


class CutflowElement(processor.AccumulatorABC):
    """Class to represent individual rows of a cutflow table"""

    def __init__(self, cut, all_cuts, cutflow, weights, is_first_element=False):
        """Create each cutflow table row"""
        self.cut = cut
        self.cutflow = cutflow
        self.n_evts = ak.sum(weights)
        self.is_first_element = is_first_element

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

    def calculate_fractions(self):
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
                previous_element = self.cutflow.flow[self.cutflow.flow.index(self) - 1]
                self.f_mar = self.n_all / previous_element.n_all
            except ZeroDivisionError:
                self.f_mar = 0.0
