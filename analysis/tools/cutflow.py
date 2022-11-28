"""Module to define the Cutflow and CutflowElement classes"""

# python
from tabulate import tabulate
# columnar analysis
from coffea import processor
import awkward as ak
from coffea.analysis_tools import PackedSelection


class Cutflow(processor.AccumulatorABC):
    """Class to represent the number of events that pass each cut in a selection

    Cutflow currently stores and can print tables of the following values:
    - n_ind: number of events that pass each cut individually
    - n_all: number of events that pass the logical AND of the current and all preceding cuts
    - f_ind: fraction of events that pass each cut individually
    - f_mar: fraction of events passing all preceding cuts that pass the current cut
    - f_all: fraction of events that pass the logical AND of the current and all preceding cuts
    """

    def __init__(self, all_cuts, selection, weights):
        """Make Cutflow, starting with 'No selection' row"""
        self.all_cuts = all_cuts # PackedSelection of all relevant cuts
        self.selection = selection # list of cut names to apply
        self.weights = weights # array of event weights
        self.flow = [CutflowElement("No selection", self, first_element=True)]
        for cut in selection:
            self.flow.append(CutflowElement(cut, self))

    def identity(self):
        """Create additive identity Cutflow to allow accumlator behavior"""
        identity_all_cuts = PackedSelection()
        for cut in self.selection:
            identity_all_cuts.add(cut, ak.full_like(self.weights, False, dtype=bool))
        identity_selection = self.selection
        identity_weights = ak.zeros_like(self.weights)
        return Cutflow(identity_all_cuts, identity_selection, identity_weights)

    def add(self, other):
        """Add two cutflows"""
        for i in range(len(self.flow)):
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

    def __init__(self, cut, cutflow, first_element=False):
        """Create each cutflow table row"""
        self.cut = cut
        self.cutflow = cutflow
        self.n_evts = ak.sum(cutflow.weights)
        self.first_element = first_element

        if first_element:
            self.n_ind = self.n_evts
            self.n_all = self.n_evts
        else:
            cumulative_cuts = self.cutflow.selection[:self.cutflow.selection.index(cut) + 1]
            self.n_ind = ak.sum(self.cutflow.weights[self.cutflow.all_cuts.all(cut)])
            self.n_all = ak.sum(self.cutflow.weights[self.cutflow.all_cuts.all(*cumulative_cuts)])

    def identity(self):
        """Create additive identity CutflowElement"""
        identity = CutflowElement(self.cut, self.cutflow, self.first_element)
        identity.n_events = 0.0
        identity.n_ind = 0.0
        identity.n_all = 0.0
        return identity

    def add(self, other):
        """Add two CutflowElements"""
        self.n_evts = self.n_evts + other.n_evts
        self.n_ind = self.n_ind + other.n_ind
        self.n_all = self.n_all + other.n_all

    def calculate_fractions(self):
        """Calculate individual, cumulative, and marginal fractional cutflow values"""
        # only calculate if fractions have not already been calculated
        if self.first_element:
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
