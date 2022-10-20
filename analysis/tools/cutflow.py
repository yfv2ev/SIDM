"""Module to define the Cutflow class"""

from tabulate import tabulate
import awkward as ak


class Cutflow:
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
        self.all_cuts = all_cuts
        self.selection = selection
        self.weights = weights
        self.flow = [self.CutflowElement("No selection", self)]
        for cut in selection:
            self.flow.append(self.CutflowElement(cut, self))

    def print_table(self, fraction=False):
        """Print simple cutflow table to stdout"""
        if fraction:
            data = [[e.cut, 100*e.f_ind, 100*e.f_mar, 100*e.f_all] for e in self.flow]
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


    class CutflowElement:
        """Class to represent individual rows of a cutflow table"""

        def __init__(self, cut, cutflow):
            """Create each cutflow table row"""
            self.cut = cut
            n_evts = ak.sum(cutflow.weights)

            if not hasattr(cutflow, "flow"):
                self.n_ind = n_evts
                self.n_all = n_evts
                self.f_ind = 1.0
                self.f_mar = 1.0
                self.f_all = 1.0
            else:
                cumulative_cuts = cutflow.selection[:cutflow.selection.index(cut) + 1]
                self.n_ind = ak.sum(cutflow.weights[cutflow.all_cuts.all(cut)])
                self.n_all = ak.sum(cutflow.weights[cutflow.all_cuts.all(*cumulative_cuts)])
                self.f_ind = self.n_ind / n_evts
                self.f_all = self.n_all / n_evts
                try:
                    self.f_mar = self.n_all / cutflow.flow[-1].n_all
                except ZeroDivisionError:
                    self.f_mar = 0.0
