"""Module to define the Cutflow class"""

from tabulate import tabulate


class Cutflow:
    """Class to represent the number of events that pass each cut in a selection

    Cutflow currently stores the number of events that pass each cut individually and the
    number of events that pass the logical AND of the current and all preceding cuts.
    """

    def __init__(self, all_cuts, selection):
        """Make Cutflow, starting with 'No selection' row"""
        self.flow = [self.CutflowElement("No selection", all_cuts, selection)]
        previous_element = self.flow[0]
        for cut in selection:
            self.flow.append(self.CutflowElement(cut, all_cuts, selection, previous_element))
            previous_element = self.flow[-1]

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

        def __init__(self, cut, all_cuts, selection, previous_element=None):
            """Create each cutflow table row"""
            self.cut = cut
            n_evts = len(all_cuts.all(selection[0]))

            if cut == "No selection":
                self.n_ind = n_evts
                self.n_all = n_evts
                self.f_ind = 1.0
                self.f_mar = 1.0
                self.f_all = 1.0
            else:
                cumulative_cuts = selection[:selection.index(cut) + 1]
                self.n_ind = list(all_cuts.all(cut)).count(True)
                self.n_all = list(all_cuts.all(*cumulative_cuts)).count(True)
                self.f_ind = self.n_ind / float(n_evts)
                self.f_all = self.n_all / float(n_evts)
                try:
                    self.f_mar = self.n_all / float(previous_element.n_all)
                except ZeroDivisionError:
                    self.f_mar = 0.0
