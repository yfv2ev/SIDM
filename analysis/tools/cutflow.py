"""Module to define the Cutflow class"""

from tabulate import tabulate


class Cutflow:
    """Class to represent the number of events that pass each cut in a selection

    Cutflow currently stores the number of events that pass each cut individually and the
    number of events that pass the logical AND of the current and all preceding cuts.
    """

    def __init__(self, selection):
        """Make Cutflow, starting with 'No selection' row"""
        self.flow = [self.CutflowElement("No selection", selection)]
        previous_element = self.flow[0]
        for cut in selection.names:
            self.flow.append(self.CutflowElement(cut, selection, previous_element))
            previous_element = self.flow[-1]

    def print_table(self, fraction=False):
        """Print simple cutflow table to stdout"""
        if fraction:
            data = [[e.cut, 100*e.f_ind, 100*e.f_all, 100*e.f_mar] for e in self.flow]
            headers = [
                "cut name",
                "individual cut %",
                "all cuts %",
                "marginal %",
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

        def __init__(self, cut, selection, previous_element=None):
            """Create each cutflow table row"""
            self.cut = cut
            n_evts = len(selection.all(selection.names[0]))

            if cut == "No selection":
                self.n_ind = n_evts
                self.n_all = n_evts
                self.f_ind = 1.0
                self.f_all = 1.0
                self.f_mar = 1.0
            else:
                all_cuts = selection.names[:selection.names.index(cut) + 1]
                self.n_ind = list(selection.all(cut)).count(True)
                self.n_all = list(selection.all(*all_cuts)).count(True)
                self.f_ind = self.n_ind / float(n_evts)
                self.f_all = self.n_all / float(n_evts)
                self.f_mar = self.n_all / float(previous_element.n_all)
 