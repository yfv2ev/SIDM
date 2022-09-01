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
        for cut in selection.names:
            self.flow.append(self.CutflowElement(cut, selection))

    def print_table(self):
        """Print simple cutflow table to stdout"""
        data = [[e.cut, e.n_ind, e.n_all] for e in self.flow]
        headers = ["cut name", "passing individual cut", "passing all cuts"]
        print(tabulate(data, headers))


    class CutflowElement:
        """Class to represent individual rows of a cutflow table"""

        def __init__(self, cut, selection):
            """Create each cutflow table row"""
            self.cut = cut
            if cut == "No selection":
                n_evts = len(selection.all(selection.names[0]))
                self.n_ind = n_evts
                self.n_all = n_evts
            else:
                all_cuts = selection.names[:selection.names.index(cut) + 1]
                self.n_ind = list(selection.all(cut)).count(True)
                self.n_all = list(selection.all(*all_cuts)).count(True)
 