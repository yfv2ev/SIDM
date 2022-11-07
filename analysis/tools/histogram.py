"""Module to define the Histogram and Axis classes"""

# columnar analysis
import hist
import awkward as ak


class Histogram:
    """Class to represent histograms

    Histogram mostly exists so that hist.Hists and the appropriate filling arguments can be
    defined in one place.
    """

    def __init__(self, axes, storage="weight", weight_key=None):
        self.axes = axes
        self.storage = storage
        self.weight_key = weight_key
        self.hist = None

    def make_hist(self, channels=None):
        """Build associated hist.Hist

        Perform outside __init__ because channels aren't known until runtime.
        """

        # optionally add channels axis to hist
        if channels is not None:
            channel_axis = hist.axis.StrCategory(channels, name="channel")
            self.axes = [Axis(channel_axis, lambda objs: objs["ch"])] + self.axes

        axes = [a.axis for a in self.axes]
        self.hist = hist.Hist(*axes, storage=self.storage)

    def fill(self, objs, wgts):
        """Fill associated hist.Hist"""
        fill_args = {a.name : a.fill_func(objs) for a in self.axes}
        fill_args["weight"] : wgts[self.weight_key]
        self.hist.fill(**fill_args)


class Axis:
    """Class to represent histogram axes

    Axis just bundles together hist.axis objects and functions to fill them.
    """

    def __init__(self, axis, fill_func):
        self.axis = axis
        self.name = self.axis.name
        self.fill_func = fill_func
