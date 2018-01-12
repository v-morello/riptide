##### Standard imports #####
import json

##### Non-standard imports #####
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

##### Local imports #####
from . import Candidate



class CandidatePlotLayout(object):
    """ Base class for any candidate plot layout. All the properties MUST be
    implemented in any derived class. """
    @property
    def default_figsize(self):
        raise NotImplementedError

    @property
    def default_dpi(self):
        raise NotImplementedError

    @property
    def axSubints(self):
        raise NotImplementedError

    @property
    def axProfile(self):
        raise NotImplementedError


class CandidatePlotDefaultLayout(CandidatePlotLayout):
    """ """
    default_figsize = (16, 7)
    default_dpi = 80

    tableLayout = {
        "topMargin": 0.05,
        "leftMargin": 0.05,
        "rightMargin": 0.05,
        "lineHeight": 0.10,
        "columnWidths" : [
            ("name", 1),
            ("value", 1),
            ("unit", 1),
            ]
        }

    def __init__(self, figsize=None, dpi=None):
        self.grid = gs.GridSpec(6, 8)
        self.figsize = self.default_figsize if figsize is None else figsize
        self.dpi = self.default_dpi if dpi is None else dpi

    @property
    def axSubints(self):
        return plt.subplot(self.grid[:4, 2:5])

    @property
    def axProfile(self):
        return plt.subplot(self.grid[4:, 2:5])

    @property
    def axDMCurve(self):
        return plt.subplot(self.grid[:2, 5:])

    @property
    def axPeriodCurve(self):
        return plt.subplot(self.grid[2:4, 5:])

    @property
    def axWidthCurve(self):
        return plt.subplot(self.grid[4:, 5:])

    @property
    def axTable(self):
        return plt.subplot(self.grid[:, :2])




def plot_text(text, ax, x, y, horizontalalignment='right'):
    ax.text(x, y, text,
        horizontalalignment=horizontalalignment, verticalalignment='center',
        family='monospace', transform=ax.transAxes, fontsize=12)

class Entry(object):
    """ """
    def __init__(self, name, value, unit=""):
        self.name = name
        self.value = value
        self.unit = unit

    def __str__(self):
        dic = {"name": self.name, "value": self.value, "unit": self.unit}
        return "Entry " + str(dic)

    def __repr__(self):
        return str(self)

class Table(object):
    """ """
    columnOffsets = {
        "name"  : 0.20,
        "value" : 0.85,
        "unit"  : 0.88
        }
    topMargin = 0.01
    lineHeight = 0.030

    def __init__(self):
        self.entries = []

    def add_entry(self, entry):
        self.entries.append(entry)

    def skip_line(self):
        self.entries.append(None)

    def plot(self, ax):
        # Hide all axes
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.patch.set_visible(False)
        for spine in ax.spines.values():
            spine.set_visible(False)

        # y = vertical offset from bottom
        y = 1.0 - self.topMargin
        for entry in self.entries:
            if entry is not None:
                plot_text(entry.name, ax, self.columnOffsets['name'], y)
                plot_text(entry.value, ax, self.columnOffsets['value'], y)
                plot_text(entry.unit, ax, self.columnOffsets['unit'], y, 'left')
            y -= self.lineHeight



class CandidatePlot(object):
    """
    """
    _layouts = {
        'default': CandidatePlotDefaultLayout,
        }

    def __init__(self, cand, layout='default', figsize=None, dpi=None):
        self.cand = cand
        self.layout = None
        self.figure = None
        self.shift_bins = 0

        self._set_layout(layout, figsize=figsize, dpi=dpi)
        self._compute_shift_bins()
        self._plot()

    def _set_layout(self, layout, figsize=None, dpi=None):
        if layout in self._layouts:
            self.layout = self._layouts[layout](figsize=figsize, dpi=dpi)
        else:
            msg = "Parameter \"layout\" must be one of: {!s}".format(self._layouts.keys())
            raise ValueError(msg)

        self.figure = plt.figure(
            figsize=self.layout.figsize,
            dpi=self.layout.dpi)

    def _compute_shift_bins(self):
        prof = self.cand.subints.normalised_profile
        self.shift_bins = prof.size//2 - prof.argmax()

    def _plotSubints(self):
        ax = self.layout.axSubints
        tobs = self.cand.metadata['tobs']
        extent = [0, 1, tobs, 0]
        ax.imshow(
            np.roll(self.cand.subints.data, self.shift_bins, axis=1),
            aspect='auto',
            cmap=plt.cm.Greys,
            extent=extent)
        plt.xticks([])
        plt.ylabel("Time (s)")
        plt.title("Sub-Integrations")

    def _plotProfile(self):
        ax = self.layout.axProfile
        prof = np.roll(self.cand.subints.normalised_profile, self.shift_bins)
        ax.bar(
            np.arange(0, prof.size),
            prof,
            width=1,
            color='#303030')
        plt.xlim(-0.5, prof.size-0.5)
        plt.title("Integrated Profile")

    def _plotDMCurve(self):
        ax = self.layout.axDMCurve
        cv = self.cand.dm_curve
        ax.plot(cv.trials, cv.snr, marker='o', markersize=3)
        plt.xlim(cv.trials[0], cv.trials[-1])
        plt.grid(linestyle=':')
        plt.title("DM Curve")

    def _plotPeriodCurve(self):
        ax = self.layout.axPeriodCurve
        cv = self.cand.period_curve
        ax.plot(cv.trials, cv.snr)
        plt.xlim(cv.trials[0], cv.trials[-1])
        plt.grid(linestyle=':')
        plt.title("Period Curve")

    def _plotWidthCurve(self):
        pass

    def _plotTable(self):
        table = Table()
        table.add_entry(Entry("Source", self.cand.metadata['source_name']))

        coord = self.cand.metadata['skycoord']
        ra = coord.ra.to_string(unit='hour', precision=2, pad=2)
        dec = coord.dec.to_string(unit='degree', precision=2, pad=2)
        glon_deg = coord.galactic.l.deg
        glat_deg = coord.galactic.b.deg
        tobs = self.cand.metadata['tobs']

        table.add_entry(Entry("RA", ra))
        table.add_entry(Entry("Dec", dec))
        table.add_entry(Entry("glon", "{:.3f}".format(glon_deg), "deg"))
        table.add_entry(Entry("glat", "{:.3f}".format(glat_deg), "deg"))
        table.add_entry(Entry("Tobs", "{:.2f}".format(tobs), "s"))
        table.skip_line()

        period = self.cand.metadata['best_period']
        table.add_entry(Entry("Period", "{:.6f}".format(1000.0 * period), "ms"))

        dm = self.cand.metadata['best_dm']
        table.add_entry(Entry("DM", "{:.3f}".format(dm)))

        snr = self.cand.metadata['best_snr']
        table.add_entry(Entry("S/N", "{:.1f}".format(snr)))
        table.plot(self.layout.axTable)

    def _plot(self):
        self._plotProfile()
        self._plotSubints()
        self._plotDMCurve()
        self._plotPeriodCurve()
        self._plotTable()
        self.figure.tight_layout()

    def display(self):
        self.figure.show()

    def saveimg(self, fname):
        self.figure.savefig(fname, dpi=self.dpi)
