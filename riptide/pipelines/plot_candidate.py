##### Standard imports #####
import json

##### Non-standard imports #####
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from astropy.time import Time

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

    @property
    def axDMCurve(self):
        raise NotImplementedError

    @property
    def axPeriodCurve(self):
        raise NotImplementedError

    @property
    def axWidthCurve(self):
        raise NotImplementedError

    @property
    def axTable(self):
        raise NotImplementedError


class CandidatePlotDefaultLayout(CandidatePlotLayout):
    """ """
    default_figsize = (16, 5)
    default_dpi = 80

    def __init__(self, figsize=None, dpi=None):
        self.grid = gs.GridSpec(6, 10)
        self.figsize = self.default_figsize if figsize is None else figsize
        self.dpi = self.default_dpi if dpi is None else dpi

    @property
    def axSubints(self):
        return plt.subplot(self.grid[:4, 3:7])

    @property
    def axProfile(self):
        return plt.subplot(self.grid[4:, 3:7])

    @property
    def axDMCurve(self):
        return plt.subplot(self.grid[:3, 7:])

    @property
    def axPeriodCurve(self):
        return plt.subplot(self.grid[3:, 7:])

    # @property
    # def axWidthCurve(self):
    #     return plt.subplot(self.grid[4:, 7:])

    @property
    def axTable(self):
        return plt.subplot(self.grid[:, :3])


def plot_text(text, ax, x, y, horizontalalignment='right'):
    ax.text(x, y, text,
        horizontalalignment=horizontalalignment, verticalalignment='center',
        family='monospace', transform=ax.transAxes, fontsize=11)


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
        "name"  : 0.25,
        "value" : 0.87,
        "unit"  : 0.89
        }
    topMargin = 0.022
    lineHeight = 0.040

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
    """ Plot a Candidate object

    Parameters
    ----------
    cand: Candidate
        Candidate object to plot
    layout: str, optional
        Name of the plot layout (default: 'default')
    debase: bool, optional
        If True, subtract from the folded profile its median value before
        plotting it. This should result in an approximately zero off-pulse
        mean. (default: True)
    bar: bool, optional
        If True, plot the folded profile as a bar chart, otherwise plot
        it as a line. (default: True)
    figsize: tuple, optional
        Tuple (width, height) of the matplotlib.Figure
    dpi: float, optional
        dpi resultion of the matplotlib.Figure
    """
    _layouts = {
        'default': CandidatePlotDefaultLayout,
        }

    def __init__(self, cand, layout='default', debase=True, bar=True, figsize=None, dpi=None):
        self.cand = cand
        self.layout = None
        self.debase = False
        self.bar = bar
        self.figure = None
        self.shift_bins = 0
        self.layout = None

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
        plt.title("{:d} Sub-Integrations".format(self.cand.subints.nsubs))

    def _plotProfile(self):
        ax = self.layout.axProfile
        prof = np.roll(self.cand.subints.normalised_profile, self.shift_bins)
        if self.debase:
            prof -= np.median(prof)

        if self.bar:
            ax.bar(
                np.arange(0, prof.size),
                prof,
                width=1,
                color='#303030')
        else:
            ax.plot(
                np.arange(0, prof.size),
                prof,
                lw=1.0,
                color='#303030')
        plt.xlim(-0.5, prof.size-0.5)
        plt.ylabel("S/N")
        plt.title("Integrated Profile")

    def _plotDMCurve(self):
        ax = self.layout.axDMCurve
        cv = self.cand.dm_curve
        ax.plot(cv.trials, cv.snr, marker='o', markersize=3, color='#FF502A')

        # Call xlim() only if we have more than 1 DM trial, just to avoid a
        # warning
        if len(cv.trials) > 1:
            plt.xlim(cv.trials[0], cv.trials[-1])
        plt.grid(linestyle=':')
        plt.xlabel('DM (pc $\mathrm{cm}^{-3}$)')
        plt.ylabel("S/N")

    def _plotPeriodCurve(self):
        ax = self.layout.axPeriodCurve
        cv = self.cand.period_curve

        dp = cv.trials - self.cand.metadata['best_period']
        dp_max = abs(dp).max()
        if dp_max < 1e-3:
            factor, unit = 1e6, '$\mu$s'
        elif 1e-3 <= dp_max < 1.0:
            factor, unit = 1e3, 'ms'
        else:
            factor, unit = 1.0, 's'

        dp *= factor

        ax.plot(dp, cv.snr)
        plt.xlim(dp[0], dp[-1])
        plt.grid(linestyle=':')
        plt.xlabel('Delta Period ({:s})'.format(unit))
        plt.ylabel("S/N")

    def _plotWidthCurve(self):
        ax = self.layout.axWidthCurve
        cv = self.cand.width_curve

        width = self.cand.metadata['best_width']
        ducy = self.cand.metadata['best_ducy']

        # Duty cycle trials expressed in percent
        ducy_trials = cv.trials / (width / ducy) * 100.0

        ax.plot(ducy_trials, cv.snr, marker='o', markersize=3, color='#1E8449')
        plt.xlim(ducy_trials[0], ducy_trials[-1])
        plt.xscale('log')
        plt.grid(which='minor', linestyle=':')
        plt.grid(which='major', linestyle='-')
        plt.xlabel('Duty Cycle (%)')
        plt.ylabel("S/N")

    def _plotTable(self):
        table = Table()
        table.add_entry(Entry("Source", self.cand.metadata['source_name']))

        coord = self.cand.metadata['skycoord']
        ra = coord.ra.to_string(unit='hour', precision=2, pad=2)
        dec = coord.dec.to_string(unit='degree', precision=2, pad=2)
        glon_deg = coord.galactic.l.deg
        glat_deg = coord.galactic.b.deg


        table.add_entry(Entry("RA", ra))
        table.add_entry(Entry("Dec", dec))
        table.add_entry(Entry("glon", "{:.3f}".format(glon_deg), "deg"))
        table.add_entry(Entry("glat", "{:.3f}".format(glat_deg), "deg"))
        table.skip_line()

        mjd = self.cand.metadata['mjd']
        tobs = self.cand.metadata['tobs']
        # The precision=0 argument ensures that when printing in UTC format,
        # we do NOT get fractional seconds.
        time = Time(mjd, scale='utc', format='mjd', precision=0)
        utc = time.isot

        table.add_entry(Entry("MJD", "{:.6f}".format(mjd)))
        table.add_entry(Entry("UTC", utc))
        table.add_entry(Entry("Tobs", "{:.2f}".format(tobs), "s"))
        table.skip_line()

        period = self.cand.metadata['best_period']
        table.add_entry(Entry("Period", "{:.6f}".format(1000.0 * period), "ms"))

        dm = self.cand.metadata['best_dm']
        table.add_entry(Entry("DM", "{:.3f}".format(dm), "$\mathrm{pc}\,\mathrm{cm}^{-3}$"))

        ducy = self.cand.metadata['best_ducy']
        table.add_entry(Entry("Duty Cycle", "{:.1f}".format(ducy * 100.0), "%"))

        snr = self.cand.metadata['best_snr']
        table.add_entry(Entry("S/N", "{:.1f}".format(snr)))
        table.plot(self.layout.axTable)

    def _plot(self):
        self._plotProfile()
        self._plotSubints()
        self._plotDMCurve()
        self._plotPeriodCurve()
        #self._plotWidthCurve()
        self._plotTable()
        self.figure.tight_layout()

    def display(self):
        self.figure.show()

    def saveimg(self, fname):
        self.figure.savefig(fname, dpi=self.layout.dpi)