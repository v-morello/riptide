import logging
from collections import namedtuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import astropy.units as uu
from astropy.time import Time


log = logging.getLogger('riptide.candidate')


class Candidate(object):
    """
    Final data product of the riptide pipeline

    Attributes
    ----------
    params : dict
        Dictionary with best-fit parameters of the signal:
        period, freq, dm, width, ducy, snr

    tsmeta : Metadata
        Metadata of the TimeSeries object (DM trial) in which the Candidate was
        found to have the highest S/N, and from which it was folded

    peaks : pandas.DataFrame
        A pandas DataFrame with the attributes of the periodogram peaks
        associated to the Candidate

    subints : ndarray
        A two-dimensional numpy array with shape (num_subints, num_bins)
        containing the folded sub-integrations

    profile : ndarray
        Folded profile as a one-dimensional numpy array, normalised such that
        the background noise standard deviation is 1, and the mean of the
        profile is zero

    dm_curve : tuple
        Tuple of numpy arrays (dm, snr) containing respectively the sequence
        of DM trials, and corresponding best S/N value across all trial widths
    """
    def __init__(self, params, tsmeta, peaks, subints):
        self.params = params
        self.tsmeta = tsmeta
        self.peaks = peaks
        self.subints = subints

    def to_dict(self):
        """ Convert to dictionary for serialization """
        return {
            'params': self.params,
            'tsmeta': self.tsmeta,
            'peaks': self.peaks,
            'subints': self.subints
        }

    @property
    def profile(self):
        if self.subints.ndim == 1:
            return self.subints
        return self.subints.sum(axis=0)

    @property
    def dm_curve(self):
        # NOTE: copy() works around a bug in pandas 0.23.x and earlier
        # https://stackoverflow.com/questions/53985535/pandas-valueerror-buffer-source-array-is-read-only
        # TODO: consider requiring pandas 0.24+ in the future
        df = self.peaks.copy().groupby('dm').max()
        return df.index.values, df.snr.values

    @classmethod
    def from_pipeline_output(cls, ts, peak_cluster, bins, subints=1):
        """
        Method used by the pipeline to produce a candidate from intermediate
        data products. 
        
        subints can be an int or None. None means pick the number of subints
        that fit inside the data.

        If 'subints' is too large to fit in the data, then
        this function will call TimeSeries.fold() with subints=None.
        """
        centre = peak_cluster.centre
        P0 = centre.period

        if subints is not None and subints * P0 >= ts.length:
            msg = (
                f"Period ({P0:.3f}) x requested subints ({subints:d}) exceeds time series length "
                f"({ts.length:.3f}), setting subints = full periods that fit in the data"
            )
            log.debug(msg)
            subints = None

        subints_array = ts.fold(centre.period, bins, subints=subints)
        return cls(centre.summary_dict(), ts.metadata, peak_cluster.summary_dataframe(), subints_array)

    @classmethod
    def from_dict(cls, items):
        """ De-serialize from dictionary """
        return cls(items['params'], items['tsmeta'], items['peaks'], items['subints'])

    def plot(self, figsize=(18, 4.5), dpi=80):
        """
        Create a plot of the candidate

        Parameters
        ----------
        figsize : tuple
            figsize argument passed to plt.figure()
        dpi : int
            dpi argument passed to plt.figure()

        Returns
        -------
        fig : matplotlib.Figure
        """
        fig = plt.figure(figsize=figsize, dpi=dpi)
        plot_candidate(self)
        return fig

    def show(self, **kwargs):
        """
        Create a plot of the candidate and display it. Accepts the same keyword
        arguments as plot().
        """
        self.plot(**kwargs)
        plt.show()

    def savefig(self, fname, **kwargs):
        """
        Create a plot of the candidate and save it as PNG under the specified 
        file name. Accepts the same keyword arguments as plot().
        """
        self.plot(**kwargs)
        plt.savefig(fname)
        plt.close()

    def __str__(self):
        name = type(self).__name__
        return f"{name}({self.params})"

    def __repr__(self):
        return str(self)


TableEntryBase = namedtuple('TableEntry', ['name', 'value', 'formatter', 'unit'])


class TableEntry(TableEntryBase):
    def plot(self, X, y, **kwargs):
        """
        X : list
            list of X coordinates for each column
        y : float
            Y coordinate of the line
        """
        assert(len(X) == 3)
        fmt = "{{:{}}}".format(self.formatter)
        plt.text(X[0], y, self.name, **kwargs)
        plt.text(X[1], y, fmt.format(self.value), ha='right', **kwargs)
        plt.text(X[2], y, self.unit, **kwargs)


def plot_table(params, tsmeta):
    """
    """
    plt.axis('off')
    coord = tsmeta['skycoord']
    ra_hms = coord.ra.to_string(unit=uu.hour, sep=':', precision=2, pad=True)
    dec_hms = coord.dec.to_string(unit=uu.deg, sep=':', precision=2, pad=True)

    # TODO: Check that the scale is actually UTC in the general case
    # PRESTO, SIGPROC and other packages may not have the same
    # date/time standard
    obsdate = Time(tsmeta['mjd'], format='mjd', scale='utc', precision=0)

    blank = TableEntry(name='', value='', formatter='s', unit='')

    entries = [
        TableEntry(name='Period', value=params['period'] * 1000.0, formatter='.6f', unit='ms'),
        TableEntry(name='DM', value=params['dm'], formatter='.2f', unit='pc cm$^{-3}$'),
        TableEntry(name='Width', value=params['width'], formatter='d', unit='bins'),
        TableEntry(name='Duty cycle', value=params['ducy'] * 100.0, formatter='.2f', unit='%'),
        TableEntry(name='S/N', value=params['snr'], formatter='.1f', unit=''),
        blank,
        TableEntry(name='Source', value=tsmeta['source_name'], formatter='s', unit=''),
        TableEntry(name='RA', value=ra_hms, formatter='s', unit=''),
        TableEntry(name='Dec', value=dec_hms, formatter='s', unit=''),
        TableEntry(name='MJD', value=obsdate.mjd, formatter='.6f', unit=''),
        TableEntry(name='UTC', value=obsdate.iso, formatter='s', unit=''),
    ]

    y0 = 0.94  # Y coordinate of first line
    dy = 0.105  # line height
    X = [0.0, 0.80, 0.84] # Coordinate of columns name, value, unit

    for ii, entry in enumerate(entries):
        entry.plot(X, y0 - ii * dy, family='monospace')


def plot_dm_curve(dm, snr):
    dm_min = dm.min()
    dm_max = dm.max()
    plt.plot(dm, snr, color='r', marker='o', markersize=3)

    # Avoid matplotlib warning when calling xlim() with two equal values 
    if dm_min == dm_max:
        plt.xlim(dm_min - 0.5, dm_min + 0.5)
    else:
        plt.xlim(dm_min, dm_max)
    plt.grid(linestyle=':')
    plt.xlabel("DM (pc cm$^{-3}$)")
    plt.ylabel("Best S/N")


def plot_subints(X, T):
    """
    X : ndarray
        Sub-integrations array, shape = (nsubs, nbins)
    T : float
        Integration time in seconds
    """
    __, nbins = X.shape

    X = np.hstack((X, X[:, :nbins//2]))
    __, nbins_ext = X.shape

    plt.imshow(
        X, 
        cmap='Greys', interpolation='nearest', aspect='auto',
        extent=[-0.5, nbins_ext-0.5, T, 0] # Note: t = 0 is at the top of the plot
        )
    plt.fill_between([nbins, nbins_ext], [0, 0], [T, T], color='b', alpha=0.08)
    plt.xlim(-0.5, nbins_ext-0.5)
    plt.ylabel("Time (seconds)")
    plt.title("1.5 Periods of Signal")


def plot_profile(P):
    """
    P : profile normalised to unit background noise variance
    """
    nbins = len(P)
    P = np.concatenate((P, P[:nbins//2]))
    nbins_ext = len(P)

    plt.bar(range(nbins_ext), P - np.median(P), width=1, color='#404040')

    ymin, ymax = plt.ylim()
    plt.fill_between([nbins, nbins_ext], [ymin, ymin], [ymax, ymax], color='b', alpha=0.08)
    plt.ylim(ymin, ymax)

    plt.xlim(-0.5, nbins_ext-0.5)
    plt.xlabel("Phase bin")
    plt.ylabel("Normalised amplitude")


def plot_candidate(cand):
    """
    Plot candidate on the current figure
    """
    # https://matplotlib.org/tutorials/intermediate/gridspec.html
    nrows, ncols = 2, 7
    gs = GridSpec(nrows, ncols, figure=plt.gcf())

    plt.subplot(gs[:1, 2:])
    plot_subints(cand.subints, cand.tsmeta['tobs'])

    plt.subplot(gs[1:, 2:])
    plot_profile(cand.subints.sum(axis=0))

    plt.subplot(gs[:1, :2])
    plot_table(cand.params, cand.tsmeta)

    plt.subplot(gs[1:, :2])
    plot_dm_curve(*cand.dm_curve)

    plt.tight_layout()


if __name__ == '__main__':
    from riptide import load_json
    c = load_json("/home/vince/work/time_series/J1119-7936/candidate_0004.json")
    c.show()

    c = load_json("/home/vince/work/time_series/J1855+0307/candidate_0000.json")
    c.show()

    plt.show()