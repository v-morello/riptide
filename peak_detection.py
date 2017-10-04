import warnings

##### Non-standard imports #####
import numpy as np
import matplotlib.pyplot as plt
import pandas

# NOTE: the spacing of FFA period trials in frequency space is nearly constant
# Variations in spacing are mostly caused by the discrete downsampling steps

# NOTE: one segment should cover at the very least 1 DFT bin worth of trials.
# This is to ensure some level of
# statistical independence between consecutive bins
# 10 DFT bins or more seems much better
def segment(periods, tobs, segment_dftbins_length=10.0):
    """ Break the period trials array into equal-sized segments, returning
    their boundaries and other information as a pandas.DataFrame
    """
    # Number of DFT bin indexes spanned by the period trials
    dbi_range = tobs / periods[0] - tobs / periods[-1]
    nseg = int(dbi_range / segment_dftbins_length)
    slen = len(periods) // nseg

    boundaries = [
        (iseg * slen, (iseg + 1) * slen)
        for iseg in range(nseg)
        ]

    # Last segment must end at the last period trial
    last = boundaries[-1]
    boundaries[-1] = (last[0], len(periods))
    boundaries = np.asarray(boundaries)

    data = pandas.DataFrame(
        columns=['istart', 'iend', 'imid', 'pmid', 'logpmid']
        )
    data['istart'] = boundaries[:, 0]
    data['iend'] = boundaries[:, 1]
    data['imid'] = (data['istart'] + data['iend']) // 2
    data['pmid'] = periods[data['imid']]
    data['logpmid'] = np.log(data['pmid'])
    return data

def segment_stats(snrs, boundaries):
    """ Compute the median and robust standard deviation of each Periodogram
    segment. """
    percentiles = (25, 50, 75)
    stats = pandas.DataFrame(
        [ np.percentile(snrs[s:e], percentiles) for s, e in boundaries[['istart', 'iend']].as_matrix() ],
        columns=['p25', 'median', 'p75']
        )
    stats['sigma'] = (stats['p75'] - stats['p25']) / 1.349
    stats['pmid'] = boundaries['pmid']
    stats['logpmid'] = boundaries['logpmid']
    return stats


# NOTE: threshold will be a 2nd degree polynomial function of x = log(period)
# This function returns the polynomial coefficients in x
def threshold_function_dynamic(stats, snr_min=6.5, nsigma=6.5, polydeg=2):
    x = stats['logpmid']
    y = stats['median'] + nsigma * stats['sigma']
    poly = np.poly1d( np.polyfit(x, y, polydeg) )
    def func(period):
        return np.maximum(poly(np.log(period)), snr_min)
    return func, poly.coefficients

def threshold_function_static(snr_min, polydeg=2):
    def func(period):
        thr = np.empty(len(period))
        thr.fill(snr_min)
        return thr
    polyco = np.zeros(polydeg + 1)
    return func, polyco


class Peak(object):
    """ """
    def __init__(self, periods, snrs, width, dm):
        self._periods = periods
        self._snrs = snrs
        imax = snrs.argmax()
        self._best_period = periods[imax]
        self._best_snr = snrs[imax]
        self._width = width
        self._dm = dm

    @property
    def periods(self):
        return self._periods

    @property
    def snrs(self):
        return self._snrs

    @property
    def width(self):
        return self._width

    @property
    def dm(self):
        return self._dm

    @property
    def best_period(self):
        return self._best_period

    @property
    def best_snr(self):
        return self._best_snr

    def plot(self):
        delta_period_us = 1.0e6 * (self.periods - self.best_period)
        best_period_ms = 1.0e3 * self.best_period
        plt.plot(delta_period_us, self.snrs, marker='.', markersize=4)
        plt.title('P0 = {0:.6f} ms, Width = {1:d} bins, DM = {2:.3f}'.format(best_period_ms, self.width, self.dm))
        plt.xlabel('Delta P0 (us)')
        plt.ylabel('S/N')
        plt.grid(linestyle=':')
        plt.tight_layout()

    def display(self, figsize=(6, 4), dpi=100):
        plt.figure(figsize=figsize, dpi=dpi)
        self.plot()
        plt.show()

    def __str__(self):
        dm_str = 'None' if self.dm is None else '{0:.3f}'.format(self.dm)
        return 'Peak [P0 = {p.best_period:.9e}, W = {p.width:3d}, DM = {dm_str:s}, S/N = {p.best_snr:6.2f}]'.format(p=self, dm_str=dm_str)

    def __repr__(self):
        return str(self)



# TODO: docstring for this
def iterslices(indices):
    # Special case where there is only one peak
    if not len(indices):
        yield slice(None, None)
        raise StopIteration

    # Multiple peak case
    else:
        yield slice(None, indices[0])
        for ii, jj in zip(indices[:-1], indices[1:]):
            yield slice(ii, jj)
        yield slice(indices[-1], None)
        raise StopIteration


def find_peaks_single(pgram, iwidth, boundaries, min_segments=8, snr_min=6.5, nsigma=6.5, peak_clustering_radius=0.20, polydeg=2):
    periods = pgram.periods
    snrs = pgram.snrs[:, iwidth]
    width = pgram.widths[iwidth]
    tobs = pgram.tobs
    dm = pgram.metadata['dm']

    stats = segment_stats(snrs, boundaries)

    # NOTE: tfunc is a function of period only. tfunc has one parameter expected
    # to be a numpy array
    if len(boundaries) < min_segments:
        #warnings.warn('find_peaks_single(): not enough segments for dynamic threshold fitting, applying constant S/N threshold')
        tfunc, polyco = threshold_function_static(snr_min=snr_min, polydeg=polydeg)
    else:
        tfunc, polyco = threshold_function_dynamic(stats, snr_min=snr_min, nsigma=nsigma, polydeg=polydeg)

    # This is our selection threshold as a function of *period*
    # (NOT log(period) this time)
    threshold = tfunc(periods)

    # Now find peaks: sets of period trials that lie above the threshold
    # Two points whose dft bin indexes lie within 'peak_clustering_radius'
    # of each other are considered part of the same peak
    peaks = []

    significant_mask = snrs > threshold
    significant_indices = np.where(significant_mask)[0]
    if len(significant_indices):
        # DFT bin indexes corresponding to the significant periodogram points
        dbi = tobs / periods[significant_mask]

        # Whenever the difference between two consecutive elements of 'dbi'
        # exceeds 'peak_clustering_radius', this marks the end of a peak and the
        # start of a new one
        dbi_diff = np.diff(dbi)
        diffmask = np.abs(dbi_diff) > peak_clustering_radius
        ibreaks = np.where(diffmask)[0] + 1

        for sl in iterslices(ibreaks):
            peak_indices = significant_indices[sl]
            current_peak = Peak(
                periods[peak_indices],
                snrs[peak_indices],
                width,
                dm
                )
            peaks.append(current_peak)

    return stats, polyco, threshold, peaks


def find_peaks(pgram, segment_dftbins_length=10.0, min_segments=8, snr_min=6.5, nsigma=6.5, peak_clustering_radius=0.20, polydeg=2):
    ### Cut period trials into segments
    boundaries = segment(pgram.periods, pgram.tobs, segment_dftbins_length=segment_dftbins_length)

    all_peaks = []
    polyco_tracker = {}
    stats_tracker = {}

    for iwidth, width in enumerate(pgram.widths, start=0):
        stats, polyco, threshold, peaks = find_peaks_single(
            pgram,
            iwidth,
            boundaries,
            min_segments=min_segments,
            snr_min=snr_min,
            nsigma=nsigma,
            peak_clustering_radius=peak_clustering_radius,
            polydeg=polydeg
            )
        all_peaks = all_peaks + peaks
        polyco_tracker[width] = polyco
        stats_tracker[width] = stats
    return all_peaks, stats_tracker, polyco_tracker
