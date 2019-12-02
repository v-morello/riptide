import logging
import typing
from math import ceil

import numpy as np

from riptide.clustering import cluster1d
from riptide.timing import timing


log = logging.getLogger('riptide.peak_detection')


class Peak(typing.NamedTuple):
    """ 
    A simple NamedTuple with the essential parameters of a peak found 
    in a Periodogram
    """
    period: float
    freq: float
    width: int
    ducy: float  # duty cycle = width / average_folding_bins
    iw: int      # width trial index
    ip: int      # period trial index
    snr: float
    dm: float

    def summary_dict(self):
        """ 
        Returns a minimal dictionary of attributes to be written as CSV
        by the pipeline 
        """
        attrs = ('period', 'freq', 'dm', 'width', 'ducy', 'snr')
        return {a: getattr(self, a) for a in attrs}


def segment_stats(f, s, T, segwidth=5.0):
    """
    Cut a periodogram in consecutive, equal-sized segments with a 
    frequency span equal to segwidth / T, and return the centre frequencies,
    median S/N and robust S/N standard deviation of all segments.
    
    This information is then used to fit a sensible peak selection threshold
    as a function of frequency.

    Parameters
    -----------
    f : ndarray
        Frequencies in Hz
    S : ndarray
        Signal-to-noise ratios for a single width trial
    T : float
        Integration time in seconds
    segwidth : float
        Frequency segment width expressed in units of 1/T

    Returns
    -------
    fc : ndarray
        Median frequency of the segments
    smed : ndarray
        Median S/N of the segments
    sstd : ndarray
        S/N standard deviation of the segments, measured from the interquartile
        range of the segment's S/N distribution (stddev = IQR / 1.349)
    """
    w = segwidth / T
    #log.debug("Segment width (Hz): {:.6f}".format(w))

    # NOTE: the spacing of frequency trials is almost constant
    m = ceil(abs(f[-1] - f[0]) / w) # number of segments
    #log.debug("Segments: {:d}".format(m))

    p = len(f) // m # number of complete segments
    #log.debug("Points/segment: {:d}".format(p))

    n = m * p # effective number of elements
    f = f[:n]
    s = s[:n]

    fc = np.median(f.reshape(m, p), axis=1)
    s25, smed, s75 = np.percentile(s.reshape(m, p), (25, 50, 75), axis=-1)
    sstd = (s75 - s25) / 1.349
    return fc, smed, sstd


def fit_threshold(fc, tc, polydeg=2):
    """
    Fit a polynomial in log(f) to the selection threshold control points
    (fc, tc)

    Parameters
    ----------
    fc : ndarray
        Frequency of the control points
    tc : ndarray
        Value of the selection threshold at the control frequencies
    polydeg : int
        Degree of the log(f) polynomial to fit

    Returns
    -------
    poly : np.poly1d
        Polynomial in log(f) that represents the selection threshold as a
        function of frequency
    """
    coeffs = np.polyfit(np.log(fc), tc, polydeg)
    return np.poly1d(coeffs)


def find_peaks_single(f, s, T, smin=6.0, segwidth=5.0, nstd=7.0, minseg=10, polydeg=2, clrad=0.1):
    """
    Find peaks in a single pulse width trial. Returns a list of array indices
    that correspond to peak centres
    """
    peak_indices = []

    # Control points
    fc, smed, sstd = segment_stats(f, s, T, segwidth=segwidth)
    sc = smed + nstd * sstd

    # Selection threshold: polynomial in log(f)
    if len(fc) >= minseg:
        poly = fit_threshold(fc, sc, polydeg=polydeg)
        polyco = poly.coefficients
    else: # constant threshold if not enough points for fit
        polyco = [smin]
        poly = np.poly1d(polyco)

    # Selected frequencies and frequency indices
    dynthr = poly(np.log(f))
    mask = (s > dynthr) & (s > smin)
    indices = np.where(mask)[0]
    fsel = f[indices]

    clusters = cluster1d(fsel, clrad / T)
    for cl in clusters:
        ix = indices[cl]
        ipeak = s[ix].argmax()
        ipeak = ix[ipeak]
        peak_indices.append(ipeak)
    return peak_indices, polyco


@timing
def find_peaks(pgram, smin=6.0, segwidth=5.0, nstd=6.0, minseg=10, polydeg=2, clrad=0.1):
    """
    Identify significant peaks in a periodogram using a dynamically fitted
    S/N selection threshold. The fitting involves the following procedure for
    each pulse width trial separately:

    1. Cut the frequency range covered by the periodogram in segments
       of length 1 / T_obs
    2. Get the median S/N 'm' and robust S/N standard deviation 's' of each 
       segment. The dynamic selection threshold for that segment should be
       t = m + nstd x s
    3. Fit a polynomial in log(f) to the control points (f_i, t_i) thus 
       obtained
    4. Any point whose S/N exceeds both the dynamic threshold and the value 
       'smin' are considered significant
    5. Cluster these points. Two points are in the same peak if their
       trial frequencies are within clrad / T_obs of each other. All such
       clusters constitute a Peak.
    
    Parameters
    ----------
    pgram : Periodogram
        Input periodogram to search for peaks
    smin : float, optional
        Minimum S/N that a peak must exceed, in addition to having to exceed 
        the dynamic selection threshold
    segwidth : float, optional
        Width of a frequency segment in units of 1 / T_obs
    nstd : float, optional
        See above for an explanation
    minseg : float, optional
        Minimum number of segments below which only a static selection
        threshold is applied
    polydeg : float, optional
        Degree of polynomial in log(f)
    clrad : float, optional
        Clustering radius in frequency space, in units of 1 / T_obs

    Returns
    -------
    peaks: list
        List of Peak objects
    polycos: dict
        Dictionary of polynomial coefficients {iw: polyco}
        where 'iw' is the width trial index, and 'polyco' a list of polynomial
        coefficients in log(f). These can be passed to np.poly1d
    """
    f = pgram.freqs
    T = pgram.tobs
    dm = pgram.metadata['dm']
    bins_avg = pgram.bins_avg

    peaks = []
    polycos = {}
    for iw, width in enumerate(pgram.widths):
        s = pgram.snrs[:, iw].astype(float)
        cur_peak_indices, cur_polycos = find_peaks_single(
            f, s, T, 
            smin=smin, segwidth=segwidth, nstd=nstd, minseg=minseg, polydeg=polydeg, clrad=clrad
        )
        for ipeak in cur_peak_indices:
            peak_freq = f[ipeak]

            # NOTE: type enforcement is necessary, otherwise some Peak members
            # have np.float32 type which causes trouble down the line
            # NOTE 2: dm can be None on fake time series
            peak = Peak(
                freq=float(peak_freq), period=float(1.0/peak_freq), width=int(width), 
                ducy=float(width/bins_avg), iw=int(iw), ip=int(ipeak), snr=float(s[ipeak]), 
                dm=dm)
            #log.debug(peak)
            peaks.append(peak)
        polycos[iw] = cur_polycos
    return peaks, polycos
