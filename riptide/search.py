import numpy as np

### Local module imports
import riptide.libcpp as libcpp
from .ffautils import generate_width_trials
from .periodogram import Periodogram
from .timing import timing


@timing
def ffa_search(tseries, period_min=1.0, period_max=30.0, fpmin=8, bins_min=240, bins_max=260,
    ducy_max=0.20, wtsp=1.5, deredden=True, rmed_width=4.0, rmed_minpts=101, already_normalised=False):
    """ 
    Run a FFA search of a single TimeSeries object, producing its periodogram.

    Parameters
    ----------
    tseries : TimeSeries
        The time series object to search
    period_min : float
        Minimum period to search in seconds
    period_max : float
        Maximum period to search in seconds
    fpmin : int
        Minimum number of signal periods that must fit in the data. In other
        words, place a cap on period_max equal to DATA_LENGTH / fpmin
    bins_min : int
        Minimum number of phase bins in the folded data. Higher values
        provide better duty cycle resolution. As the code searches longer trial
        periods, the data are iteratively downsampled so that the number of 
        phase bins remains between bins_min and bins_max
    bins_max : int
        Maximum number of phase bins in the folded data. Must be strictly
        larger than bins_min, approx. 10% larger is a good choice
    wtsp : float
        Multiplicative factor between consecutive boxcar width trials. The
        smallest width is always 1 phase bin, and the sequence of width
        trials is generated with the formula:
        W(n+1) = max( floor(wtsp x W(n)), W(n) + 1 )
        wtsp = 1.5 gives the following sequence of width trials (in number of
        phase bins): 1, 2, 3, 4, 6, 9, 13, 19 ...
    ducy_max : float
        Maximum duty cycle to optimally search. Limits the maximum width of the
        boxcar matched filters applied to any given profile. 
        Example: on a 300 phase bin profile, ducy_max = 0.2 means that no 
        boxcar filter of width > 60 bins will be applied
    deredden : bool
        Subtract red noise from the time series before searching
    rmed_width : float
        The width of the running median filter to subtract from the input data
        before processing, in seconds 
    rmed_minpts : int
        The running median is calculated of a time scrunched version of the
        input data to save time: rmed_minpts is the minimum number of
        scrunched samples that must fit in the running median window
        Lower values make the running median calculation less accurate but
        faster, due to allowing a higher scrunching factor
    already_normalised : bool
        Assume that the data are already normalised to zero mean and unit 
        standard deviation

    Returns
    -------
    ts : TimeSeries
        The de-reddened and normalised time series that was actually searched
    pgram : Periodogram
        The output of the search, which contains among other thigs a 2D array 
        representing S/N as a function of trial period and trial width.
    """
    ### Prepare data: deredden then normalise IN THAT ORDER
    if deredden:
        tseries = tseries.deredden(rmed_width, minpts=rmed_minpts)
    if not already_normalised:
        tseries = tseries.normalise()

    widths = generate_width_trials(bins_min, ducy_max=ducy_max, wtsp=wtsp)
    periods, foldbins, snrs = libcpp.periodogram(
        tseries.data, tseries.tsamp, widths, 
        period_min, period_max, bins_min, bins_max
        )
    pgram = Periodogram(widths, periods, foldbins, snrs, metadata=tseries.metadata)
    return tseries, pgram
