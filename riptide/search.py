import numpy as np

### Local module imports
from .libffa import LibFFA
from .processing_plan import ProcessingPlan
from .periodogram import Periodogram
from .timing import timing


def num_period_trials(plan):
    """ Compute the total number of period trials of a ProcessingPlan. This
    is required to allocate output arrays of the exact right size for the
    FFA C code."""
    num = 0
    for __, step in plan.steps.iterrows():
        ds = step.dsfactor
        ns = int(plan.nsamp / ds)
        b = np.arange(int(step.bins_min), int(step.bins_max))
        num += (ns // b).sum()
    return num


def prune_periodogram(periods, snrs):
    """ Remove redundant period trials in periodogram output arrays. """
    mask = np.ones(periods.size, dtype=bool)
    # When transforming an m x b input array with the FFA, shifts close to m may correspond
    # to periods slightly larger than b+1. Those period trials are therefore redundant, as they
    # will be re-performed when FFA transforming with b' = b + 1.

    # Reverse period trials array, makes the problem simpler:
    # we want to mask values in 'rp' to make it strictly decreasing.
    rp = periods[::-1]
    diff = np.diff(rp)

    # Indices such that the next element in 'rp' has a strictly larger period trial
    indices = np.where(diff > 0)[0]

    for ix in indices:
        for jj in range(ix+1, rp.size):
            mask[jj] = False
            if rp[jj] < rp[ix]:
                break

    mask = mask[::-1]
    return periods[mask], snrs[mask]


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
    plan : ProcessingPlan
        The processing plan followed
    pgram : Periodogram
        The output of the search, containing a 2D array representing S/N
        as a function of trial period and trial width.
    """
    ### Prepare data: deredden then normalise IN THAT ORDER
    if deredden:
        tseries = tseries.deredden(rmed_width, minpts=rmed_minpts)
    if not already_normalised:
        tseries = tseries.normalise()

    plan = ProcessingPlan.create(
        tseries.nsamp, tseries.tsamp,
        period_min=period_min, period_max=period_max, fpmin=fpmin,
        bins_min=bins_min, bins_max=bins_max, ducy_max=ducy_max, wtsp=wtsp
        )

    ### Prepare the call to the main C function
    npt = num_period_trials(plan)

    # Input arrays
    dsfactors = plan.steps.dsfactor.values
    bins_min = plan.steps.bins_min.values
    bins_max = plan.steps.bins_max.values

    # Allocate output arrays
    periods = np.zeros(npt, dtype=np.float32)
    snrs = np.zeros(npt * plan.widths.size, dtype=np.float32)

    #### The main C function that does all the work #####
    LibFFA.py_periodogram(
        tseries.data, tseries.nsamp,
        dsfactors, bins_min, bins_max, dsfactors.size,
        plan.widths, plan.widths.size,
        periods,
        snrs
        )

    # periods is in expressed in number of samples up to this point
    # NOTE: 'periods' does not come out perfectly sorted
    periods *= tseries.tsamp

    # Reshape S/N trials to 2D
    snrs = snrs.reshape(periods.size, plan.widths.size)

    # Remove redundant period trials, i.e. those coresponding to excessively high shifts
    # in some FFA transforms
    snrs = snrs.reshape(periods.size, plan.widths.size)
    periods, snrs = prune_periodogram(periods, snrs)

    pgram = Periodogram(plan, periods, plan.widths, snrs, metadata=tseries.metadata)
    return tseries, plan, pgram
