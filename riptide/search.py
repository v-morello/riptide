import numpy as np

### Local module imports
from .libffa import LibFFA
from .processing_plan import ProcessingPlan
from .periodogram import Periodogram


def num_period_trials(plan):
    """ Compute the total number of period trials of a ProcessingPlan. This
    is required to allocate output arrays of the exact right size for the
    FFA C code."""
    num = 0
    for istep, step in plan.steps.iterrows():
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



def ffa_search(tseries, rmed_width=4.0, rmed_minpts=101, period_min=1.0, period_max=30.0, fpmin=8, bins_min=240, bins_max=260, ducy_max=0.20, wtsp=1.5, threads=1):
    """ Run a FFA search of a TimeSeries.

    Parameters:
    -----------
        tseries: TimeSeries
            The time series object to search
        [... and lots of kwargs]

    Returns:
    --------
        ts: TimeSeries
            The de-reddened and normalised time series that was actually searched
        plan: ProcessingPlan
            The processing plan followed
        pgram: Periodogram
            The output of the search, containing a 2D array representing S/N
            as a function of trial period and trial width.
    """
    ### Prepare data: deredden then normalise IN THAT ORDER
    ### 'ts' is the post-processed TimeSeries that will be searched
    ts = tseries.deredden(rmed_width, minpts=rmed_minpts)
    ts.normalise(inplace=True, correct_autocov=True)

    plan = ProcessingPlan.create(
        ts.nsamp, ts.tsamp,
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
        ts.data, ts.nsamp,
        dsfactors, bins_min, bins_max, dsfactors.size,
        plan.widths, plan.widths.size,
        threads,
        periods,
        snrs
        )

    # periods is in expressed in number of samples up to this point
    # NOTE: 'periods' does not come out perfectly sorted
    periods *= ts.tsamp

    # Reshape S/N trials to 2D
    snrs = snrs.reshape(periods.size, plan.widths.size)

    # Remove redundant period trials, i.e. those coresponding to excessively high shifts
    # in some FFA transforms
    snrs = snrs.reshape(periods.size, plan.widths.size)
    periods, snrs = prune_periodogram(periods, snrs)

    pgram = Periodogram(plan, periods, plan.widths, snrs, metadata=tseries.metadata)
    return ts, plan, pgram
