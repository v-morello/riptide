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
        ns = int(plan.nsamp_ds / ds)
        b = np.arange(int(step.bins_min), int(step.bins_max))
        num += (ns // b).sum()
    return num


def ffa_search(tseries, rmed_width=4.0, rmed_minpts=101, period_min=1.0, period_max=30.0, fpmin=8, bins_min=240, bins_max=260, ducy_max=0.20, wtsp=1.5, threads=1):
    """ Run a FFA search of a TimeSeries. 
    
    Parameters:
    -----------

    Returns:
    --------
    """
    plan = ProcessingPlan(
        tseries.nsamp, tseries.tsamp, 
        period_min=period_min, period_max=period_max, fpmin=fpmin, 
        bins_min=bins_min, bins_max=bins_max, ducy_max=ducy_max, wtsp=wtsp
        )
    
    ### Prepare data: downsample, deredden, normalise IN THAT ORDER
    ### 'ts' is the post-processed TimeSeries that will be searched
    ts = tseries.downsample(plan.ds_ini, inplace=False)
    ts.deredden(rmed_width, minpts=rmed_minpts, inplace=True)
    ts.normalise(inplace=True)
    
    ### Prepare the call to the main C function
    npt = num_period_trials(plan)
    
    # Input arrays 
    dsfactors = plan.steps.dsfactor.as_matrix()
    bins_min = plan.steps.bins_min.as_matrix()
    bins_max = plan.steps.bins_max.as_matrix()

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
    
    pgram = Periodogram(periods, plan.widths, snrs)
    return ts, plan, pgram
