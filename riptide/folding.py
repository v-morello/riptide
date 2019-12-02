import numpy as np

from riptide.libffa import downsample


def downsample_vertical(X, factor):
    m, __ = X.shape

    if not factor > 1:
        raise ValueError("factor must be > 1")
    if not factor < m:
        raise ValueError("factor must be strictly smaller than the number of input lines")

    Y = np.ascontiguousarray(X.T)
    out = np.asarray([downsample(arr, factor) for arr in Y])
    return np.ascontiguousarray(out.T)


def fold(ts, period, bins, subints=None):
    """
    Fold TimeSeries at given period

    Parameters
    ----------
    ts : TimeSeries
        Input time series to fold
    period : float
        Period in seconds
    bins : int
        Number of phase bins
    subints : int or None, optional
        Number of desired sub-integrations. If None, the number of 
        sub-integrations will be the number of full periods that fit in 
        the data (default: None)

    Returns
    -------
    folded : ndarray
        The folded data as a numpy array. If subints > 1, it has a shape
        (subints, bins). Otherwise it is a 1D array with 'bins' elements.
    
    Raises
    ------
    ValueError: if the data cannot be folded with the requested parameters,
    e.g. bin width is shorter than sampling time, or subint length is shorter
    than requested period
    """
    if period > ts.length:
        raise ValueError("Period exceeds data length")

    tbin = period / bins
    if not tbin > ts.tsamp:
        raise ValueError("Bin width is shorter than sampling time")

    if subints is not None:
        subints = int(subints)
        if not subints >= 1:
            raise ValueError("subints must be >= 1 or None")

        full_periods = ts.length / period
        if subints > full_periods:
            raise ValueError(f"subints ({subints}) exceeds the number of signal periods that fit in the data ({full_periods})")

    factor = tbin / ts.tsamp
    tsdown = ts.downsample(factor)
    m = tsdown.nsamp // bins
    nsamp_eff = m * bins
    
    folded = tsdown.data[:nsamp_eff].reshape(m, bins)
    folded *= (m * factor) ** -0.5

    if subints == 1 or m == 1:
        return folded.sum(axis=0)
    elif subints is None:
        return folded
    elif subints == m:
        return folded
    else:
        # vertical downsampling factor
        vf = m / subints
        return downsample_vertical(folded, vf)
