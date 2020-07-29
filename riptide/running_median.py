import numpy as np
from riptide.libcpp import running_median


def time_scrunch(data, factor):
    factor = int(factor)
    if factor <= 1:
        return data
    N = (data.size // factor) * factor
    return data[:N].reshape(-1, factor).mean(axis=1)


def fast_running_median(data, width_samples, min_points=101):
    """
    Compute an approximate running median of data over large window sizes.

    Parameters
    ----------
    data : ndarray
        Input data
    width : int
        Required width of the running median window in number of samples
    min_points : int
        The running median is calculated of a time scrunched version of the
        input data to save time: min_points is the minimum number of
        scrunched samples that must fit in the running median window.
        Lower values make the running median calculation less accurate but
        faster, due to allowing a higher scrunching factor.
    """
    scrunch_factor = int(max(1, width_samples / float(min_points)))
    scrunched_data = time_scrunch(data, scrunch_factor)
    rmed_lores = running_median(scrunched_data, min_points)
    x_lores = np.arange(scrunched_data.size) * scrunch_factor + 0.5 * (scrunch_factor - 1)
    return np.interp(np.arange(data.size), x_lores, rmed_lores)
