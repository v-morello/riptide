import numpy as np
import riptide.libcpp


def running_median(x, width_samples):
    """
    Computes the running median of data with the specified window size.

    Parameters
    ----------
    x : ndarray
        One dimensional input data.
    width_samples : int
        The width of the running median window in number of elements.
        It must be an odd number smaller than the input data length, otherwise
        ValueError is raised.

    Returns
    -------
    rmed : ndarray
        The running median of 'x'.

    Notes
    -----
    The C++ running median code internally pads both ends of the arrray with
    the edge values.

    If the input array is not contiguous in memory, a temporary contiguous copy
    is made and passed to the C++ function (which only accepts C-contiguous
    arrays). Otherwise no performance hit is incurred.

    See Also
    --------
    fast_running_median : an approximate running median that runs much faster
        with large window sizes (> 100 elements).
    """
    return riptide.libcpp.running_median(np.ascontiguousarray(x), width_samples)


def scrunch(data, factor):
    """
    Reduce the resolution of data by adding consecutive elements together
    """
    factor = int(factor)
    N = (data.size // factor) * factor
    return data[:N].reshape(-1, factor).mean(axis=1)


def fast_running_median(data, width_samples, min_points=101):
    """
    Compute an approximate running median of data over large window sizes. The
    idea is to downsample the data (if necessary), call running_median() on it
    and linearly interpolate it back to the original resolution.

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
        NOTE: 'min_points' must be an odd number.

    See Also
    --------
    running_median : an exact running median but slower for large window sizes
    """
    if not (min_points % 2):
        raise ValueError("min_points must be an odd number")
    scrunch_factor = int(max(1, width_samples / float(min_points)))

    if scrunch_factor == 1:
        return running_median(data, width_samples)

    scrunched_data = scrunch(data, scrunch_factor)
    rmed_lores = running_median(scrunched_data, min_points)
    x_lores = np.arange(scrunched_data.size) * scrunch_factor + 0.5 * (
        scrunch_factor - 1
    )
    return np.interp(np.arange(data.size), x_lores, rmed_lores)
