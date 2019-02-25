import numpy as np

def generate_width_trials(nbins, ducy_max=0.20, wtsp=1.5):
    widths = []
    w = 1
    wmax = int(max(1, ducy_max * nbins))
    while w <= wmax:
        widths.append(w)
        w = int(max(w + 1, wtsp * w))
    return np.asarray(widths)

def get_period_offset(shift, nprofs, nbins):
    """ Compute the period offset (in number of samples) from 'nbins'
    corresponding to a given top-to-bottom right shift in the FFA transform."""
    return (shift * nbins * 1.0) / (nprofs * nbins - shift)

def get_period(shift, nprofs, nbins):
    """ Compute the signal period (in number of samples) corresponding
    to a given top-to-bottom right shift in the FFA transform."""
    return get_period_offset(shift, nprofs, nbins) + nbins

def autocov(x, s):
    """ Autocovariance of time series, as a function of lag s.
    In particular, autocov(x,0) = var(x)"""
    n = x.size
    # NOTE: use float64 accumulator to avoid saturation issues when the
    # data have large values
    return np.multiply(x[:n-s], x[s:], dtype=np.float64).sum() / (n-s-1)
