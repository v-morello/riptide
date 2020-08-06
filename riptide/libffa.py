### Standard library imports
import os
import ctypes

### Non-standard imports
import numpy as np
import numpy.ctypeslib as npct
from numpy import log, sin, cos, exp, pi

### Local imports
from .ffautils import generate_width_trials
import riptide.libcpp as libcpp


def generate_signal(nsamp, period, phi0=0.5, ducy=0.02, amplitude=10.0, stdnoise=1.0):
    """ Generate a time series containing a periodic signal with a von Mises
    pulse profile. This function is useful for test purposes.

    Parameters
    ----------
    nsamp : int
        Number of samples to generate.
    period : float
        Period in number of samples.
    phi0 : float, optional
        Initial pulse phase in number of periods.
    ducy : float, optional
        Duty cycle of the pulse, i.e. the ratio FWHM / Period
    amplitude : float, optional
        True amplitude of the signal as defined in the reference paper.
        The *expectation* of the S/N of the generated signal is
            S/N_true = amplitude / stdnoise,
        assuming that a matched filter with the exact shape of the pulse is 
        employed to measure S/N (here: von Mises with given duty cycle). 
        riptide employs boxcar filters in the search, which results in a slight
        S/N loss. See the reference paper for details. 

        A further degradation will be observed on bright signals, because
        they bias the estimation of the mean and standard deviation of the 
        noise in a blind search.
    stdnoise : float, optional
        Standard deviation of the background noise. If set to 0, a noiseless 
        signal is generated.

    Returns
    -------
    tseries : ndarray (1D, float)
        Output time series.
    """
    # von mises parameter
    kappa = log(2.0) / (2.0 * sin(pi*ducy/2.0)**2)

    # Generate pulse train
    phase_radians = (np.arange(nsamp, dtype=float) / period - phi0) * (2 * pi)
    signal = exp(kappa*(cos(phase_radians) - 1.0))

    # Normalise to unit L2-norm, then scale by amplitude
    scale_factor = amplitude * (signal ** 2).sum() ** -0.5
    signal *= scale_factor

    # Add noise
    if stdnoise > 0.0:
        noise = np.random.normal(size=nsamp, loc=0.0, scale=stdnoise)
    else:
        noise = 0.0

    tseries = signal + noise
    return tseries


def ffa2(data):
    """ 
    Compute the FFA transform of a two-dimensional input

    Parameters
    ----------
    data : ndarray (2D)
        Input time series data in two-dimensional form with shape (m, p),
        where m is the number of signal periods and p the number of phase bins.

    Returns
    -------
    transform : ndarray (2D)
        The FFA transform of 'data', as a float32 2D array of shape (m, p)

    See Also
    --------
    ffafreq : trial frequencies in the output transform
    ffaprd : trial periods in the output transform
    """
    return libcpp.ffa2(data)


def ffa1(data, p):
    """ 
    Compute the FFA transform of a one-dimensional input (time series)
    at base period p

    Parameters
    ----------
    data : ndarray (1D)
        Input time series data. If N is the total number of samples in the
        data, the last N % p samples are ignored, as they do not form a
        complete pulse period
    p : int
        Base period of the transform, in number of samples

    Returns
    -------
    transform : ndarray (2D)
        The FFA transform of 'data', as a float32 2D array of shape (m, p), 
        where m is the number of complete pulse periods in the data

    See Also
    --------
    ffafreq : trial frequencies in the output transform
    ffaprd : trial periods in the output transform
    """
    if not data.ndim == 1:
        raise ValueError("input data must be one-dimensional")
    if not (isinstance(p, int) and p > 0):
        raise ValueError("p must be an integer > 1")
    if p > data.size:
        raise ValueError("p must be smaller than the total number of samples")
    m = data.size // p
    return ffa2(data[:m*p].reshape(m, p))


def ffafreq(N, p, dt=1.0):
    """
    Returns the trial frequencies that correspond to every folded profile
    in the FFA output.

    Parameters
    ----------
    N : int
        Total number of samples in the input data
    p : int
        Base period of the FFA transform in number of samples
    dt : float, optional
        Sampling time

    Returns
    -------
    freqs : ndarray
        Array with m elements containing the sequence of trial frequencies
        in the FFA output
    """
    if not (isinstance(N, int) and N > 0):
        raise ValueError("N must be a strictly positive integer")

    if not (isinstance(p, int) and p > 1):
        raise ValueError("p must be an integer > 1")

    if not N >= p:
        raise ValueError("p must be smaller than (or equal to) N")

    if not dt > 0:
        raise ValueError("dt must be strictly positive")

    f0 = 1.0 / p
    m = N // p
    if m == 1:
        f = np.asarray([f0])
    else:
        s = np.arange(m)
        f = (f0 - s / (m-1.0) * f0**2)
    f /= dt
    return f


def ffaprd(N, p, dt=1.0):
    """
    Returns the trial periods that correspond to every folded profile
    in the FFA output.

    Parameters
    ----------
    N : int
        Total number of samples in the input data
    p : int
        Base period of the FFA transform in number of samples
    dt : float, optional
        Sampling time

    Returns
    -------
    periods : ndarray
        Array with m elements containing the sequence of trial periods
    """
    return 1.0 / ffafreq(N, p, dt=dt)


def boxcar_snr(data, widths, stdnoise=1.0):
    """ 
    Compute the S/N ratio of pulse profile(s) for a range of 
    boxcar width trials.

    Parameters
    ----------
    data : ndarray
        Input profile(s). Can be of any shape, but the last axis
        must be pulse phase.
    widths : ndarray, 1D
        Trial pulse widths, expressed in number of phase bins.
    stdnoise : float
        Standard deviation of the background noise in all profiles.

    Returns
    -------
    snr : ndarray
        Output with the same shape as data, with an additional axis
        which represents trial pulse width index.
    """
    widths = np.asarray(widths, dtype=np.uint64)

    # Number of bins is the length of the last axis
    b = data.shape[-1]

    # Input to C++ function must be 2D
    cppinput = data.reshape(-1, b).astype(np.float32)
    m = cppinput.shape[0]
    snr = libcpp.snr2(cppinput, widths, stdnoise)
    shape = list(data.shape[:-1]) + [widths.size]
    return snr.reshape(shape)


def downsample(data, factor):
    """ Downsample an array by a real-valued factor.

    Parameters
    ----------
    data : array_like
        Time series data to downsample.
    factor : float
        Downsampling factor.

    Returns
    -------
    out : ndarray, float32
        Downsampled data.
    """
    return libcpp.downsample(data, factor)
