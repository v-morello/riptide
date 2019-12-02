### Standard library imports
import os
import ctypes

### Non-standard imports
import numpy as np
import numpy.ctypeslib as npct
from numpy import log, sin, cos, exp, pi

### Local imports
from .ffautils import generate_width_trials

###############################################################################

def load_libffa():
    # this file's dir/name
    fdir, __ = os.path.split(os.path.abspath(__file__))

    # load the library, using numpy mechanisms
    lib = npct.load_library(
        os.path.join(fdir, 'c_src', 'libffa'), '.'
        )

    ### Set argument and return types of the library's functions
    lib.py_transform.restype = None
    lib.py_transform.argtypes = [
        npct.ndpointer(dtype=np.float32, ndim=2, flags=('C_CONTIGUOUS', 'ALIGNED')),
        ctypes.c_size_t,
        ctypes.c_size_t,
        npct.ndpointer(dtype=np.float32, ndim=2, flags=('C_CONTIGUOUS', 'ALIGNED', 'WRITEABLE')),
        ]

    lib.get_snr_2d.restype = None
    lib.get_snr_2d.argtypes = [
        npct.ndpointer(dtype=np.float32, ndim=2, flags=('C_CONTIGUOUS', 'ALIGNED')),
        ctypes.c_size_t,
        ctypes.c_size_t,
        npct.ndpointer(dtype=np.int64, ndim=1, flags=('C_CONTIGUOUS', 'ALIGNED')),
        ctypes.c_size_t,
        ctypes.c_double,
        npct.ndpointer(dtype=np.float32, ndim=2, flags=('C_CONTIGUOUS', 'ALIGNED', 'WRITEABLE')),
        ]

    lib.downsample.restype = None
    lib.downsample.argtypes = [
        npct.ndpointer(dtype=np.float32, ndim=1, flags=('C_CONTIGUOUS', 'ALIGNED')), # input
        ctypes.c_size_t, # input size
        ctypes.c_double, # downsampling factor
        npct.ndpointer(dtype=np.float32, ndim=1, flags=('C_CONTIGUOUS', 'ALIGNED', 'WRITEABLE')) # output
        ]

    lib.py_periodogram.restype = None
    lib.py_periodogram.argtypes = [
        npct.ndpointer(dtype=np.float32, ndim=1, flags=('C_CONTIGUOUS', 'ALIGNED')), # input
        ctypes.c_size_t, # input size
        npct.ndpointer(dtype=float, ndim=1, flags=('C_CONTIGUOUS', 'ALIGNED')), # sequence: dsfactor
        npct.ndpointer(dtype=int, ndim=1, flags=('C_CONTIGUOUS', 'ALIGNED')), # sequence: bins_min
        npct.ndpointer(dtype=int, ndim=1, flags=('C_CONTIGUOUS', 'ALIGNED')), # sequence: bins_max
        ctypes.c_size_t, # number of ProcessingPlan steps
        npct.ndpointer(dtype=int, ndim=1, flags=('C_CONTIGUOUS', 'ALIGNED')), # sequence: widths
        ctypes.c_size_t, # number of width trials
        npct.ndpointer(dtype=np.float32, ndim=1, flags=('C_CONTIGUOUS', 'ALIGNED', 'WRITEABLE')), # period trials (output)
        npct.ndpointer(dtype=np.float32, ndim=1, flags=('C_CONTIGUOUS', 'ALIGNED', 'WRITEABLE')), # output S/N
        ]
    return lib


class LibFFA(object):
    """ A class that wraps all C functions as static class members """
    _lib = load_libffa()
    py_transform = _lib.py_transform
    get_snr_2d = _lib.get_snr_2d
    downsample = _lib.downsample
    py_periodogram = _lib.py_periodogram

###############################################################################

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
    tseries: ndarray (1D, float)
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
    if not data.ndim == 2:
        raise ValueError("input data must be two-dimensional")

    m, p = data.shape
    output = np.zeros(shape=(m,p), dtype=np.float32)
    LibFFA.py_transform(data.astype(np.float32), m, p, output)
    return output


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
    freqs: ndarray
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
    periods: ndarray
        Array with m elements containing the sequence of trial periods
    """
    return 1.0 / ffafreq(N, p, dt=dt)


def get_snr(data, stdnoise=1.0):
    """ 
    Mainly for test purposes. Compute the S/N ratio of pulse profile(s) for
    a range of boxcar width trials.

    Parameters
    ----------
    data: ndarray
        Input profile(s). Can be of any shape, as long as the last axis
        is pulse phase.
    stdnoise: float
        Standard deviation of the background noise in all profiles.

    Returns
    -------
    snr: ndarray
        Output with the same shape as data, except for the last axis
        which represents trial pulse width index.
    widths: ndarray, 1D
        Trial pulse widths, an array that has the same length as the last
        axis of 'snr'.
    """
    # Number of bins is the length of the last axis
    b = data.shape[-1]

    # Input to C function must be 2D
    cinput = data.reshape(-1, b).astype(np.float32)
    m = cinput.shape[0]

    # Width trials
    widths = generate_width_trials(b, ducy_max=0.5, wtsp=1.5)
    nw = widths.size

    # Prepare output
    out = np.zeros(shape=(m, nw), dtype=np.float32)
    LibFFA.get_snr_2d(cinput, m, b, widths, nw, stdnoise**2.0, out)

    # Reshape output properly
    shape = list(data.shape[:-1]) + [nw]
    out = out.reshape(shape)
    return out, widths


def downsample(data, factor):
    """ Downsample an array by a real-valued factor.

    Parameters
    ----------
    data: array_like
        Time series data to downsample.
    factor: float
        Downsampling factor.

    Returns
    -------
    out: ndarray, float32
        Downsampled data.
    """
    if not factor > 1:
        raise ValueError('factor must be > 1')
    outsize = int(len(data) / float(factor))
    out = np.zeros(outsize, dtype=np.float32)
    LibFFA.downsample(np.asarray(data, dtype=np.float32), len(data), factor, out)
    return out
