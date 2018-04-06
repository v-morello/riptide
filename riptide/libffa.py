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
    fdir, fname = os.path.split(os.path.abspath(__file__))

    # load the library, using numpy mechanisms
    lib = npct.load_library(
        os.path.join(fdir, '..', 'c_src', 'libffa'), '.'
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
        ctypes.c_size_t,
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
        ctypes.c_size_t, # number of threads for S/N evaluation
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

def generate_signal(nsamp, period, phi0=0.0, ducy=0.01, amplitude=1.0, stdnoise=1.0):
    """ Generate a time series containing a periodic signal with a von Mises
    pulse profile. This function is mostly for test purposes.

    Parameters:
    -----------
        nsamp: int
            Number of samples to generate.
        period: float
            Period in number of samples.
        phi0: float
            Initial pulse phase in number of periods.
        ducy: float
            Duty cycle of the pulse.
        amplitude: float
            Pulse amplitude.
        stdnoise: float
            Standard deviation of the background noise.

    Returns:
    --------
        tseries: ndarray (1D, float)
           Output time series.
    """
    # von mises parameter
    kappa = log(2.0) / (2.0 * sin(pi*ducy/2.0)**2)

    # Generate pulse train
    phase_radians = (np.arange(nsamp, dtype=float) / period - phi0) * (2 * pi)
    signal = amplitude * exp(kappa*(cos(phase_radians) - 1.0))

    # Add noise
    if stdnoise > 0.0:
        noise = np.random.normal(size=nsamp, loc=0.0, scale=stdnoise)
    else:
        noise = 0.0

    tseries = signal + noise
    return tseries


def ffa_transform_2d(data):
    """ Compute the FFA transform of an array of profiles.

    Parameters:
    -----------
        data: ndarray (2D)
            Input time series data in two-dimensional shape (m, b), where m
            is the number of signal periods and b the number of phase bins.

    Returns:
    --------
        transform: ndarray (2D)
            The FFA transform of 'data', as a 2D array of shape (m, b).
            Line number k corresponds to a summation path with a top-to-bottom
            right shift of k bins.
    """
    m, b = data.shape
    output = np.zeros(shape=(m,b), dtype=np.float32)
    LibFFA.py_transform(data.astype(np.float32), m, b, output)
    return output


def ffa_transform_1d(data, pnum):
    """ Compute the FFA transform of a time series at a given integer
    period number. The algorithm only considers a number of time
    samples that is a multiple of 'pnum', and therefore ignores the
    last nsamp % pnum samples.

    Parameters:
    -----------
        data: ndarray (1D)
            Input time series data. Trailing samples with indices larger
            than the largest multiple of pnum are ignored.
        pnum: int
            Signal period in samples.

    Returns:
    --------
        transform: ndarray (2D)
            The FFA transform of 'data', as a 2D array of shape
            (m, pnum), where m is the number of complete periods of the signal.
            Line number k corresponds to a summation path with a top-to-bottom
            right shift of k bins.
    """
    m = data.size // pnum
    return ffa_transform_2d(data[:m*pnum].reshape(m, pnum))


def get_snr(data, stdnoise=1.0, threads=1):
    """ Compute the S/N ratio of pulse profile(s) for a range of boxcar width
    trials.

    Parameters:
    -----------
        data: ndarray
            Input profile(s). Can be of any shape, as long as the last axis
            is pulse phase.
        stdnoise: float
            Standard deviation of the background noise in all profiles.
        threads: int
            Number of OpenMP threads to use. Each threads processes a block of
            profiles.

    Returns:
    --------
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
    LibFFA.get_snr_2d(cinput, m, b, widths, nw, stdnoise**2.0, threads, out)

    # Reshape output properly
    shape = list(data.shape[:-1]) + [nw]
    out = out.reshape(shape)
    return out, widths


def downsample(data, factor):
    """ Downsample an array by a real-valued factor.

    Parameters:
    -----------
        data: array_like
            Time series data to downsample.
        factor: float
            Downsampling factor.

    Returns:
    --------
        out: ndarray, float32
            Downsampled data.
    """
    if not factor > 1:
        raise ValueError('factor must be > 1')
    outsize = int(len(data) / float(factor))
    out = np.zeros(outsize, dtype=np.float32)
    LibFFA.downsample(np.asarray(data, dtype=np.float32), len(data), factor, out)
    return out
