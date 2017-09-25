##### Non-standard imports #####
import numpy as np

##### Local imports #####
from .running_median import fast_running_median
from .libffa import downsample, generate_signal
from .reading import PrestoInf


class TimeSeries(object):
    """ Container for time series data to be searched with the FFA. """
    def __init__(self, data, tsamp, copy=False):
        """ Create a new TimeSeries object, a container passed to the FFA search code.

        Parameters:
        -----------
            data: array_like
                Time series to search.
            tsamp: float
                Sampling time of data.
            copy: bool
                If set to True, the resulting time series will hold a new copy of data,
                otherwise it only holds a reference to it. (default: False)
        """
        if copy:
            self.data = np.asarray(data, dtype=np.float32).copy()
        else:
            self.data = np.asarray(data, dtype=np.float32)
        self.tsamp = float(tsamp)

    def normalise(self, inplace=False):
        """ Normalize to zero mean and unit variance. if 'inplace' is False,
        a new TimeSeries object with the normalized data is returned. """
        m = self.data.mean()
        s = self.data.std()
        if inplace:
            self.data = (self.data - m) / s
        else:
            return TimeSeries((self.data - m) / s, self.tsamp)

    def deredden(self, width, minpts=101, inplace=False):
        """ Remove red noise. """
        width_samples = int(round(width / self.tsamp))
        rmed = fast_running_median(self.data, width_samples, minpts)
        if inplace:
            self.data -= rmed
        else:
            return TimeSeries(self.data - rmed, self.tsamp)

    def downsample(self, factor, inplace=False):
        """ Downsample by a real-valued factor. """
        if inplace:
            self.data = downsample(self.data, factor)
            self.tsamp *= factor
        else:
            return TimeSeries(downsample(self.data, factor), factor * self.tsamp)

    @classmethod
    def generate(cls, length, tsamp, period=1.0, phi0=0.0, ducy=0.02, amplitude=1.0, stdnoise=1.0):
        """ Generate a time series containing a periodic signal with a von Mises
        pulse profile.

        Parameters:
        -----------
            length: float
                Length of the data in seconds.
            tsamp: float
                Sampling time in seconds.
            period: float
                Period in seconds.
            phi0: float
                Initial pulse phase in number of periods.
            ducy: float
                Duty cycle of the pulse.
            amplitude: float
                Pulse amplitude.
            stdnoise: float
                Standard deviation of the background noise.
        """
        nsamp = int(length / tsamp + 0.5)
        period_samples = period / tsamp
        data = generate_signal(nsamp, period_samples, phi0=phi0, ducy=ducy, amplitude=amplitude, stdnoise=stdnoise)
        return cls(data, tsamp, copy=False)

    @classmethod
    def from_numpy_array(cls, array, tsamp, copy=False):
        return cls(data, tsamp, copy=copy)

    @classmethod
    def from_binary(cls, fname, tsamp, dtype=np.float32):
        data = np.fromfile(fname, dtype=dtype)
        return cls(data, tsamp, copy=False)

    @classmethod
    def from_npy_file(cls, fname, tsamp):
        data = np.load(fname)
        return cls(data, tsamp, copy=False)

    @classmethod
    def from_presto_inf(cls, fname):
        inf = PrestoInf(fname)
        return cls(inf.load_data(), tsamp=inf.tsamp)

    @classmethod
    def from_sigproc(cls, fname):
        raise NotImplementedError

    @property
    def nsamp(self):
        return self.data.size

    @property
    def length(self):
        return self.nsamp * self.tsamp

    @property
    def tobs(self):
        return self.length

    def __str__(self):
        name = type(self).__name__
        out = '{name} {{nsamp = {x.nsamp:d}, tsamp = {x.tsamp:.4e}, tobs = {x.length:.3f}}}'.format(name=name, x=self)
        return out

    def __repr__(self):
        return str(self)

    def save_hdf5(self, fname):
        # TODO: implement this
        raise NotImplementedError

    @classmethod
    def load_hdf5(cls, fname):
        # TODO: implement this
        raise NotImplementedError
