##### Non-standard imports #####
import numpy as np
import h5py

##### Local imports #####
from .running_median import fast_running_median
from .libffa import downsample, generate_signal
from .reading import PrestoInf
from .metadata import Metadata


class TimeSeries(object):
    """ Container for time series data to be searched with the FFA. """
    def __init__(self, data, tsamp, metadata=None, copy=False):
        """ Create a new TimeSeries object, a container passed to the FFA search code.

        Parameters:
        -----------
            data: array_like
                Time series to search.
            tsamp: float
                Sampling time of data.
            metadata: Metadata object
                Optional Metadata object describing the observation from which
                the data originate. (default: None)
            copy: bool
                If set to True, the resulting time series will hold a new copy of data,
                otherwise it only holds a reference to it. (default: False)
        """
        if copy:
            self.data = np.asarray(data, dtype=np.float32).copy()
        else:
            self.data = np.asarray(data, dtype=np.float32)
        self.tsamp = float(tsamp)
        self.metadata = metadata if metadata is not None else Metadata({})

    def normalise(self, inplace=False):
        """ Normalize to zero mean and unit variance. if 'inplace' is False,
        a new TimeSeries object with the normalized data is returned. """
        m = self.data.mean()
        s = self.data.std()
        if inplace:
            self.data = (self.data - m) / s
        else:
            return TimeSeries((self.data - m) / s, self.tsamp, metadata=self.metadata)

    def deredden(self, width, minpts=101, inplace=False):
        """ Remove red noise. """
        width_samples = int(round(width / self.tsamp))
        rmed = fast_running_median(self.data, width_samples, minpts)
        if inplace:
            self.data -= rmed
        else:
            return TimeSeries(self.data - rmed, self.tsamp, metadata=self.metadata)

    def downsample(self, factor, inplace=False):
        """ Downsample by a real-valued factor. """
        if inplace:
            self.data = downsample(self.data, factor)
            self.tsamp *= factor
        else:
            return TimeSeries(downsample(self.data, factor), factor * self.tsamp, metadata=self.metadata)

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
        nsamp = int(round(length / tsamp))
        period_samples = period / tsamp
        data = generate_signal(nsamp, period_samples, phi0=phi0, ducy=ducy, amplitude=amplitude, stdnoise=stdnoise)
        metadata = Metadata({
            'source_name': 'fake',
            'signal_shape': 'Von Mises',
            'signal_period': period,
            'signal_initial_phase': phi0,
            'signal_duty_cycle': ducy,
            })
        return cls(data, tsamp, copy=False, metadata=metadata)

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
        metadata = Metadata.from_presto_inf(inf)
        return cls(inf.load_data(), tsamp=inf.tsamp, metadata=metadata)

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
        with h5py.File(fname, 'w') as fobj:
            # Create a group to store metadata, as attribute of said group
            self.metadata._save_to_hdf5_file(fobj)

            # Create a data group for the time series data itself
            data_group = fobj.create_group('data')
            data_group.attrs.update({
                'tsamp': self.tsamp
                })
            data_group.create_dataset('time_series', data=self.data, dtype=np.float32)

    @classmethod
    def load_hdf5(cls, fname):
        with h5py.File(fname, 'r') as fobj:
            metadata = Metadata._from_hdf5_file(fobj)
            data_group = fobj['data']
            data = data_group['time_series'].value
            tsamp = data_group.attrs['tsamp']
        return cls(data, tsamp, metadata=metadata, copy=False)
