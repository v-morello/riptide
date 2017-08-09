import numpy as np

### Local imports
from .running_median import fast_running_median
from .libffa import downsample

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
        """ Normalize to zero mean and unit variance. if 'inplace' is False, a new TimeSeries
        object with the normalized data is returned. """
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
    def from_binary(cls, fname, tsamp, dtype=np.float32):
        data = np.fromfile(fname, dtype=dtype)
        return cls(data, tsamp, copy=False)
    
    @classmethod
    def from_npy(cls, fname, tsamp):
        data = np.load(fname)
        return cls(data, tsamp, copy=False)
    
    @classmethod
    def from_presto_inf(cls, fname):
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
