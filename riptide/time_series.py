import copy

##### Non-standard imports #####
import numpy as np

##### Local imports #####
from .running_median import fast_running_median
from .libffa import downsample, generate_signal
from .reading import PrestoInf, SigprocHeader
from .metadata import Metadata
from .ffautils import autocov
from .folding import fold
from .timing import timing


class TimeSeries(object):
    """ Container for time series data to be searched with the FFA. """
    def __init__(self, data, tsamp, metadata=None, copy=False):
        """ Create a new TimeSeries object, a container passed to the FFA search
        code. In most cases, classmethods should be used instead.

        Parameters
        ----------
        data : array-like
            Time series data to search.
        tsamp : float
            Sampling time of data in seconds.
        metadata : Metadata or dict, optional
            Optional Metadata object / dict describing the observation from which
            the data originate
        copy : bool, optional
            If set to True, the resulting time series will hold a new copy of data,
            otherwise it only holds a reference to it
        """
        if copy:
            self._data = np.asarray(data, dtype=np.float32).copy()
        else:
            self._data = np.asarray(data, dtype=np.float32)
        self._tsamp = float(tsamp)
        self.metadata = Metadata(metadata) if metadata is not None else Metadata({})

        # Carrying a tobs attribute is quite practical in later stages of
        # the pipeline (peak detection in periodograms)
        self.metadata['tobs'] = self.length

    @property
    def data(self):
        """ numpy array holding the time series data, in float32 format. """
        return self._data

    @property
    def tsamp(self):
        """ Sampling time in seconds. """
        return self._tsamp

    def copy(self):
        """ Returns a new copy of the TimeSeries """
        return copy.deepcopy(self)

    def normalise(self, inplace=False):
        """ Normalise to zero mean and unit variance. if 'inplace' is False,
        a new TimeSeries object with the normalized data is returned.

        Parameters
        ----------
        inplace : bool, optional
            If set to True, perform the operation in-place, otherwise, a new
            TimeSeries object is returned

        Returns
        -------
        out : TimeSeries or None
            The normalised TimeSeries, if 'inplace' was set to False
        """
        # NOTE: use float64 accumulator to avoid saturation issues when the
        # data have large values
        m = self.data.mean(dtype=np.float64)
        v = self.data.var(dtype=np.float64)
        norm = v ** 0.5
        
        if inplace:
            self._data = (self.data - m) / norm
        else:
            return TimeSeries((self.data - m) / norm, self.tsamp, metadata=self.metadata)


    @timing
    def deredden(self, width, minpts=101, inplace=False):
        """ Subtract from the data an aproximate running median. To save time,
        this running median is computed on a downsampled copy of the data, then
        upsampled back to the original resolution and finally subtracted from
        the original data.

        Parameters
        ----------
        width : float
            Width of the running median window in seconds.
        minpts : int, optional
            Downsample the data so that the width of the running median window
            is equal to 'minpoints' samples. The running median will be computed
            on that downsampled copy of the data, and then upsampled
        inplace : bool, optional
            If set to True, perform the operation in-place, otherwise, a new
            TimeSeries object is returned

        Returns:
        --------
        out: TimeSeries or None
            The de-redended TimeSeries, if 'inplace' was set to False
        """
        width_samples = int(round(width / self.tsamp))
        rmed = fast_running_median(self.data, width_samples, minpts)
        if inplace:
            self._data -= rmed
        else:
            return TimeSeries(self.data - rmed, self.tsamp, metadata=self.metadata)

    def downsample(self, factor, inplace=False):
        """ Downsample data by a real-valued factor, by grouping and adding
        together consecutive samples (or fractions of samples).

        Parameters
        ----------
        factor : float
            Downsampling factor
        inplace : bool, optional
            If set to True, perform the operation in-place, otherwise, a new
            TimeSeries object is returned

        Returns
        -------
        out : TimeSeries or None
            The downsampled TimeSeries, if 'inplace' was set to False
        """
        if inplace:
            self._data = downsample(self.data, factor)
            self._tsamp *= factor
        else:
            return TimeSeries(downsample(self.data, factor), factor * self.tsamp, metadata=self.metadata)

    def fold(self, period, bins, subints=None):
        """
        Fold TimeSeries at given period.

        Parameters
        ----------
        period : float
            Period in seconds
        bins : int
            Number of phase bins
        subints: int or None, optional
            Number of desired sub-integrations. If None, the number of 
            sub-integrations will be the number of full periods that fit in 
            the data

        Returns
        -------
        folded : ndarray
            The folded data as a numpy array. Has two dimensions if subints
        """
        return fold(self, period, bins, subints=subints)

    @classmethod
    def generate(cls, length, tsamp, period, phi0=0.5, ducy=0.02, amplitude=10.0, stdnoise=1.0):
        """ Generate a time series containing a periodic signal with a von Mises
        pulse profile, and some background white noise (optional).

        Parameters
        ----------
        length : float
            Length of the data in seconds.
        tsamp : float
            Sampling time in seconds.
        phi0 : float, optional
            Initial pulse phase in number of periods
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
            Standard deviation of the background noise (default: 1.0).
            If set to 0, a noiseless signal is generated.

        Returns
        -------
        tseries : ndarray (1D, float)
            Output time series.
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
        """ Create a new TimeSeries from a numpy array (or array-like).

        Parameters
        ----------
        array : array-like
            The time series data.
        tsamp : float
            Sampling time of the data in seconds.
        copy : bool, optional
            If set to True, the resulting time series will hold a new copy of
            'array', otherwise it only holds a reference to it
        
        Returns
        -------
        out: TimeSeries
            TimeSeries object.
        """
        return cls(array, tsamp, copy=copy)

    @classmethod
    def from_binary(cls, fname, tsamp, dtype=np.float32):
        """ Create a new TimeSeries from a raw binary file, containing the
        time series data without any header or footer. This will work as long
        as the data can be loaded with numpy.fromfile().

        Parameters
        ----------
        fname : str
            File name to load.
        tsamp : float
            Sampling time of the data in seconds.
        dtype : numpy data type, optional
            Data type of the file

        Returns:
        --------
        out: TimeSeries
            TimeSeries object.
        """
        data = np.fromfile(fname, dtype=dtype)
        return cls(data, tsamp, copy=False)

    @classmethod
    def from_npy_file(cls, fname, tsamp):
        """ Create a new TimeSeries from a .npy file, written with numpy.save().

        Parameters
        ----------
        fname : str
            File name to load.
        tsamp : float
            Sampling time of the data in seconds.

        Returns
        -------
        out : TimeSeries
            TimeSeries object.
        """
        data = np.load(fname)
        return cls(data, tsamp, copy=False)

    @classmethod
    @timing
    def from_presto_inf(cls, fname):
        """ Create a new TimeSeries from a .inf file written by PRESTO. The
        associated .dat file must be in the same directory.

        Parameters
        ----------
        fname : str
            File name to load.

        Returns
        -------
        out: TimeSeries
            TimeSeries object.
        """
        inf = PrestoInf(fname)
        metadata = Metadata.from_presto_inf(inf)
        return cls(inf.load_data(), tsamp=inf['tsamp'], metadata=metadata)

    @classmethod
    @timing
    def from_sigproc(cls, fname, extra_keys={}):
        """ Create a new TimeSeries from a file written by SIGPROC's dedisperse
        routine.

        Parameters
        ----------
        fname : str
            File name to load.
        extra_keys : dict, optional
            Optional {key: type} dictionary. Use it to specify how to parse any
            non-standard keys that could be found in the header, or even to
            override the data type of standard keys.

            Example:
                {
                'telescope_diameter' : float,
                'num_trusses': int,
                'planet_name': str
                }

        Returns
        -------
        out : TimeSeries
            TimeSeries object.
        """
        sig = SigprocHeader(fname, extra_keys=extra_keys)

        # This call checks if the file contains a dedispersed time series
        # in either 8-bit or 32-bit format
        # For 8-bit data, the signedness is specified via the 'signed' boolean key
        metadata = Metadata.from_sigproc(sig, extra_keys=extra_keys)

        # Load time series data
        with open(fname, 'rb') as fobj:
            fobj.seek(sig.bytesize)
            if metadata['nbits'] == 8:
                dtype = np.int8 if metadata['signed'] else np.uint8
                # Don't forget to cast to float32 after reading !
                data = np.fromfile(fobj, dtype=dtype).astype(np.float32)
            else: # assume float32
                data = np.fromfile(fobj, dtype=np.float32)

        return cls(data, tsamp=sig['tsamp'], metadata=metadata)

    @property
    def nsamp(self):
        """ Number of samples in the data. """
        return self.data.size

    @property
    def length(self):
        """ Length of the data in seconds """
        return self.nsamp * self.tsamp

    @property
    def tobs(self):
        """ Length of the data in seconds. Alias of property 'length'. """
        return self.length

    def __str__(self):
        name = type(self).__name__
        out = '{name} {{nsamp = {x.nsamp:d}, tsamp = {x.tsamp:.4e}, tobs = {x.length:.3f}}}'.format(name=name, x=self)
        return out

    def __repr__(self):
        return str(self)

    @classmethod
    def from_dict(cls, items):
        return cls(items['data'], items['tsamp'], metadata=items['metadata'], copy=False)

    def to_dict(self):
        return {
            'data': self.data,
            'tsamp': self.tsamp,
            'metadata': self.metadata
        }

