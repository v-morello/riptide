##### Non-standard imports #####
import numpy as np

##### Local module imports #####
from .libffa import downsample, ffa_transform_2d



def fold_array(data, psamp, nbins):
    """ Fold data into a two-dimensional array of single pulses.

    The input data are downsampled by a real-valued factor so that one signal
    period contains exactly the specified number of phase bins.

    Parameters:
    -----------
        data: ndarray or array-like
            Input time series data (one dimensional).
        psamp: float
            Folding period in number of samples.
        nbins: int
            Number of phase bins in every output single pulse. Cannot exceed
            psamp.

    Returns:
    --------
        single_pulses: ndarray
            Folded output. A two-dimensional array of shape
            (num_pulses, num_bins).
    """
    # Input sanity checks
    nbins = int(round(nbins))
    psamp = float(psamp)

    if nbins > psamp:
        raise ValueError("\'psamp\' ({0:.3f}) must be larger or equal to \'nbins\' ({1:d})".format(psamp, nbins))

    # Downsample input data and generate single pulse stack with requested
    # number of phase bins
    tim = downsample(
        np.asarray(data, dtype=np.float32),
        psamp / nbins
        )
    nperiods = tim.size // nbins
    single_pulses = tim[:nperiods*nbins].reshape(nperiods, nbins)
    return single_pulses


def downsample_lines(data, nlines):
    """ Downsample a two-dimensional array along its first (vertical) axis.

    Parameters:
    -----------
        data: ndarray
            Input data (two dimensional).
        nlines: int
            Requested number of lines in the downsampled output.

    Returns:
    --------
        output: ndarray
            Downsampled output.
    """
    if data.ndim != 2:
        raise ValueError('Input data must be a 2D array')

    nx, ny = data.shape

    nlines = int(round(nlines))
    nlines = min(nx, nlines)
    nlines = max(1, nlines)

    dsf = nx / float(nlines)

    # Create the output array in transposed order, this is more practical
    # to run the downsample() function
    output = np.zeros(shape=(ny, nlines), dtype=np.float32)

    # Transpose and put the Y axis in contiguous memory order
    data = data.T.copy()
    for ibin in range(ny):
        output[ibin] = downsample(data[ibin], dsf)
    output = output.T
    return output


class SubIntegrations(object):
    """ Container for folded data. """
    def __init__(self, subints, period, orig_nsamp, orig_tsamp):
        """ This should not be called directly. Use classmethods instead. """
        # This is enough information to compute all the properties of a
        # SubIntegrations object
        self._data = subints
        self._period = period
        self._orig_nsamp = orig_nsamp
        self._orig_tsamp = orig_tsamp

    @classmethod
    def _check_input(cls, data, tsamp, period, nbins, nsubs):
        nsamp = data.size
        length = nsamp * tsamp
        if not nbins > 1:
            raise ValueError('nbins must be strictly larger than 1')
        if nsubs is not None and not nsubs >= 1:
            raise ValueError('nsubs must be larger or equal to 1, or None')
        if period > length:
            raise ValueError('Folding period ({0:.3e}) exceeds the length of the data ({1:.3e})'.format(period, length))
        if period < tsamp * nbins:
            nbins_max = int(period / tsamp)
            raise ValueError('Cannot resolve one period with this many phase bins (requested = {0:d}, max. acceptable = {1:d})'.format(nbins, nbins_max))
        if nsubs is not None and nsubs * period > length:
            nsubs_max = int(length / period)
            raise ValueError('Requested number of sub-integrations is too high (requested = {0:d}, max. acceptable = {1:d})'.format(nsubs, nsubs_max))

    @property
    def data(self):
        return self._data

    @property
    def period(self):
        return self._period

    @property
    def orig_tsamp(self):
        """ Sampling time of the original time series data """
        return self._orig_tsamp

    @property
    def orig_nsamp(self):
        """ Number of samples in the original time series data, including
        incomplete signal periods. """
        return self._orig_nsamp

    @property
    def orig_length(self):
        """ Length of the original time series data, including incomplete
        signal periods. """
        return self.orig_nsamp * self.orig_tsamp

    @property
    def nsubs(self):
        return self._data.shape[0]

    @property
    def nbins(self):
        return self._data.shape[1]

    @property
    def npulses(self):
        """ Total number of complete pulses in the input data. """
        return int(self.orig_length / self.period)

    @property
    def tbin(self):
        """ Duration of a phase bin. """
        return self.period / self.nbins

    @property
    def tsub(self):
        """ Length of a sub-integration """
        return self.npulses * self.period / self.nsubs

    @property
    def normalised_profile(self):
        """ Summed profile, normalised to the same variance as the original
        time series data. """
        # Factor by which the original time series was downsampled
        # (by the fold_array() function)
        dsfactor = self.tbin/ self.orig_tsamp
        return self.data.sum(axis=0) * (self.npulses * dsfactor) ** -0.5

    def get_ffa_transform(self):
        return ffa_transform_2d(self.data)

    @classmethod
    def from_numpy_array(cls, data, tsamp, period, nbins=128, nsubs=None):
        """ Fold data into a two-dimensional sub-integrations array, wrapped by
        a SubIntegrations object.

        Parameters:
        -----------
            data: ndarray
                Input time series data to fold.
            tsamp: float
                Sampling time of the input data.
            period: float
                Signal period, same unit as tsamp.
            nbins: int
                Number of phase bins (default: 128)
            nsubs: int
                If set to None, one sub-integration per complete signal period
                will be created, and no downsampling along the time axis occurs.
                (default: None)
        """
        ### This may raise a number of ValueErrors if the input is not right
        cls._check_input(data, tsamp, period, nbins, nsubs)

        # NOTE: fold_array() ignores samples that correspond to the last
        # incomplete pulse in the input data. This means that 'pulses' may
        # have a slightly non-zero mean even if mean(tseries) = 0
        pulses = fold_array(data, period / tsamp, nbins)
        npulses = pulses.shape[0]

        if (nsubs is not None) and (nsubs < npulses):
            # The call to copy() is to put the data in contiguous memory order
            subints = downsample_lines(pulses, nsubs).copy()
        else:
            subints = pulses
        return cls(subints, period=period, orig_nsamp=data.size, orig_tsamp=tsamp)


    @classmethod
    def from_time_series(cls, tseries, period, nbins=128, nsubs=None):
        """ Fold TimeSeries into a two-dimensional sub-integrations array.

        Parameters:
        -----------
            tseries: TimeSeries
                Input time series data to fold
            period: float
                Signal period, same unit as tseries.tsamp
            nbins: int
                Number of phase bins (default: 128)
            nsubs: int
                If set to None, one sub-integration per signal period will be
                created, and no downsampling along the time axis occurs.
                (default: None)

        Returns:
        --------
            subints: SubIntegrations
                The folded and appropriately downsampled output, wrapped by
                a SubIntegrations object.
        """
        return cls.from_numpy_array(tseries.data, tseries.tsamp, period, nbins=nbins, nsubs=nsubs)
