##### Standard imports #####
import pprint

##### Non-standard imports #####
import h5py
from astropy.coordinates import SkyCoord
import astropy.units as uu

##### Local module imports #####
from .reading import PrestoInf


class Metadata(dict):
    """ Carries information about an observation across all data products
    (TimeSeries, Periodogram, etc.)"""

    _required_attrs = {
        'source_name' : str,
        'skycoord' : SkyCoord,
        'dm': float,
        'mjd': float,
        }

    _HDF5_group_name = 'metadata'


    def __init__(self, attrs):
        """ Create new Metadata from a dictionary of attributes.

        The 'attrs' dictionary is expected to contain a certain number of required
        keys with a specific type, at the moment these are:
        source_name, skycoord, dm, mjd

        If any of the above keys are NOT present in 'attrs', they will be set to
        None in the Metadata object. If they are present, they must have the
        correct type or a ValueError will be raised.

        All the other input attributes must have types compatible with HDF5,
        otherwise errors may be raised down the line when saving to / loading
        from HDF5.
        """
        super().__init__(attrs)
        self._finalise()

    def _finalise(self):
        """ Check that required keys are either absent, in which case they will
        be set to None, or present with the correct type """
        req = self._required_attrs
        for key in req.keys():
            val = self.get(key, None)
            if val is None:
                self[key] = val
            elif not isinstance(val, req[key]):
                msg = 'Metadata key \'{k:s}\' must have type \'{t:s}\' instead of \'{ti:s}\''.format(
                    k=key,
                    t=req[key].__name__,
                    ti=type(val).__name__
                    )
                raise ValueError(msg)

    @classmethod
    def from_presto_inf(cls, inf):
        """ Create Metadata object from PRESTO .inf file or PrestoInf object

        Parameters:
        -----------
            inf: PrestoInf object, or str
                PrestoInf object, or path to a PRESTO .inf file.
        """
        # Interpret 'inf' as a file path if it is a string
        if type(inf) == str:
            inf = PrestoInf(inf)

        attrs = inf.parsed_attrs.copy()
        attrs['dm'] = inf.dm
        attrs['skycoord'] = SkyCoord(inf.raj, inf.decj, unit=(uu.hour, uu.deg), frame='icrs')
        attrs['source_name'] = inf.source
        attrs['mjd'] = inf.mjd
        return cls(attrs)

    @classmethod
    def from_sigproc(cls, sig, extra_keys={}):
        """ Create Metadata object from SIGPROC dedispersed time series file,
        or SigprocHeader object.
        """
        # Interpret 'sig' as a file path if it is a string
        if type(sig) == str:
            sig = SigprocHeader(sig, extra_keys=extra_keys)

        # Make sure this is a 32-bit dedispersed time series
        if sig['nbits'] != 32:
            raise ValueError('Only 32-bit data is currently supported. File \'{0:s}\' contains {1:d}-bit data.'.format(sig.fname, sig['nbits']))
        if sig['nchans'] > 1:
            raise ValueError('File \'{0:s}\' contains multi-channel data (nchans = {1:d}), instead of a dedispersed time series'.format(sig.fname, sig['nchans']))

        attrs = dict(sig).copy()
        attrs['dm'] = attrs.get('refdm', None)
        attrs['skycoord'] = sig.skycoord
        attrs['source_name'] = attrs.get('source_name', None)
        attrs['mjd'] = attrs.get('tstart', None)
        return cls(attrs)

    def __str__(self):
        return 'Metadata\n%s' % pprint.pformat(dict(self))

    def __repr__(self):
        return str(self)

    def _to_hdf5_attributes(self):
        """ Convert to a dictionary of attributes ready to be stored in an HDF5
        group. Special care must be taken with the SkyCoord object, that we
        convert into 'rajd' and 'decjd' keys. Any Metadata key whose value is
        None is not present in the output dictionary. """
        attrs = dict(self)
        skycoord = attrs.get('skycoord', None)

        # Replace 'skycoord' key by 'rajd' and 'decjd'
        if skycoord is not None:
            attrs['rajd'] = skycoord.ra.deg
            attrs['decjd'] = skycoord.dec.deg
            attrs.pop('skycoord', None)

        # Remove keys that have a value of None
        # We use a list and not a generator expression here, because the
        # generator expression would change as we iterate through it and
        # delete keys from the dictionary.
        keys_to_remove = [key for key in attrs.keys() if attrs[key] is None]
        for key in keys_to_remove:
            attrs.pop(key, None)
        return attrs

    @classmethod
    def _from_hdf5_attributes(cls, h5attrs):
        """ Private method used to load from HDF5. """
        rajd = h5attrs.get('rajd', None)
        decjd = h5attrs.get('decjd', None)
        if (rajd is None) or (decjd is None):
            skycoord = None
        else:
            skycoord = SkyCoord(rajd, decjd, unit=(uu.deg, uu.deg))

        # Build initialisation dict, remove rajd and decjd keys that where
        # solely used to store the SkyCoord object into HDF5
        attrs = {}
        attrs.update(h5attrs)
        attrs['skycoord'] = skycoord
        attrs.pop('rajd', None)
        attrs.pop('decjd', None)
        return cls(attrs)

    def _save_to_hdf5_file(self, h5file):
        """ Create a metadata group in given HDF5.File object, and write
        all metadata as attributes of that group. """
        mgroup = h5file.create_group(self._HDF5_group_name)
        mgroup.attrs.update(self._to_hdf5_attributes())

    @classmethod
    def _from_hdf5_file(cls, h5file):
        """ Create a Metadata object from the attributes stored in the metadata
        group of given HDF5.File. """
        return cls._from_hdf5_attributes(h5file[cls._HDF5_group_name].attrs)
