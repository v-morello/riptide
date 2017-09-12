##### Non-standard imports #####
import h5py
from astropy.coordinates import SkyCoord
import astropy.units as uu

##### Local module imports #####
from .reading import PrestoInf


class Metadata(object):
    """ Carries information about an observation across all data products
    (TimeSeries, Periodogram, etc.)"""

    # NOTE: Hacky workaround to be able to store python None values
    # as HDF5 attributes. Normally we should be able to do this with
    # an h5py.Empty object (the clean way) but that refuses to work.
    __NONE_VAL = '__None'

    def __init__(self, source=None, mjd=None, dm=None, skycoord=None):
        """ Create a new Metadata object, that carries information about a
        data source across all data products. All arguments are optional and
        set to 'None' by default.

        Parameters:
        -----------
            source: str
                Source name.
            dm: float
                Source dispersion measure in pc / cm^3.
            mjd: float
                MJD of the observation. Whether this represent the start
                or the midpoint of the observation is up to the user.
            skycoord: astropy.SkyCoord
                Coordinates of the source.
        """
        ### Source name
        if source is not None:
            self.source = str(source)
        else:
            self.source = None

        ### MJD
        if mjd is not None:
            self.mjd = float(mjd)
        else:
            self.mjd = None

        ### DM
        if dm is not None:
            self.dm = float(dm)
        else:
            self.dm = None

        ### Coordinates, stored as an astropy SkyCoord object
        if skycoord is not None:
            if type(skycoord) == SkyCoord:
                self.skycoord = skycoord
            else:
                raise ValueError(
                    '\'skycoord\' must be an astropy.SkyCoord object, instead of \'%s\'' % type(skycoord).__name__
                )
        else:
            self.skycoord = None

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

        skycoord = SkyCoord(inf.raj, inf.decj, unit=(uu.hour, uu.deg))
        return cls(
            source=inf.source,
            mjd=inf.mjd,
            dm=inf.dm,
            skycoord=skycoord
            )

    @classmethod
    def from_sigproc(cls):
        raise NotImplementedError

    def to_dict(self):
        if self.skycoord is None:
            rajd = None
            decjd = None
        else:
            rajd = self.skycoord.ra.deg
            decjd = self.skycoord.dec.deg

        return {
            'source': self.source,
            'mjd'   : self.mjd,
            'dm'    : self.dm,
            'rajd'  : rajd,
            'decjd' : decjd
            }

    def __str__(self):
        return 'Metadata %s' % self.to_dict()

    def __repr__(self):
        return str(self)

    @classmethod
    def _from_dict(cls, attrs):
        """ Private method used to load from HDF5. """
        rajd = attrs['rajd']
        decjd = attrs['decjd']
        if (rajd is None) or (decjd is None):
            skycoord = None
        else:
            skycoord = SkyCoord(rajd, decjd, unit=(uu.deg, uu.deg))
        return cls(
            source=attrs['source'],
            mjd=attrs['mjd'],
            dm=attrs['dm'],
            skycoord=skycoord
            )

    def save_hdf5(self, fname):
        """ Save Metadata as a standalone HDF5 file. Mostly for test
        purposes. """
        with h5py.File(fname, 'w') as fobj:
            metadata_group = fobj.create_group('metadata')
            # Convert to dict then replace 'None' values by
            # __NONE_VAL
            attrs = self.to_dict()
            for key, val in attrs.items():
                if val is None:
                    attrs[key] = self.__NONE_VAL
            metadata_group.attrs.update(attrs)

    @classmethod
    def load_hdf5(cls, fname):
        """ Load Metadata from an HDF5 file. Mostly for test purposes. """
        with h5py.File(fname, 'r') as fobj:
            metadata_group = fobj['metadata']
            # The line below must be put within the 'with' statement, otherwise
            # 'metadata_group' becomes inaccessible
            attrs = dict(metadata_group.attrs.items())

            # Replace __NONE_VAL by python None
            for key, val in attrs.items():
                if val == cls.__NONE_VAL:
                    attrs[key] = None
        return cls.from_dict(attrs)
