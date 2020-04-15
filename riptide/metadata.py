import os
import pprint
import json

from astropy.coordinates import SkyCoord
import astropy.units as uu
from schema import Schema, And, Or, Optional
from .reading import PrestoInf, SigprocHeader


SCHEMA_ITEMS = {
    Optional('source_name') : Or(str, None),
    Optional('skycoord') : Or(SkyCoord, None),
    Optional('dm'): Or(And(float, lambda x: x >= 0), None),
    Optional('mjd'): Or(And(float, lambda x: x >= 0), None),
    Optional('tobs'): Or(And(float, lambda x: x > 0), None),
    Optional('fname'): Or(str, None),

    # Accept any extra keys of type string with JSON-serializable values
    Optional(str): json.dumps
    }

SCHEMA = Schema(SCHEMA_ITEMS, ignore_extra_keys=True)


class Metadata(dict):
    """ 
    A dict subclass that carries information about an observation across all
    data products (TimeSeries, Periodogram, etc.)

    The 'attrs' dictionary can only have keys of type str and json-serializable
    values (there are some exceptions, see below). There are also reserved keys
    which, if present, must match the criteria below:

    - source_name: str
    - skycoord: astropy.coordinates.Skycoord
    - dm: float, positive
    - mjd: float, positive
    - tobs: float, strictly positive
    - fname: str

    If any of the above keys are NOT present, they will be set to
    None in the Metadata object.
    """
    def __init__(self, items={}):
        SCHEMA.validate(items)
        super(Metadata, self).__init__(items)
        
        for k in SCHEMA_ITEMS:
            if isinstance(k.schema, str):
                self.setdefault(k.schema, None)

    @classmethod
    def from_presto_inf(cls, inf):
        """ 
        Create Metadata object from PRESTO .inf file or PrestoInf object

        Parameters
        ----------
        inf : PrestoInf or str
            PrestoInf object or path to a PRESTO .inf file
        """
        # Interpret 'inf' as a file path if it is a string
        if type(inf) == str:
            inf = PrestoInf(inf)

        attrs = dict(inf)
        attrs['skycoord'] = inf.skycoord
        attrs['fname'] = os.path.realpath(inf.fname)
        attrs['tobs'] = attrs['tsamp'] * attrs['nsamp']
        return cls(attrs)

    @classmethod
    def from_sigproc(cls, sh, extra_keys={}):
        """ 
        Create Metadata object from SIGPROC dedispersed time series file,
        or SigprocHeader object.

        Parameters
        ----------
        sh : SigprocHeader or str
            SigprocHeader object or path to a PRESTO .inf file
        """
        # Interpret 'sh' as a file path if it is a string
        if type(sh) == str:
            sh = SigprocHeader(sh, extra_keys=extra_keys)

        if sh['nchans'] > 1:
            raise ValueError(f"File {sh.fname!r} contains multi-channel data (nchans = {sh['nchans']}), instead of a dedispersed time series")

        # Make sure this is a 32-bit dedispersed time series
        # We support either: 32-bit float data, or 8-bit data but only if signedness is specified in the header
        nbits = sh['nbits']
        if not nbits in {8, 32}:
            raise ValueError(f"Only 8-bit and 32-bit SIGPROC data are supported. File {sh.fname!r} contains {nbits}-bit data")
        if nbits == 8 and 'signed' not in sh:
            raise ValueError(f"SIGPROC Header says this is 8-bit data, but does not specify its signedness via the 'signed' key")

        attrs = dict(sh).copy()
        attrs['dm'] = attrs.get('refdm', None)
        attrs['skycoord'] = sh.skycoord
        attrs['source_name'] = attrs.get('source_name', None)
        attrs['mjd'] = attrs.get('tstart', None)
        attrs['fname'] = os.path.realpath(sh.fname)
        attrs['tobs'] = sh.tobs
        return cls(attrs)

    def to_dict(self):
        return dict(self)

    @classmethod
    def from_dict(cls, items):
        return cls(items)

    def __str__(self):
        return 'Metadata %s' % pprint.pformat(dict(self))

    def __repr__(self):
        return str(self)
