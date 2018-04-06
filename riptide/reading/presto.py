"""
Module to parse PRESTO time series stored as .inf/.dat file pairs.
"""
### Standard library imports
import os
import glob
import numpy

### Non-standard imports
from astropy.coordinates import SkyCoord
import astropy.units as uu


#####################################################

# NOTE: all of these 'incorporate' function must have the same
# signature (PrestoInf, key, value, value_type) even if some
# of the arguments end up being useless.

def _incorporate(inf, key, val, vtype):
    setattr(inf, key, vtype(val))
    inf.parsed_keys.add(key)

def _incorporate_bool(inf, key, val, vtype):
    setattr(inf, key, bool(int(val)))
    inf.parsed_keys.add(key)

def _incorporate_onoff(inf, key, val, vtype):
    # Create empty list attribute to hold on/off bin pairs
    # if it does not exist already
    if not hasattr(inf, key):
        setattr(inf, key, [])
        inf.parsed_keys.add(key)

    start, end = map(int, map(str.strip, val.split(',')))
    onoff = getattr(inf, key)
    onoff.append((start, end))

#####################################################

class PrestoInf(object):
    """ """
    ### Parsing plan ###
    # here the keys are substrings identifying a line in a .inf file
    # the values are tuples (attribute_name, action, attribute_type)
    # used to decide how to incorporate the parsed data into the
    # final object
    _parsing_plan = {
        'data file name': ('basename', _incorporate, str),
        'telescope': ('telescope', _incorporate, str),
        'instrument': ('instrument', _incorporate, str),
        'object being observed': ('source', _incorporate, str),
        'j2000 right ascension': ('raj', _incorporate, str),
        'j2000 declination': ('decj', _incorporate, str),
        'observed by': ('observer', _incorporate, str),
        'epoch of observation': ('mjd', _incorporate, float),
        'barycentered': ('barycentered', _incorporate_bool, bool),
        'number of bins': ('nsamp', _incorporate, int),
        'width of each': ('tsamp', _incorporate, float),
        'any breaks': ('breaks', _incorporate_bool, bool),
        'type of observation': ('obstype', _incorporate, str),
        'beam diameter': ('bdiam', _incorporate, float),
        'dispersion measure': ('dm', _incorporate, float),
        'central freq of low channel': ('fbot', _incorporate, float),
        'total bandwidth': ('bandwidth', _incorporate, float),
        'number of channels': ('nchan', _incorporate, int),
        'channel bandwidth': ('cbw', _incorporate, float),
        'data analyzed by': ('analyst', _incorporate, str),

        # Special treatment for on/off bin pair entries
        # Here, in tuples (start, end) are appended to the specified attribute,
        # which is a list. The third element of the tuple is ignored here.
        'on/off bin pair': ('onoff', _incorporate_onoff, None),
        }

    def __init__(self, fname):
        """ Load PRESTO time series data into a convenient object. 'fname' can be either
        without or with the suffix .inf
        """
        if not fname.lower().endswith('.inf'):
            fname = fname + '.inf'

        self.fname = os.path.abspath(fname)
        self.parsed_keys = set()

        # Read .inf file into a dictionary {description (str): value (str)}
        dvdict = self._read(fname)

        # Incorporate values as attributes, following parsing plan
        self._incorporate_dvdict(dvdict)

    @classmethod
    def _read(cls, fname, delimiter='='):
        """ Read .inf file into dictionary {description: str_value}"""
        dvdict = {}
        with open(fname, 'r') as fobj:
            for line in fobj.readlines():
                try:
                    chunks = list(map(str.strip, line.split(delimiter)))
                    descr = chunks[0]
                    value = chunks[-1]
                    dvdict[descr] = value
                except Exception as err:
                    pass
        return dvdict

    @classmethod
    def _find_plan(cls, descr):
        return next(
            plan
            for partial_descr, plan in cls._parsing_plan.items()
            if partial_descr.lower() in descr.lower()
            )

    def _incorporate_dvdict(self, dvdict):
        cls = type(self)
        for descr, val in dvdict.items():
            try:
                attr_name, func, attr_type = cls._find_plan(descr)
            except StopIteration:
                #print('No parsing action matching \'%s\' was found' % descr)
                continue
            func(self, attr_name, val, attr_type)

    def load_data(self):
        """ Returns the associated time series data as a numpy float32 array. """
        return numpy.fromfile(self.data_fname, dtype=numpy.float32)

    @property
    def parsed_attrs(self):
        return {
            key: getattr(self, key)
            for key in self.parsed_keys
            }

    @property
    def skycoord(self):
        return SkyCoord(self.raj, self.decj, unit=(uu.hour, uu.deg), frame='icrs')

    @property
    def data_fname(self):
	# NOTE: second argument of rsplit() is 'maxsplit'
        return self.fname.rsplit('.', 1)[0] + '.dat'

    def __str__(self):
        cls = type(self).__name__
        return '%s %s' % (cls, self.basename)

    def __repr__(self):
        return str(self)
