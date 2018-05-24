#
# PRESTO .inf file parsing module
# (c) Vincent Morello, 2018
#
import os
import json

import numpy
from astropy.coordinates import SkyCoord
import astropy.units as uu

# The original C code that reads/writes .inf files can be found here:
# https://github.com/scottransom/presto/blob/master/src/ioinf.c
# The official way to parse is based on line order.

def str2bool(str):
    """ Convert string to boolean. """
    return int(str) != 0

presto_inf_parsing_plan = [
    ('basename', str),
    ('telescope', str),
    ('instrument', str),
    ('source_name', str),
    ('raj', str),
    ('decj', str),
    ('observer', str),
    ('mjd', float),
    ('barycentered', str2bool),
    ('nsamp', int),
    ('tsamp', float),
    ('breaks', str2bool),
    ('obstype', str),
    ('fov', float),
    ('dm', float),
    ('fbot', float),
    ('bandwidth', float),
    ('nchan', int),
    ('cbw', float),
    ('analyst', str),
    ('notes', str)
    ]

def split_lines(lines):
    """ Split lines in a .inf file into description/value pairs. """
    sep = '= ' # NOTE: the trailing space is important
    pairs = []
    extra_lines = []
    for line in lines:
        try:
            descr, value = map(str.strip, line.split(sep))
            pairs.append((descr, value))
        except:
            extra_lines.append(line)
    return pairs, extra_lines

def parse_pairs(pairs, extra_lines):
    """ Parse the output of split_lines() into a dictionary. """
    onoff_pairs = pairs[12:-8]
    keyval_pairs = pairs[:12] + pairs[-8:]

    # "Additional notes" at the end of the file
    # We append that to the key-value pair list and parse it as any other
    notes = '\n'.join(extra_lines[1:]).strip()
    keyval_pairs.append(('notes', notes))

    # Parsed key-value pairs as dictionary
    items = {}
    for pair, plan_step in zip(keyval_pairs, presto_inf_parsing_plan):
        descr, value = pair
        keyname, keytype = plan_step
        items[keyname] = keytype(value)
    return items

def inf2dict(text):
    """ Parse the text of a PRESTO .inf file into a dictionary. """
    lines = text.strip().split('\n')
    pairs, extra_lines = split_lines(lines)
    return parse_pairs(pairs, extra_lines)

class PrestoInf(dict):
    """ Parse PRESTO's .inf files that contain dedispersed time series
    metadata. """
    def __init__(self, fname):
        self._fname = os.path.realpath(fname)
        with open(fname, 'r') as fobj:
            items = inf2dict(fobj.read())
        super(PrestoInf, self).__init__(items)

    @property
    def fname(self):
        """ Absolute path to original file. """
        return self._fname

    @property
    def data_fname(self):
        """ Path to the associated .dat file """
        # NOTE: second argument of rsplit() is 'maxsplit'
        return self.fname.rsplit('.', 1)[0] + '.dat'

    @property
    def skycoord(self):
        """ astropy.SkyCoord object with the coordinates of the source. """
        return SkyCoord(self['raj'], self['decj'], unit=(uu.hour, uu.degree))

    def load_data(self):
        """ Returns the associated time series data as a numpy float32 array. """
        return numpy.fromfile(self.data_fname, dtype=numpy.float32)
