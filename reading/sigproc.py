""" Read dedispersed time series from SIGPROC. """
##### Standard imports #####
import os
import struct

##### Non-standard imports #####
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as uu


# Standard SIGPROC keys and associated data types
sigproc_keys_database = {
    'source_name': str,
    'machine_id': int,
    'telescope_id': int,
    'src_raj': float,
    'src_dej': float,
    'az_start': float,
    'za_start': float,
    'data_type': int,
    'refdm': float,
    'fch1': float,
    'nchans': int,
    'nbits' : int,
    'tstart': float,
    'tsamp' : float,
    'nifs' : int,
    'barycentric': int,
    }
sigproc_header_start_flag = 'HEADER_START'
sigproc_header_end_flag = 'HEADER_END'

def read_str(fobj):
    """ Read string from open binary file object. """
    size, = struct.unpack('i', fobj.read(4))
    return fobj.read(size).decode(encoding='utf-8')

def read_attribute(fobj, keydb):
    """ Read SIGPROC {key, value} pair from open binary file object. """
    key = read_str(fobj)
    if key == sigproc_header_end_flag:
        return key, None

    if key not in keydb:
        raise KeyError(
            'Type of SIGPROC header attribute \'%s\' is unknown, please specify its type explicitly.' % key
            )
    else:
        atype = keydb[key]

    if atype == str:
        val = read_str(fobj)
    elif atype == int:
        val, = struct.unpack('i', fobj.read(4))
    elif atype == float:
        val, = struct.unpack('d', fobj.read(8))
    else:
        raise ValueError('Key \'%s\' has unsupported type \'%s\'' % (key, atype))
    return key, val

def read_all_attributes(fobj, keydb):
    """
    Returns:
    --------
        attrs: dict
            Dictionary of attributes.
        size: int
            Size of the header in bytes
    """
    attrs = {}
    while True:
        key, val = read_attribute(fobj, keydb)
        if key == sigproc_header_end_flag:
            break
        attrs[key] = val
    return attrs, fobj.tell()

def read_sigproc_header(fname, extra_attributes={}):
    """ Read SIGPROC header from a file.

    Parameters:
    -----------
        fname: str
            File name to read.
        extra_attributes: dict
            Optional {key: type} dictionary, specifying how to parse any
            non-standard keys that could be found in the header.

    Returns:
    --------
        header: dict
            Dictionary containing the SIGPROC header attributes.
        bytesize: int
            Size of the header in bytes.
    """
    keydb = sigproc_keys_database

    # Add any extra keys to header key database
    if extra_attributes:
        keydb = sigproc_keys_database.copy()
        keydb.update(extra_attributes)

    with open(fname, 'rb') as fobj:
        flag = read_str(fobj) # should be HEADER_START
        errmsg = 'File starts with \'{0:s}\' flag instead of the expected \'{1:s}\''.format(flag, sigproc_header_start_flag)
        assert flag == sigproc_header_start_flag, errmsg
        header, bytesize = read_all_attributes(fobj, keydb)
    return header, bytesize

def parse_float_coord(f):
    """ Parse coordinate in SIGPROC's own decimal floating point,
    to either hours (RA) or degrees (Dec).
    """
    sign = np.sign(f)
    x = abs(f)
    hh, x = divmod(x, 10000.)
    mm, ss = divmod(x, 100.)
    return sign * (hh + mm / 60.0 + ss / 3600.0)


class SigprocTimeSeries(object):
    """ Read and manipulate dedispersed time series written by SIGPROC. """
    def __init__(self, fname, extra_attributes={}):
        self._fname = os.path.abspath(fname)
        (self._header, self._header_bytesize) = read_sigproc_header(self._fname, extra_attributes)
        if self.header['nbits'] != 32:
            raise ValueError('Only 32-bit data is currently supported (this is {0:d}-bit data)'.format(self.header['nbits']))
        if self.header['nchans'] > 1:
            raise ValueError('This appears to be filterbank data (nchans = {0:d}), instead of a dedispersed time series'.format(self.header['nchans']))

    @property
    def fname(self):
        return self._fname

    @property
    def header_bytesize(self):
        return self._header_bytesize

    @property
    def header(self):
        return self._header

    @property
    def skycoord(self):
        rajd = parse_float_coord(self.header['src_raj'])
        dejd = parse_float_coord(self.header['src_dej'])
        return SkyCoord(rajd, dejd, unit=(uu.hour, uu.degree), frame='icrs')

    def load_data(self):
        with open(self.fname, 'rb') as fobj:
            fobj.seek(self.header_bytesize)
            data = np.fromfile(fobj, dtype=np.float32)
        return data
