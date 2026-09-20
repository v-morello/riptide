"""Read dedispersed time series from SIGPROC."""

import os
import struct

import astropy.units as uu
import numpy as np
from astropy.coordinates import SkyCoord

# SIGPROC keys and associated data types
# Copied from Ewan Barr's sigpyproc
# NOTE:
# Values of type int are assumed to be stored as 32-bit in the header
# Values of type float are assumed to be stored as 64-bit (C double)
# Values of type bool are assumed to be stored as an unsigned char (8-bit)
sigproc_keydb = {
    "filename": str,
    "telescope_id": int,
    "telescope": str,
    "machine_id": int,
    "data_type": int,
    "rawdatafile": str,
    "source_name": str,
    "barycentric": int,
    "pulsarcentric": int,
    "az_start": float,
    "za_start": float,
    "src_raj": float,
    "src_dej": float,
    "tstart": float,
    "tsamp": float,
    "nbits": int,
    "nsamples": int,
    "fch1": float,
    "foff": float,
    "fchannel": float,
    "nchans": int,
    "nifs": int,
    "refdm": float,
    "flux": float,
    "period": float,
    "nbeams": int,
    "ibeam": int,
    "hdrlen": int,
    "pb": float,
    "ecc": float,
    "asini": float,
    "orig_hdrlen": int,
    "new_hdrlen": int,
    "sampsize": int,
    "bandwidth": float,
    "fbottom": float,
    "ftop": float,
    "obs_date": str,
    "obs_time": str,
    "accel": float,
    "signed": bool,
}


# These flags mark the boundaries of the header in a SIGPROC data file
HEADER_START = "HEADER_START"
HEADER_END = "HEADER_END"


def read_str(fobj):
    """Read string from open binary file object."""
    (size,) = struct.unpack("i", fobj.read(4))
    # NOTE: reading a string in a binary file via fobj.read() returns
    # a str in python2, but a bytes in python3. In python3 we need to called
    # bytes.decode() to get a str, but in python2 calling decode() gives us
    # a unicode object that must be cast to str again.
    # The code below works with both python2 and python3.
    return str(fobj.read(size).decode())


def read_attribute(fobj, keydb):
    """Read SIGPROC {key, value} pair from open binary file object."""
    key = read_str(fobj)
    if key == HEADER_END:
        return key, None

    atype = keydb.get(key, None)
    if atype is None:
        errmsg = (
            f"Type of SIGPROC header attribute '{key:s}' is unknown, please specify it"
        )
        raise KeyError(errmsg)

    if atype is str:
        val = read_str(fobj)
    elif atype is int:
        (val,) = struct.unpack("i", fobj.read(4))
    elif atype is float:
        (val,) = struct.unpack("d", fobj.read(8))
    elif atype is bool:
        (val,) = struct.unpack("B", fobj.read(1))  # B = unsigned char
        val = bool(val)
    else:
        errmsg = f"Key '{key:s}' has unsupported type '{atype:s}'"
        raise ValueError(errmsg)
    return key, val


def read_sigproc_header(fobj, extra_keys=None):
    """Read a SIGPROC header from an open file object.

    Parameters
    ----------
    fobj : file
        Open file object to read.
    extra_keys : dict
        Optional {key: type} dictionary, specifying how to parse any
        non-standard keys that could be found in the header

    Returns
    -------
    header : dict
        Dictionary containing the SIGPROC header attributes
    bytesize : int
        Size of the header in bytes
    """
    if extra_keys is None:
        extra_keys = {}
    keydb = sigproc_keydb

    # Add any extra keys to header key database
    if extra_keys:
        keydb = sigproc_keydb.copy()
        keydb.update(extra_keys)

    # Read HEADER_START flag
    fobj.seek(0)
    flag = read_str(fobj)
    errmsg = (
        f"File starts with '{flag:s}' flag instead of the expected '{HEADER_START:s}'"
    )
    assert flag == HEADER_START, errmsg

    # Read all header attributes
    attrs = {}
    while True:
        key, val = read_attribute(fobj, keydb)
        if key == HEADER_END:
            break
        attrs[key] = val

    return attrs, fobj.tell()


def parse_float_coord(f):
    """Parse a SIGPROC decimal coordinate.

    Convert it to either hours (RA) or degrees (Dec).
    """
    sign = np.sign(f)
    x = abs(f)
    hh, x = divmod(x, 10000.0)
    mm, ss = divmod(x, 100.0)
    return sign * (hh + mm / 60.0 + ss / 3600.0)


class SigprocHeader(dict):
    """Wrap information carried by a SIGPROC file header.

    This object behaves like a dictionary.
    """

    def __init__(self, fname, extra_keys=None):
        if extra_keys is None:
            extra_keys = {}
        self._fname = os.path.abspath(fname)
        with open(self.fname, "rb") as fobj:
            (attrs, self._bytesize) = read_sigproc_header(fobj, extra_keys)
        super().__init__(attrs)

    @property
    def fname(self):
        """Absolute path to original file."""
        return self._fname

    @property
    def bytesize(self):
        """Number of bytes occupied by the header in the original file."""
        return self._bytesize

    @property
    def bytes_per_sample(self):
        """Return the number of bytes in each sample."""
        return self["nchans"] * self["nbits"] // 8

    @property
    def nsamp(self):
        """Number of samples in the data."""
        return (os.path.getsize(self.fname) - self.bytesize) // self.bytes_per_sample

    @property
    def tobs(self):
        """Total length of the data in seconds."""
        return self.nsamp * self["tsamp"]

    @property
    def skycoord(self):
        """astropy.SkyCoord object with the coordinates of the source."""
        rajd = parse_float_coord(self["src_raj"])
        dejd = parse_float_coord(self["src_dej"])
        return SkyCoord(rajd, dejd, unit=(uu.hour, uu.degree), frame="icrs")
