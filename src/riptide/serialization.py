import base64
import importlib
import json

import astropy.units as uu
import numpy as np
import pandas
from astropy.coordinates import SkyCoord


# NOTE: Using importlib avoids placing "import riptide" at the top of the file,
# which solves circular import issues
def get_riptide_version():
    """Return the installed riptide version."""
    riptide = importlib.import_module("riptide")
    return riptide.__version__


# NOTE: it is implicitly assumed that any JSON-serializable riptide class is in
# the *base* riptide module
def get_class(clsname):
    """Return a class exported by the base riptide module."""
    riptide = importlib.import_module("riptide")
    return getattr(riptide, clsname)


class JSONEncoder(json.JSONEncoder):
    """Encode riptide and scientific Python objects as JSON."""

    def default(self, obj):  # noqa: PLR0911
        """Encode an object not handled by the base JSON encoder."""
        # NOTE: this method is called only for types not supported by default
        # Since Metadata is a dict (supported by default), it *never* gets called
        # for Metadata objects
        # TODO: find a way to handle Metadata correctly in the future

        # See: https://stackoverflow.com/questions/27909658/json-encoder-and-decoder-for-complex-numpy-arrays
        if isinstance(obj, np.ndarray):
            b64_bytes = base64.b64encode(np.ascontiguousarray(obj).data)
            b64_str = b64_bytes.decode()
            return {
                "__type__": "numpy.ndarray",
                "data": b64_str,
                "dtype": str(obj.dtype),
                "shape": obj.shape,
            }

        # Handle numpy numeric types
        if isinstance(obj, np.integer):
            return int(obj)

        if isinstance(obj, np.floating):
            return float(obj)

        if isinstance(obj, pandas.DataFrame):
            return {
                "__type__": "pandas.DataFrame",
                "values": self.default(obj.values),
                "columns": list(obj.columns),
            }

        if isinstance(obj, SkyCoord):
            return {
                "__type__": "astropy.SkyCoord",
                "rajd": obj.icrs.ra.deg,
                "decjd": obj.icrs.dec.deg,
                "frame": "icrs",
            }

        # If we reach that point, we assume that anything with a to_dict() method
        # is a riptide serializable object
        if hasattr(obj, "to_dict"):
            items = obj.to_dict()
            items["__type__"] = type(obj).__name__

            # Version the object properly
            if hasattr(obj, "version") and obj.version:
                items["__version__"] = obj.version
            else:
                items["__version__"] = get_riptide_version()
            return items

        return super().default(obj)


def object_hook(items):
    """Decode a serialized object represented as a dictionary."""
    if "__type__" not in items:
        return items

    typename = items["__type__"]

    if typename == "numpy.ndarray":
        b64_bytes = items["data"].encode()
        data = base64.b64decode(b64_bytes)
        return np.frombuffer(data, items["dtype"]).reshape(items["shape"])

    # NOTE: decoding is done from the deepest nodes of the tree first
    # Which means items['values'] is already a numpy array here
    if typename == "pandas.DataFrame":
        return pandas.DataFrame(items["values"], columns=items["columns"])

    if typename == "astropy.SkyCoord":
        ra = items["rajd"]
        dec = items["decjd"]
        frame = items["frame"]
        return SkyCoord(ra * uu.deg, dec * uu.deg, frame=frame)

    # If typename is not one of the above, then we assume it is
    # a riptide serializable object, and its dict-encoded form is
    # expected to have a __version__ key
    cls = get_class(typename)
    obj = cls.from_dict(items)

    # Handle version
    if "__version__" in items:
        version = items["__version__"]
    else:
        version = get_riptide_version()
    obj.version = version
    return obj


def from_json(s):
    """Decode a JSON string containing a riptide object or composition."""
    return json.loads(s, object_hook=object_hook)


def to_json(obj, **kwargs):
    """
    Encode an object as JSON, including riptide objects or compositions.

    Any keyword arguments are passed to json.dumps().
    """
    kwargs = dict(kwargs)
    kwargs["cls"] = JSONEncoder
    kwargs.setdefault("indent", 4)
    return json.dumps(obj, **kwargs)


def load_json(fname):
    """Load a JSON file containing a riptide object or composition."""
    with open(fname) as f:
        return from_json(f.read())


def save_json(fname, obj, **kwargs):
    """
    Save a riptide object or composition to a JSON file.

    Any keyword arguments are passed to json.dumps().
    """
    with open(fname, "w") as f:
        f.write(to_json(obj, **kwargs))
