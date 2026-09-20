from __future__ import annotations

import json
import os
import pprint
from collections.abc import Mapping

from astropy.coordinates import SkyCoord
from pydantic import BaseModel, ConfigDict, Field, model_validator

from .reading import PrestoInf, SigprocHeader


class _MetadataModel(BaseModel):
    """Validate reserved metadata fields and JSON-compatible extras."""

    model_config = ConfigDict(
        arbitrary_types_allowed=True,
        extra="allow",
    )

    source_name: str | None = None
    skycoord: SkyCoord | None = None
    dm: float | None = Field(default=None, ge=0)
    mjd: float | None = Field(default=None, ge=0)
    tobs: float | None = Field(default=None, gt=0)
    fname: str | None = None

    @model_validator(mode="before")
    @classmethod
    def validate_extra_values(cls, values):
        """Validate metadata mappings and arbitrary extra values."""
        if not isinstance(values, Mapping):
            raise ValueError("metadata must be a mapping")

        for key, value in values.items():
            if not isinstance(key, str):
                raise ValueError("metadata keys must be strings")
            if key not in cls.model_fields:
                try:
                    json.dumps(value)
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        f"metadata value for extra key {key!r} is not JSON serializable"
                    ) from exc
        return values


class Metadata(dict):
    """
    Carry observation metadata across all data products.

    Metadata is a dict subclass used by TimeSeries, Periodogram, and other
    data products.

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

    def __init__(self, items=None):
        if items is None:
            items = {}
        validated = _MetadataModel.model_validate(items)
        super().__init__(validated.model_dump(mode="python"))

    @classmethod
    def from_presto_inf(cls, inf):
        """
        Create Metadata object from PRESTO .inf file or PrestoInf object.

        Parameters
        ----------
        inf : riptide.reading.PrestoInf or str
            PrestoInf object or path to a PRESTO .inf file
        """
        # Interpret 'inf' as a file path if it is a string
        if type(inf) is str:
            inf = PrestoInf(inf)

        attrs = dict(inf)
        attrs["skycoord"] = inf.skycoord
        attrs["fname"] = os.path.realpath(inf.fname)
        attrs["tobs"] = attrs["tsamp"] * attrs["nsamp"]
        return cls(attrs)

    @classmethod
    def from_sigproc(cls, sh, extra_keys=None):
        """
        Create Metadata from a SIGPROC dedispersed time series file.

        The input may also be a SigprocHeader object.

        Parameters
        ----------
        sh : riptide.reading.SigprocHeader or str
            SigprocHeader object or path to a PRESTO .inf file
        """
        # Interpret 'sh' as a file path if it is a string
        if extra_keys is None:
            extra_keys = {}
        if type(sh) is str:
            sh = SigprocHeader(sh, extra_keys=extra_keys)

        if sh["nchans"] > 1:
            raise ValueError(
                f"File {sh.fname!r} contains multi-channel data "
                f"(nchans = {sh['nchans']}), instead of a dedispersed time series"
            )

        # Make sure this is a 32-bit dedispersed time series
        # We support either 32-bit float data or 8-bit data with signedness
        # specified in the header.
        nbits = sh["nbits"]
        if nbits not in {8, 32}:
            raise ValueError(
                "Only 8-bit and 32-bit SIGPROC data are supported. "
                f"File {sh.fname!r} contains {nbits}-bit data"
            )
        if nbits == 8 and "signed" not in sh:  # noqa: PLR2004
            raise ValueError(
                "SIGPROC Header says this is 8-bit data, but does not specify "
                "its signedness via the 'signed' key"
            )

        attrs = dict(sh).copy()
        attrs["dm"] = attrs.get("refdm", None)
        attrs["skycoord"] = sh.skycoord
        attrs["source_name"] = attrs.get("source_name", None)
        attrs["mjd"] = attrs.get("tstart", None)
        attrs["fname"] = os.path.realpath(sh.fname)
        attrs["tobs"] = sh.tobs
        return cls(attrs)

    def to_dict(self):
        """Return metadata as a dictionary."""
        return dict(self)

    @classmethod
    def from_dict(cls, items):
        """Create Metadata from a dictionary."""
        return cls(items)

    def __str__(self):
        """Return a human-readable representation."""
        return f"Metadata {pprint.pformat(dict(self))}"

    def __repr__(self):
        """Return the developer representation."""
        return str(self)
