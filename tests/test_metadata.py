from astropy.coordinates import SkyCoord
from pydantic import ValidationError
from pytest import raises

from riptide import Metadata

EXPECTED_DM = 10.0


def test_metadata_validates_reserved_fields_and_defaults():
    """Validate reserved fields while retaining dictionary behavior."""
    skycoord = SkyCoord("12h00m00s", "-30d00m00s")
    metadata = Metadata(
        {"dm": EXPECTED_DM, "skycoord": skycoord, "extra": [1, "two"]}
    )

    assert metadata["dm"] == EXPECTED_DM
    assert metadata["skycoord"] is skycoord
    assert metadata["extra"] == [1, "two"]
    assert metadata["source_name"] is None
    assert metadata["mjd"] is None
    assert metadata["tobs"] is None
    assert metadata["fname"] is None


def test_metadata_rejects_invalid_values():
    """Reject invalid reserved fields and non-serializable extras."""
    with raises(ValidationError):
        Metadata({"dm": -1.0})

    with raises(ValidationError):
        Metadata({"extra": object()})

    with raises(ValidationError):
        Metadata({1: "not a valid metadata key"})
