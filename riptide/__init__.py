# NOTE: best to place this at the top in case we want to import
# it in other files
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

### Major classes
from .time_series import TimeSeries
from .periodogram import Periodogram
from .metadata import Metadata
from .candidate import Candidate

### Major functions
from .search import ffa_search

from .libffa import (
    ffa1,
    ffa2,
    ffafreq,
    ffaprd,
    generate_signal,
    downsample,
    boxcar_snr)

from .peak_detection import find_peaks

### Serialization
from .serialization import save_json, load_json

from .tests import test

__all__ = [
    'TimeSeries',
    'Periodogram',
    'Metadata',
    'Candidate',
    'ffa_search',
    'ffa1',
    'ffa2',
    'ffafreq',
    'ffaprd',
    'generate_signal',
    'downsample',
    'boxcar_snr',
    'find_peaks',
    'save_json',
    'load_json',
    'test'
    ]
