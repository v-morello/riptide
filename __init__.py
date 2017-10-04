### Major classes
from .time_series import TimeSeries
from .periodogram import Periodogram
from .metadata import Metadata
from .folding import SubIntegrations

### Major functions
from .search import ffa_search

from .libffa import (
    ffa_transform_1d,
    ffa_transform_2d,
    generate_signal,
    downsample,
    get_snr
    )

from .peak_detection import find_peaks

__all__ = [
    'TimeSeries',
    'Periodogram',
    'Metadata',
    'SubIntegrations',
    'ffa_search',
    'ffa_transform_1d',
    'ffa_transform_2d',
    'generate_signal',
    'downsample',
    'get_snr',
    'find_peaks'
    ]
