### Major classes
from .time_series import TimeSeries
from .periodogram import Periodogram
from .metadata import Metadata

### Major functions
from .search import ffa_search

from .libffa import (
    ffa_transform_1d,
    ffa_transform_2d,
    generate_signal,
    downsample
    )

__all__ = [
    'TimeSeries',
    'Periodogram',
    'Metadata',
    'ffa_search',
    'ffa_transform_1d',
    'ffa_transform_2d',
    'generate_signal',
    'downsample'
    ]
