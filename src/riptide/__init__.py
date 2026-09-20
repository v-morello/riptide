# NOTE: best to place this at the top in case we want to import
# it in other files
from ._version import version as __version__
from .candidate import Candidate
from .libffa import boxcar_snr, downsample, ffa1, ffa2, ffafreq, ffaprd, generate_signal
from .metadata import Metadata
from .peak_detection import find_peaks
from .periodogram import Periodogram
from .running_medians import fast_running_median, running_median
from .search import ffa_search
from .serialization import load_json, save_json
from .time_series import TimeSeries

__all__ = [
    "Candidate",
    "Metadata",
    "Periodogram",
    "TimeSeries",
    "__version__",
    "boxcar_snr",
    "downsample",
    "fast_running_median",
    "ffa1",
    "ffa2",
    "ffa_search",
    "ffafreq",
    "ffaprd",
    "find_peaks",
    "generate_signal",
    "load_json",
    "running_median",
    "save_json",
]
