import os
import tempfile

import numpy as np
from riptide import TimeSeries
from riptide.apps.rseek import get_parser, run_program


SIGNAL_PERIOD = 1.0
SIGNAL_FREQ = 1.0 / SIGNAL_PERIOD
DATA_TOBS = 128.0
DATA_TSAMP = 256e-6

PARSER = get_parser()
EXPECTED_COLUMNS = {"period", "freq", "width", "ducy", "dm", "snr"}
DEFAULT_OPTIONS = dict(
    Pmin=0.5, Pmax=2.0, bmin=480, bmax=520, smin=7.0, format="presto"
)

INF_TEMPLATE = """
 Data file name without suffix          =  {basename:s}
 Telescope used                         =  Parkes
 Instrument used                        =  Multibeam
 Object being observed                  =  Pulsar
 J2000 Right Ascension (hh:mm:ss.ssss)  =  00:00:01.0000
 J2000 Declination     (dd:mm:ss.ssss)  =  -00:00:01.0000
 Data observed by                       =  Kenji Oba
 Epoch of observation (MJD)             =  59000.000000
 Barycentered?           (1=yes, 0=no)  =  1
 Number of bins in the time series      =  {nsamp:d}
 Width of each time series bin (sec)    =  {tsamp:.12e}
 Any breaks in the data? (1=yes, 0=no)  =  0
 Type of observation (EM band)          =  Radio
 Beam diameter (arcsec)                 =  981
 Dispersion measure (cm-3 pc)           =  {dm:.12f}
 Central freq of low channel (Mhz)      =  1182.1953125
 Total bandwidth (Mhz)                  =  400
 Number of channels                     =  1024
 Channel bandwidth (Mhz)                =  0.390625
 Data analyzed by                       =  Space Sheriff Gavan
 Any additional notes:
    Input filterbank samples have 2 bits.
"""


def generate_data_presto(
    outdir,
    basename,
    tobs=128.0,
    tsamp=256e-6,
    period=1.0,
    dm=0.0,
    amplitude=20.0,
    ducy=0.05,
):
    """
    Generate some time series data with a fake signal, and save it in PRESTO
    inf/dat format in the specified output directory.

    Parameters
    ----------
    outdir : str
        Path to the output directory
    basename : str
        Base file name (not path) under which the .inf and .dat files
        will be saved.
    **kwargs: self-explanatory
    """
    ### IMPORTANT: seed the RNG to get reproducible results ###
    np.random.seed(0)

    ts = TimeSeries.generate(
        tobs, tsamp, period, amplitude=amplitude, ducy=ducy, stdnoise=1.0
    )
    inf_text = INF_TEMPLATE.format(
        basename=basename, nsamp=ts.nsamp, tsamp=tsamp, dm=dm
    )

    inf_path = os.path.join(outdir, f"{basename}.inf")
    dat_path = os.path.join(outdir, f"{basename}.dat")
    with open(inf_path, "w") as fobj:
        fobj.write(inf_text)
    ts.data.tofile(dat_path)


def dict2args(d):
    """
    Convert dictionary of options to command line argument list
    """
    args = []
    for k, v in d.items():
        args.append(f"--{k}")
        args.append(str(v))
    return args


def test_rseek_fakepsr():
    with tempfile.TemporaryDirectory() as outdir:
        generate_data_presto(
            outdir,
            "data",
            tobs=DATA_TOBS,
            tsamp=DATA_TSAMP,
            period=SIGNAL_PERIOD,
            dm=0.0,
            amplitude=20.0,
            ducy=0.02,
        )
        fname = os.path.join(outdir, "data.inf")
        cmdline_args = dict2args(DEFAULT_OPTIONS) + [fname]
        args = PARSER.parse_args(cmdline_args)
        df = run_program(args)

    assert df is not None
    assert set(df.columns) == EXPECTED_COLUMNS

    # Results must be sorted in decreasing S/N order
    assert all(df.snr == df.sort_values("snr", ascending=False).snr)

    # Check parameters of the top candidate
    # NOTE: these checks depend on the RNG seed and the program options
    topcand = df.iloc[0]
    assert abs(topcand.freq - SIGNAL_FREQ) < 0.1 / DATA_TOBS
    assert abs(topcand.snr - 18.5) < 0.15
    assert topcand.dm == 0
    assert topcand.width == 13


def test_rseek_purenoise():
    with tempfile.TemporaryDirectory() as outdir:
        generate_data_presto(
            outdir,
            "data",
            tobs=DATA_TOBS,
            tsamp=DATA_TSAMP,
            period=SIGNAL_PERIOD,
            dm=0.0,
            amplitude=0.0,
        )
        fname = os.path.join(outdir, "data.inf")
        cmdline_args = dict2args(DEFAULT_OPTIONS) + [fname]
        args = PARSER.parse_args(cmdline_args)
        df = run_program(args)

    assert df is None
