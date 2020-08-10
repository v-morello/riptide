import os
import tempfile

import numpy as np
from pytest import raises
from riptide.apps.rseek import get_parser, run_program
from .presto_generation import generate_data_presto


SIGNAL_PERIOD = 1.0
SIGNAL_FREQ = 1.0 / SIGNAL_PERIOD
DATA_TOBS = 128.0
DATA_TSAMP = 256e-6

PARSER = get_parser()
EXPECTED_COLUMNS = {'period', 'freq', 'width', 'ducy', 'dm', 'snr'}
DEFAULT_OPTIONS = dict(Pmin=0.5, Pmax=2.0, bmin=480, bmax=520, smin=7.0, format='presto')


def dict2args(d):
    """
    Convert dictionary of options to command line argument list
    """
    args = []
    for k, v in d.items():
        args.append(f'--{k}')
        args.append(str(v))
    return args


def test_rseek_fakepsr():
    with tempfile.TemporaryDirectory() as outdir:
        generate_data_presto(
            outdir, 'data', tobs=DATA_TOBS, tsamp=DATA_TSAMP, period=SIGNAL_PERIOD,
            dm=0.0, amplitude=20.0, ducy=0.02
        )
        fname = os.path.join(outdir, 'data.inf')
        cmdline_args = dict2args(DEFAULT_OPTIONS) + [fname]
        args = PARSER.parse_args(cmdline_args)
        df = run_program(args)

    assert df is not None
    assert set(df.columns) == EXPECTED_COLUMNS

    # Results must be sorted in decreasing S/N order
    assert all(df.snr == df.sort_values('snr', ascending=False).snr)

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
            outdir, 'data', tobs=DATA_TOBS, tsamp=DATA_TSAMP, period=SIGNAL_PERIOD,
            dm=0.0, amplitude=0.0
        )
        fname = os.path.join(outdir, 'data.inf')
        cmdline_args = dict2args(DEFAULT_OPTIONS) + [fname]
        args = PARSER.parse_args(cmdline_args)
        df = run_program(args)
    
    assert df is None
