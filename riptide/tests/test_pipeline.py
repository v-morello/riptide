import os
import glob
import tempfile
from copy import deepcopy

from pytest import raises
import yaml
import numpy as np
from riptide import TimeSeries, load_json
from riptide.pipeline.pipeline import get_parser, run_program
from riptide.pipeline.config_validation import InvalidPipelineConfig, InvalidSearchRange
from .presto_generation import generate_data_presto

# NOTE 1:
# pipeline uses multiprocessing, to get proper coverage stats we need:
# * A .coveragerc file with the following options:
# [run]
# branch = True
# parallel = True
# concurrency = multiprocessing
# * Ensure that all instances of multiprocessing.Pool() have been closed and joined, as follows:
# >> pool.close()
# >> pool.join()

# NOTE 2:
# To print logging output in full, call pytest like this:
# py.test --capture=no -o log_cli=True <FILES>

# NOTE 3:
# To get coverage stats, run this in the base riptide directory:
# coverage run -m pytest && coverage combine && coverage report -m --omit riptide/_version.py


SIGNAL_PERIOD = 1.0
DATA_TOBS = 128.0
DATA_TSAMP = 256e-6


def runner_presto_fakepsr(fname_conf, outdir):
    # Write test data
    # NOTE: generate a signal bright enough to get harmonics and thus make sure
    # that the harmonic filter gets to run
    params = [
        # (dm, amplitude, ducy)
        (0.0 , 10.0, 0.05),
        (10.0, 20.0, 0.02),
        (20.0, 10.0, 0.05)
    ]

    for dm, amplitude, ducy in params:
        basename = f"fake_DM{dm:.3f}"
        generate_data_presto(
            outdir, basename, tobs=DATA_TOBS, tsamp=DATA_TSAMP, period=SIGNAL_PERIOD, 
            dm=dm, amplitude=amplitude, ducy=ducy
        )
        
    ### Run pipeline ###
    files = glob.glob(f'{outdir}/*.inf')
    cmdline_args = ['--config', fname_conf, '--outdir', outdir] + files
    parser = get_parser()
    args = parser.parse_args(cmdline_args)
    run_program(args)

    ### Check output sanity ###
    topcand_fname = f"{outdir}/candidate_0000.json"
    assert os.path.isfile(topcand_fname)

    topcand = load_json(topcand_fname)

    # NOTE: these checks depend on the RNG seed and the pipeline config
    assert abs(topcand.params['period'] - SIGNAL_PERIOD) < 1.00e-4
    assert topcand.params['dm'] == 10.0
    assert topcand.params['width'] == 13
    assert abs(topcand.params['snr'] - 18.5) < 0.15


def runner_presto_purenoise(fname_conf, outdir):
    """
    Check that pipeline runs well even if no candidates are found
    """
    dm = 0.0
    basename = f"purenoise_DM{dm:.3f}"
    generate_data_presto(
        outdir, basename, tobs=DATA_TOBS, tsamp=DATA_TSAMP, period=SIGNAL_PERIOD, 
        dm=dm, amplitude=0.0
    )

    ### Run pipeline ###
    files = glob.glob(f'{outdir}/*.inf')
    cmdline_args = ['--config', fname_conf, '--outdir', outdir] + files
    parser = get_parser()
    args = parser.parse_args(cmdline_args)
    run_program(args)

    ### Check output sanity ###
    assert not glob.glob(f"{outdir}/*.json")
    assert not glob.glob(f"{outdir}/*.png")


def load_yaml(fname):
    with open(fname, 'r') as fobj:
        return yaml.safe_load(fobj)


def save_yaml(items, fname):
    with open(fname, 'w') as fobj:
        return yaml.safe_dump(items, fobj)


def test_pipeline_presto_fakepsr():
    # NOTE: outdir is a full path (str)
    with tempfile.TemporaryDirectory() as outdir:
        fname_conf = os.path.join(os.path.dirname(__file__), 'pipeline_config_A.yml')
        runner_presto_fakepsr(fname_conf, outdir)

    with tempfile.TemporaryDirectory() as outdir:
        fname_conf = os.path.join(os.path.dirname(__file__), 'pipeline_config_B.yml')
        runner_presto_fakepsr(fname_conf, outdir)


def test_pipeline_presto_purenoise():
    with tempfile.TemporaryDirectory() as outdir:
        fname_conf = os.path.join(os.path.dirname(__file__), 'pipeline_config_A.yml')
        runner_presto_purenoise(fname_conf, outdir)

    with tempfile.TemporaryDirectory() as outdir:
        fname_conf = os.path.join(os.path.dirname(__file__), 'pipeline_config_B.yml')
        runner_presto_purenoise(fname_conf, outdir)


def test_config_validation():
    fname_conf = os.path.join(os.path.dirname(__file__), 'pipeline_config_A.yml')
    conf_correct = load_yaml(fname_conf)

    # Wrong parameter type
    with tempfile.TemporaryDirectory() as outdir:
        conf_wrong = deepcopy(conf_correct)
        conf_wrong['dmselect']['min'] = 'LOL'
        tmp = os.path.join(outdir, 'wrong_config.yaml')
        save_yaml(conf_wrong, tmp)
        with raises(InvalidPipelineConfig):
            runner_presto_fakepsr(tmp, outdir)

    # period_min too low        
    with tempfile.TemporaryDirectory() as outdir:
        conf_wrong = deepcopy(conf_correct)
        conf_wrong['ranges'][0]['ffa_search']['period_min'] = 1.0e-9
        tmp = os.path.join(outdir, 'wrong_config.yaml')
        save_yaml(conf_wrong, tmp)
        with raises(InvalidSearchRange):
            runner_presto_fakepsr(tmp, outdir)

    # too many phase bins requested to fold candidates     
    with tempfile.TemporaryDirectory() as outdir:
        conf_wrong = deepcopy(conf_correct)
        conf_wrong['ranges'][0]['candidates']['bins'] = int(42.0e9)
        tmp = os.path.join(outdir, 'wrong_config.yaml')
        save_yaml(conf_wrong, tmp)
        with raises(InvalidSearchRange):
            runner_presto_fakepsr(tmp, outdir)

    # non-contiguous search ranges        
    with tempfile.TemporaryDirectory() as outdir:
        conf_wrong = deepcopy(conf_correct)
        conf_wrong['ranges'][0]['ffa_search']['period_max'] = 0.50042
        tmp = os.path.join(outdir, 'wrong_config.yaml')
        save_yaml(conf_wrong, tmp)
        with raises(InvalidSearchRange):
            runner_presto_fakepsr(tmp, outdir)