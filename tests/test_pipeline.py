import os
import glob
import tempfile
from copy import deepcopy

import yaml
import numpy as np
from pytest import raises
from riptide import load_json
from riptide import TimeSeries
from riptide.pipeline.pipeline import get_parser, run_program
from riptide.pipeline.config_validation import InvalidPipelineConfig, InvalidSearchRange

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
# pytest --capture=no -o log_cli=True <FILES>

# NOTE 3:
# To get coverage stats, run this in the base riptide directory:
# coverage run -m pytest && coverage combine && coverage report -m --omit src/riptide/_version.py


SIGNAL_PERIOD = 1.0
DATA_TOBS = 128.0
DATA_TSAMP = 256e-6


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


def runner_presto_fakepsr(fname_conf, outdir):
    # Write test data
    # NOTE: generate a signal bright enough to get harmonics and thus make sure
    # that the harmonic filter gets to run
    params = [
        # (dm, amplitude, ducy)
        (0.0, 10.0, 0.05),
        (10.0, 20.0, 0.02),
        (20.0, 10.0, 0.05),
    ]

    for dm, amplitude, ducy in params:
        basename = f"fake_DM{dm:.3f}"
        generate_data_presto(
            outdir,
            basename,
            tobs=DATA_TOBS,
            tsamp=DATA_TSAMP,
            period=SIGNAL_PERIOD,
            dm=dm,
            amplitude=amplitude,
            ducy=ducy,
        )

    ### Run pipeline ###
    files = glob.glob(f"{outdir}/*.inf")
    cmdline_args = ["--config", fname_conf, "--outdir", outdir] + files
    parser = get_parser()
    args = parser.parse_args(cmdline_args)
    run_program(args)

    ### Check output sanity ###
    topcand_fname = f"{outdir}/candidate_0000.json"
    assert os.path.isfile(topcand_fname)

    topcand = load_json(topcand_fname)

    # NOTE: these checks depend on the RNG seed and the pipeline config
    assert abs(topcand.params["period"] - SIGNAL_PERIOD) < 1.00e-4
    assert topcand.params["dm"] == 10.0
    assert topcand.params["width"] == 13
    assert abs(topcand.params["snr"] - 18.5) < 0.15


def runner_presto_purenoise(fname_conf, outdir):
    """
    Check that pipeline runs well even if no candidates are found
    """
    dm = 0.0
    basename = f"purenoise_DM{dm:.3f}"
    generate_data_presto(
        outdir,
        basename,
        tobs=DATA_TOBS,
        tsamp=DATA_TSAMP,
        period=SIGNAL_PERIOD,
        dm=dm,
        amplitude=0.0,
    )

    ### Run pipeline ###
    files = glob.glob(f"{outdir}/*.inf")
    cmdline_args = ["--config", fname_conf, "--outdir", outdir] + files
    parser = get_parser()
    args = parser.parse_args(cmdline_args)
    run_program(args)

    ### Check output sanity ###
    assert not glob.glob(f"{outdir}/*.json")
    assert not glob.glob(f"{outdir}/*.png")


def load_yaml(fname):
    with open(fname, "r") as fobj:
        return yaml.safe_load(fobj)


def save_yaml(items, fname):
    with open(fname, "w") as fobj:
        return yaml.safe_dump(items, fobj)


def test_pipeline_presto_fakepsr():
    # NOTE: outdir is a full path (str)
    with tempfile.TemporaryDirectory() as outdir:
        fname_conf = os.path.join(os.path.dirname(__file__), "pipeline_config_A.yml")
        runner_presto_fakepsr(fname_conf, outdir)

    with tempfile.TemporaryDirectory() as outdir:
        fname_conf = os.path.join(os.path.dirname(__file__), "pipeline_config_B.yml")
        runner_presto_fakepsr(fname_conf, outdir)


def test_pipeline_presto_purenoise():
    with tempfile.TemporaryDirectory() as outdir:
        fname_conf = os.path.join(os.path.dirname(__file__), "pipeline_config_A.yml")
        runner_presto_purenoise(fname_conf, outdir)

    with tempfile.TemporaryDirectory() as outdir:
        fname_conf = os.path.join(os.path.dirname(__file__), "pipeline_config_B.yml")
        runner_presto_purenoise(fname_conf, outdir)


def test_config_validation():
    fname_conf = os.path.join(os.path.dirname(__file__), "pipeline_config_A.yml")
    conf_correct = load_yaml(fname_conf)

    # Wrong parameter type
    with tempfile.TemporaryDirectory() as outdir:
        conf_wrong = deepcopy(conf_correct)
        conf_wrong["dmselect"]["min"] = "LOL"
        tmp = os.path.join(outdir, "wrong_config.yaml")
        save_yaml(conf_wrong, tmp)
        with raises(InvalidPipelineConfig):
            runner_presto_fakepsr(tmp, outdir)

    # period_min too low
    with tempfile.TemporaryDirectory() as outdir:
        conf_wrong = deepcopy(conf_correct)
        conf_wrong["ranges"][0]["ffa_search"]["period_min"] = 1.0e-9
        tmp = os.path.join(outdir, "wrong_config.yaml")
        save_yaml(conf_wrong, tmp)
        with raises(InvalidSearchRange):
            runner_presto_fakepsr(tmp, outdir)

    # too many phase bins requested to fold candidates
    with tempfile.TemporaryDirectory() as outdir:
        conf_wrong = deepcopy(conf_correct)
        conf_wrong["ranges"][0]["candidates"]["bins"] = int(42.0e9)
        tmp = os.path.join(outdir, "wrong_config.yaml")
        save_yaml(conf_wrong, tmp)
        with raises(InvalidSearchRange):
            runner_presto_fakepsr(tmp, outdir)

    # non-contiguous search ranges
    with tempfile.TemporaryDirectory() as outdir:
        conf_wrong = deepcopy(conf_correct)
        conf_wrong["ranges"][0]["ffa_search"]["period_max"] = 0.50042
        tmp = os.path.join(outdir, "wrong_config.yaml")
        save_yaml(conf_wrong, tmp)
        with raises(InvalidSearchRange):
            runner_presto_fakepsr(tmp, outdir)
