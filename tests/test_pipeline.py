import os
import glob
import tempfile

import numpy as np
from riptide import TimeSeries, load_json
from riptide.pipeline.pipeline import get_parser, run_program

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

SIGNAL_PERIOD = 1.0
DATA_TOBS = 128.0
DATA_TSAMP = 256e-6


def generate_data_presto(basename, dm=0.0, amplitude=20.0, ducy=0.05):
    """
    """
    ### IMPORTANT: seed the RNG to get reproducible results ###
    np.random.seed(0)

    ts = TimeSeries.generate(
        DATA_TOBS, DATA_TSAMP, SIGNAL_PERIOD, 
        amplitude=amplitude, ducy=ducy, stdnoise=1.0
    )
    inf_text = INF_TEMPLATE.format(basename=basename, nsamp=ts.nsamp, tsamp=DATA_TSAMP, dm=dm)
    return inf_text, ts.data.astype(np.float32)


def runner_presto_fakepsr(outdir):
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
        inf_path = os.path.join(outdir, basename + '.inf')
        dat_path = os.path.join(outdir, basename + '.dat')
        inf_text, data = generate_data_presto(basename, dm=dm, amplitude=amplitude, ducy=ducy)

        with open(inf_path, 'w') as f:
            f.write(inf_text)

        data.tofile(dat_path)
        
    ### Run pipeline ###
    fname_conf = os.path.join(os.path.dirname(__file__), 'pipeline_config.yml')
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
    assert abs(topcand.params['snr'] - 18.6) < 0.1


def runner_presto_purenoise(outdir):
    """
    Check that pipeline runs well even if no candidates are found
    """
    dm = 0.0
    basename = f"purenoise_DM{dm:.3f}"
    inf_text, data = generate_data_presto(basename, dm=dm, amplitude=0.0)
    inf_path = os.path.join(outdir, basename + '.inf')
    dat_path = os.path.join(outdir, basename + '.dat')
    with open(inf_path, 'w') as f:
        f.write(inf_text)
    data.tofile(dat_path)

    ### Run pipeline ###
    fname_conf = os.path.join(os.path.dirname(__file__), 'pipeline_config.yml')
    files = glob.glob(f'{outdir}/*.inf')
    cmdline_args = ['--config', fname_conf, '--outdir', outdir] + files
    parser = get_parser()
    args = parser.parse_args(cmdline_args)
    run_program(args)


def test_pipeline_presto_fakepsr():
    # NOTE: outdir is a full path (str)
    with tempfile.TemporaryDirectory() as outdir:
        runner_presto_fakepsr(outdir)


def test_pipeline_presto_purenoise():
    with tempfile.TemporaryDirectory() as outdir:
        runner_presto_purenoise(outdir)
