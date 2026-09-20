import glob
import os
import tempfile

from riptide import load_json
from riptide.pipeline.pipeline import get_parser, run_program

from .presto_generation import generate_data_presto

SIGNAL_PERIOD = 1.0
DATA_TOBS = 128.0
DATA_TSAMP = 256e-6
PERIOD_TOLERANCE = 1.00e-4
EXPECTED_DM = 10.0
EXPECTED_WIDTH = 13
SNR_TOLERANCE = 0.15


def runner_presto_fakepsr(fname_conf, outdir):
    """Run the pipeline against synthetic pulsar data."""
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

    files = glob.glob(f"{outdir}/*.inf")
    cmdline_args = ["--config", fname_conf, "--outdir", outdir] + files
    parser = get_parser()
    args = parser.parse_args(cmdline_args)
    run_program(args)

    topcand_fname = f"{outdir}/candidate_0000.json"
    assert os.path.isfile(topcand_fname)

    topcand = load_json(topcand_fname)

    # NOTE: these checks depend on the RNG seed and the pipeline config
    assert abs(topcand.params["period"] - SIGNAL_PERIOD) < PERIOD_TOLERANCE
    assert topcand.params["dm"] == EXPECTED_DM
    assert topcand.params["width"] == EXPECTED_WIDTH
    assert abs(topcand.params["snr"] - 18.5) < SNR_TOLERANCE


def runner_presto_purenoise(fname_conf, outdir):
    """Check that pipeline runs well even if no candidates are found."""
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

    files = glob.glob(f"{outdir}/*.inf")
    cmdline_args = ["--config", fname_conf, "--outdir", outdir] + files
    parser = get_parser()
    args = parser.parse_args(cmdline_args)
    run_program(args)

    assert not glob.glob(f"{outdir}/*.json")
    assert not glob.glob(f"{outdir}/*.png")


def test_pipeline_presto_fakepsr():
    """Test the pipeline with synthetic pulsar data."""
    # NOTE: outdir is a full path (str)
    with tempfile.TemporaryDirectory() as outdir:
        fname_conf = os.path.join(os.path.dirname(__file__), "pipeline_config_A.yml")
        runner_presto_fakepsr(fname_conf, outdir)

    with tempfile.TemporaryDirectory() as outdir:
        fname_conf = os.path.join(os.path.dirname(__file__), "pipeline_config_B.yml")
        runner_presto_fakepsr(fname_conf, outdir)


def test_pipeline_presto_purenoise():
    """Test the pipeline when no candidates are found."""
    with tempfile.TemporaryDirectory() as outdir:
        fname_conf = os.path.join(os.path.dirname(__file__), "pipeline_config_A.yml")
        runner_presto_purenoise(fname_conf, outdir)

    with tempfile.TemporaryDirectory() as outdir:
        fname_conf = os.path.join(os.path.dirname(__file__), "pipeline_config_B.yml")
        runner_presto_purenoise(fname_conf, outdir)
