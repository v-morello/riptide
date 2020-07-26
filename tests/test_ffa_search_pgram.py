# Test ffa_search() and Periodogram class
import tempfile

import numpy as np
import matplotlib.pyplot as plt
from pytest import raises

from riptide import TimeSeries, ffa_search, Periodogram, save_json, load_json


def test_ffa_search():
    # NOTE: we chose a length long enough so that running the 
    # 'periodogram pruning' function was actually necessary 
    # (and thus the function gets properly covered by the tests)
    length = 200.0
    tsamp = 0.001
    period = 1.0
    amplitude = 20.0
    ts = TimeSeries.generate(length, tsamp, period, amplitude=amplitude)

    bins_min = 240
    bins_max = 260
    period_min = 0.8 * period
    period_max = 1.2 * period
    tsdr, pgram = ffa_search(
        ts, 
        period_min=period_min, period_max=period_max, 
        bins_min=bins_min, bins_max=bins_max
    )

    # check trial periods are increasing
    assert all(np.maximum.accumulate(pgram.periods) ==  pgram.periods)
    assert pgram.snrs.shape == (len(pgram.periods), len(pgram.widths))
    assert pgram.metadata == ts.metadata == tsdr.metadata
    assert pgram.tobs == length
    assert all(pgram.freqs == 1.0 / pgram.periods)

    # Test that running with deredden = False and already_normalised = True
    # returns a reference to the input TimeSeries (data left untouched)
    # This is how ffa_search() is called by the pipeline
    tsdr, pgram = ffa_search(
        ts, 
        period_min=period_min, period_max=period_max, 
        bins_min=bins_min, bins_max=bins_max,
        already_normalised=True, deredden=False
    )
    assert id(tsdr) == id(ts)


    ### Periodogram serialization ###
    with tempfile.NamedTemporaryFile(suffix='.json') as f:
        save_json(f.name, pgram)
        f.flush()
        pgram_copy = load_json(f.name)
        assert np.allclose(pgram.snrs, pgram_copy.snrs)
        assert np.allclose(pgram.periods, pgram_copy.periods)
        assert np.allclose(pgram.widths, pgram_copy.widths)
        assert pgram.metadata == pgram_copy.metadata


    ### Periodogram plotting ###
    plt.switch_backend('Agg')
    fig = plt.figure(figsize=(20, 5), dpi=100)
    pgram.plot()
    with tempfile.NamedTemporaryFile(suffix='.png') as fobj:
        plt.savefig(fobj.name)
        plt.close(fig)

    # Same with iwidth = 0
    fig = plt.figure(figsize=(20, 5), dpi=100)
    pgram.plot(iwidth=0)
    with tempfile.NamedTemporaryFile(suffix='.png') as fobj:
        plt.savefig(fobj.name)
        plt.close(fig)