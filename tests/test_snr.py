from pytest import raises
import numpy as np

from riptide import boxcar_snr


def test_errors():
    cols = 32
    data = np.zeros(cols, dtype='float32')
    
    # No widths < 1
    with raises(ValueError):
        boxcar_snr(data, [0, 1])

    # No widths >= columns
    with raises(ValueError):
        boxcar_snr(data, [1, cols])

    # stdnoise must be > 0
    with raises(ValueError):
        boxcar_snr(data, [1, 2], stdnoise=-42.0)


def test_output_dims():
    cols = 32
    widths = [1, 2, 3, 5]

    # 1D input
    data = np.zeros(cols, dtype='float32')
    snr = boxcar_snr(data, widths)
    assert snr.ndim == 1
    assert snr.size == len(widths)

    # 2D input
    rows = 4
    data = np.zeros((rows, cols), dtype='float32')
    snr = boxcar_snr(data, widths)
    assert snr.ndim == 2
    assert snr.shape == (rows, len(widths))

    # 3D input
    layers = 3
    data = np.zeros((layers, rows, cols), dtype='float32')
    snr = boxcar_snr(data, widths)
    assert snr.ndim == 3
    assert snr.shape == (layers, rows, len(widths))


def test_phase_rotation_invariance():
    rows = 4
    cols = 32
    widths = [1, 2, 5, 11, 18, 31]

    data = np.random.normal(size=rows*cols).reshape(rows, cols).astype('float32')
    snr_ref = boxcar_snr(data, widths)

    for shift in range(1, cols + 1):
        snr = boxcar_snr(np.roll(data, shift, axis=-1), widths)
        assert np.allclose(snr, snr_ref)


def test_values():
    n = 64
    widths = np.arange(1, n)
    data = np.zeros(n, dtype='float32')

    for w in range(1, n): # w = true width
        data[:w] = 1.0
        snr = boxcar_snr(data, widths)
        assert snr.argmax() == w - 1

        # Best S/N should be obtained for trial width = true width = w
        # h and b are the height and baseline value of a boxcar filter
        # with zero mean and unit square sum
        h = np.sqrt((n - w) / (n * w))
        b = -w / (n - w) * h
        expected_best_snr = w * h
        assert np.allclose(snr.max(), expected_best_snr)
