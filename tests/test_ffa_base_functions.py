# Test the following:
# ffa1(), ffa2()
# ffafreq(), ffaprd()
import numpy as np
from pytest import raises

from riptide import ffa1, ffa2, ffafreq, ffaprd

# Manually calculated 8x8 test case
# This is expected to be invariant under phase rotation
# and invariant when appending columns made of zeros
FFA_IN_88 = np.array([
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 1]
]).astype(np.float32)

FFA_OUT_88 = np.array([
    [0, 0, 0, 0, 0, 0, 0, 8],
    [0, 0, 0, 0, 0, 0, 4, 4],
    [0, 0, 0, 0, 0, 2, 4, 2],
    [0, 0, 0, 0, 2, 2, 2, 2],
    [0, 0, 0, 1, 2, 2, 2, 1],
    [0, 0, 1, 2, 1, 1, 2, 1],
    [0, 1, 1, 1, 2, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1]
]).astype(np.float32)


def test_transforms():
    for shift in range(8):
        X = np.roll(FFA_IN_88, shift, axis=1)
        m, p = X.shape
        truth = np.roll(FFA_OUT_88, shift, axis=1)
        Y = ffa2(X)
        Z = ffa1(X.ravel(), p)
        assert np.allclose(Y, truth)
        assert np.allclose(Z, truth)

    for extra_cols in range(8):
        X = np.hstack([
            FFA_IN_88, np.zeros(shape=(8, extra_cols))
        ])
        m, p = X.shape
        truth = np.hstack([
            FFA_OUT_88, np.zeros(shape=(8, extra_cols))
        ])
        Y = ffa2(X)
        Z = ffa1(X.ravel(), p)
        assert np.allclose(Y, truth)
        assert np.allclose(Z, truth)

    # Test ffa2() errors
    with raises(ValueError):
        X = np.zeros(4) # 1D array
        ffa2(X)

    # Test ffa1() errors
    # Input not 1D
    with raises(ValueError):
        X = np.zeros((4, 4)) # 2D array
        ffa1(X, 4)
    
    X = np.zeros(10)

    with raises(ValueError):
        ffa1(X, X.size + 1) # p too large

    with raises(ValueError):
        ffa1(X, 4.0) # p not int
    

def test_ffafreq():
    # NOTE: ffaprd() simply does 1.0 / ffafreq(), so we only need to properly
    # cover all ffafreq() corner cases

    ### Correctness
    # See the paper for the formula: https://arxiv.org/abs/2004.03701
    m = 42
    p = 127
    N = m * p
    dt = np.pi / 1000.0 # pi milliseconds sampling time, why not
    
    s = np.arange(m, dtype=float)
    true_periods = p**2 / (p - s/(m-1.0)) * dt
    true_freqs = 1.0 / true_periods

    freqs = ffafreq(N, p, dt=dt)
    periods = ffaprd(N, p, dt=dt)
    assert np.allclose(freqs, true_freqs)
    assert np.allclose(periods, true_periods)

    # Special case where the FFA input data has only one signal period (m = 1)
    assert ffafreq(p, p, dt=dt)[0] == 1.0 / (p * dt)

    ### Errors
    with raises(ValueError):
        ffafreq(0, p, dt=dt) # N <= 0

    with raises(ValueError):
        ffafreq(np.pi, p, dt=dt) # N != int

    with raises(ValueError):
        ffafreq(N, 1, dt=dt) # p <= 1

    with raises(ValueError):
        ffafreq(N, np.pi, dt=dt) # p != int

    with raises(ValueError):
        ffafreq(N, N + 1, dt=dt) # p > n
    
    with raises(ValueError):
        ffafreq(N, p, dt=0) # dt <= 0
