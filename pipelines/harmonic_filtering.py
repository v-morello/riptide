from fractions import Fraction

def test_harmonic_relationship(A, B, max_denominator=30, snr_tol=1.50, dft_bins_tol=0.1):
    """ Test if the DetectionCluster B is a harmonic of DetectionCluster A. Note
    that the test will always evaluate to False if B has a higher S/N than A.

    Parameters:
    -----------
        A: DetectionCluster
            The brighter (assumed fundamental) detection.
        B: DetectionCluster
            The fainter (assumed harmonic) detection.
        max_denominator: int, optional
            When approximating the ratio between the frequencies of A and B
            with a rational number p/q, this is the largest value the
            denominator can take. (default: 30)
        snr_tol: float, optional
            if sA and sB are the signal-to-noise ratios of A and B, then one
            can expect to have sA / sB = sqrt(pq). For B to be considered a
            harmonic of A, we must have: sqrt(pq) <= snr_tol * sA/sB
            snr_tol should be at least 1.0, and higher values cause more
            distant harmonics to be flagged, sometimes abusively.
            (default: 1.50)
        dft_bins_tol: float, optional
            The maximum distance (in DFT bins) allowed between vA and p/q x vB,
            where vA and vB are the frequencies of A and B. If
            abs(vA - p/q x vB) exceeds that value, A and B are considered NOT
            harmonically related. dft_bins_tol should be definitely < 1.0,
            and higher values may cause some candidates to be abusively flagged
            as harmonics. (default: 0.10)

    Returns:
    --------
        test_result: bool
            Whether or not A and B are harmonically related
        fraction: fractions.Fraction
            The ratio between the frequencies of A and B (vA / vB)
    """
    tobs = A.top_detection.metadata['tobs']
    max_denominator = int(max_denominator)

    # Frequencies of A and B
    va = A.top_detection.freq
    vb = B.top_detection.freq

    # DFT bin indices of A and B
    da = va * tobs
    db = vb * tobs

    # signal-to-noise ratios
    sa = A.top_detection.snr
    sb = B.top_detection.snr

    freq_ratio = va / vb
    snr_ratio = sa / sb

    # NOTE: float() cast is to avoid period_ratio to be np.float32
    # which can't be passed as an argument to Fraction()
    fraction = Fraction(float(freq_ratio)).limit_denominator(max_denominator)
    condition_freq = abs(da - fraction * db) < dft_bins_tol

    # For a freq ratio expressed as the irreducible fraction a/b,
    # The expected S/N ratio is sqrt(ab)
    expected_snr_ratio = (fraction.numerator * fraction.denominator) ** 0.5
    condition_snr = (sa >= sb) and (expected_snr_ratio <= snr_ratio * snr_tol)


    test_result = condition_snr and condition_freq
    return test_result, fraction
