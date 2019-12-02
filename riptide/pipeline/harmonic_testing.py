import math
from fractions import Fraction
import logging


log = logging.getLogger('riptide.pipeline.harmonic_filter')


def hdiag(F, H, tobs, fmin, fmax, denom_max=100):
    """
    Calculate a number of diagnostic values to evaluate whether two sets of 
    candidate parameters are harmonically related.

    Parameters
    ----------
    F: object
        Postulated fundamental candidate parameters. Can be of any type, but
        it must have the following members:
        - freq: Frequency in Hz
        - snr: signal to noise ratio
        - ducy: duty cycle, i.e. pulse width / period
        - dm: dispersion measure in pc cm^{-3}
    H: object
        Postulated harmonic candidate parameters. Same requirements as above.
    tobs: float
        Integration time in seconds
    fmin: float
        Bottom effective observing frequency in Hz
    fmax: float
        Top effective observing frequency in Hz
    denom_max: int, optional
        Maximum allowed denominator of the harmonic fraction by which the 
        frequencies of both candidates are related. It must be limited,
        otherwise there is always a rational fraction arbitrarily close to
        the ratio of their frequencies. (default: 100)

    Returns
    -------
    dvals: dict
        Dictionary of diagnostic values
    """
    if not fmax > fmin:
        raise ValueError
    if not tobs > 0.0:
        raise ValueError

    def width(X):
        return X.ducy / X.freq

    # check_attrs(F)
    # check_attrs(H)
    slow, fast = sorted((F, H), key=lambda x: x.freq)

    # Find the fraction closest to the ratio of frequencies
    fraction = Fraction(fast.freq / slow.freq).limit_denominator(denom_max)

    # Phase distance
    # We evaluate how close are the frequencies fast.freq and fraction x slow.freq
    phase_absdiff_turns = abs(fraction * slow.freq - fast.freq) * tobs
    phase_distance = phase_absdiff_turns / fast.ducy

    # NOTE: in the end we want to return the fraction H.freq / F.freq
    # NOTE: fraction = 2 means that H is the second harmonic of F in the Fourier sense
    # (H has a frequency twice that of F)
    if H == slow:
        fraction = 1 / fraction

    # Test DM difference
    kdm = 4.15e3
    dm_absdiff = abs(F.dm - H.dm)
    dm_delay_absdiff = dm_absdiff * kdm * abs(fmin**-2 - fmax**-2)
    dm_distance = dm_delay_absdiff / min(width(F), width(H))

    # Test S/N difference
    harmonic_snr_expected = F.snr / (fraction.numerator * fraction.denominator) ** 0.5
    snr_distance = abs(H.snr - harmonic_snr_expected)

    return {
        'fraction': fraction,
        'phase_absdiff_turns': phase_absdiff_turns,
        'phase_distance' : phase_distance,
        'dm_absdiff': dm_absdiff,
        'dm_delay_absdiff': dm_delay_absdiff,
        'dm_distance': dm_distance,
        'harmonic_snr_expected': harmonic_snr_expected,
        'snr_distance': snr_distance,
    }


def htest(F, H, tobs, fmin, fmax, denom_max=100, phase_distance_max=1.0, dm_distance_max=3.0, snr_distance_max=3.0):
    """
    Test whether two sets of candidate parameters are harmonically related.
    The code first finds the closest rational fraction p/q to the ratio 
    H.freq / F.freq and then tests whether H is the plausible p/q-th harmonic
    of F. The method is *purposely* designed to under-flag rather than 
    over-flag, noting also that pipeline users can decide not to remove
    harmonics from the final candidate list.

    Parameters
    ----------
    F: object
        Postulated fundamental candidate parameters. Can be of any type, but
        it must have the following members:
        - freq: Frequency in Hz
        - snr: signal to noise ratio
        - ducy: duty cycle, i.e. pulse width / period
        - dm: dispersion measure in pc cm^{-3}
    H: object
        Postulated harmonic candidate parameters. Same requirements as above.
    tobs: float
        Integration time in seconds
    fmin: float
        Bottom effective observing frequency in Hz
    fmax: float
        Top effective observing frequency in Hz
    denom_max: int, optional
        Maximum allowed denominator of the harmonic fraction by which the 
        frequencies of both candidates are related. It must be limited,
        otherwise there is always a rational fraction arbitrarily close to
        the ratio of their frequencies. (default: 100)
    phase_distance_max: float
        Upper bound on the phase delay (in number of pulse widths) accrued 
        over 'tobs' seconds between the signal H and the hypothesised harmonic
        p/q x F. A value of 1.0 means that the harmonic relationship is 
        credible only of both trains of pulses within one pulse width of each
        other. This the proper way to measure if the frequencies H.freq and 
        p/q x F.freq are significantly close. (default: 1.0)
    dm_distance_max: float
        Upper bound on the difference between dispersion delays (expressed in 
        pulse widths) across the observing band associated to the DMs of 
        F and H. (default: 3.0)    
    snr_distance_max: float
        Upper bound on the absolute difference between the true S/N of H and
        the S/N that it should have if it was the p/q harmonic of F.
        The expected S/N of H in this case is F.snr x sqrt(p x q).
        (default: 3.0)

    Returns
    -------
    related: bool
        True if H is reported to be a harmonic of F, which happens only if all
        of the following conditions are met:
        - H.freq and p/q x F.freq are close enough
        - H.dm and F.dm are close enough
        - H.snr and F.snr x sqrt(p x q) are close enough
    fraction: fractions.Fraction
        The rational fraction p/q closest to H.freq / F.freq
    """
    dvals = hdiag(F, H, tobs, fmin, fmax, denom_max=denom_max)
    related = \
        dvals['phase_distance'] <= phase_distance_max and \
        dvals['dm_distance'] <= dm_distance_max and \
        dvals['snr_distance'] <= snr_distance_max
    fraction = dvals['fraction']
    return related, fraction
