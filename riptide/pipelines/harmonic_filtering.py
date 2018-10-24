import operator
import itertools
import json
import uuid
import math
from fractions import Fraction


class CandidateParameters(dict):
    """ Dictionary that wraps basic candidate parameters, just for the purpose
    of harmonic filtering. """
    def __init__(self, items):
        super(CandidateParameters, self).__init__(items)
        self["freq"] = 1.0 / self["period"]
        self["width"] = self["ducy"] * self["period"]
        self["is_harmonic"] = items.get("is_harmonic", False)
        self["fundamental_index"] = items.get("fundamental_index", None)
        self["fraction"] = items.get("fraction", None)

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return json.dumps(self, indent=4, default=str)


def harmonic_test(F, H, fmin=1181.0, fmax=1581.0, tobs=536.0, max_denom=50, max_distance=1.0, snr_tol=1.5):
    """ Calculate several distance measures between two candidates F and H that are suspected to
    be respectively the fundamental and one of its harmonics. """
    kdm = 4.148808e3
    r = {} # output dictionary

    # Fraction closest to the ratio nu_F / nu_H
    # or P_H / P_F
    ratio = F["freq"] / H["freq"]
    if ratio < 1.0 / max_denom:
        max_denom = int(math.ceil(1.0 / ratio))
    frac = Fraction(ratio).limit_denominator(max_denom)

    # Phase distance
    r["fraction"] = frac
    r["delta_freq"] = abs(H["freq"] * frac - F["freq"])
    r["delta_phase"] = r["delta_freq"] * tobs
    r["phase_distance"] = r["delta_phase"] / F["ducy"]

    # The DM distance is the dispersion delay between the two candidates
    # expressed in units of pulse width
    r["delta_dm_delay"] = kdm * (F["dm"] - H["dm"]) * (fmin**-2 - fmax**-2)
    r["dm_distance"] = abs(r["delta_dm_delay"]) / min(F["width"], H["width"])

    # Total distance
    r["distance"] = (r["dm_distance"]**2 + r["phase_distance"]**2) ** 0.5

    # Ratio between actual and expected S/N of the (suspected) harmonic
    r["expected_harmonic_snr"] = F["snr"] / (frac.numerator * frac.denominator) ** 0.5
    r["snr_ratio"] = max(H["snr"], r["expected_harmonic_snr"]) / min(H["snr"], r["expected_harmonic_snr"])
    r["result"] = r["distance"] < max_distance and r["snr_ratio"] < snr_tol
    return r


def flag_harmonics(params, fmin=1181.0, fmax=1581.0, tobs=536.0, max_denom=50, max_distance=1.0, snr_tol=1.5):
    """ Among a list of candidates, flag those who are harmonics of another. 

    Parameters
    ----------
    params: list or iterable
        List (or iterable) of dictionaries containing basic candidate
        parameters. Each dict must have at least the following keys:
        "period", "ducy", "dm", "snr".
    fmin: float
        Bottom observing frequency in MHz
    fmax: float
        Top observing frequency in MHz
    tobs: float
        The integration time.
    max_denom: int
        When trying to find the closest rational fraction to the ratio of
        fundamental and harmonic periods, this is the largest denominator
        that will be tested.
    max_distance: float
        Distance between candidates below which they can be declared
        harmonically related.
    snr_tol: float
        Maximum ratio between expected and true S/N ratio of the harmonic.
        Expected S/N of the harmonic is F.snr / sqrt(ab), where a/b is the
        rational number closest to the ratio of H and F's periods.

    Returns
    -------
    cparams: list
        Copy of "params", where each candidate dictionary has extra attributes:
            "is_harmonic": whether or not it is a harmonic of another candidate in the lit (called fundamental)
            "fundamental_index": index of the associated fundamental in the input list if "is_harmonic" is True, otherwise None.
            "fraction": rational fraction closest to P_harmonic / P_fundamental if "is_harmonic" is True, otherwise None.
    """
    cparams = list(map(CandidateParameters, params))
    uid = uuid.uuid4()
    for i, c in enumerate(cparams):
        c[uid] = i

    sorted_cparams = sorted(cparams, key=operator.itemgetter('snr'), reverse=True)
    for fund, harm in itertools.combinations(sorted_cparams, 2):
        if fund["is_harmonic"] or harm["is_harmonic"]:
            continue
        test = harmonic_test(
            fund, harm, 
            fmin=fmin, fmax=fmax, 
            tobs=tobs, max_denom=max_denom, 
            max_distance=max_distance, snr_tol=snr_tol)
        if test["result"]:
            harm["is_harmonic"] = True
            harm["fundamental_index"] = fund[uid]
            harm["fraction"] = test["fraction"]
    
    for c in cparams:
        del c[uid]

    return cparams