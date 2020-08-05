from schema import Schema, Use, Optional, And, Or


class InvalidSearchRange(Exception):
    pass


class InvalidPipelineConfig(Exception):
    pass


def strictly_positive(x):
    return x > 0


VALID_FORMATS = ('presto', 'sigproc')


SEARCH_RANGE_SCHEMA = Schema({
    'name': str,

    'ffa_search': {
        'period_min': And(Use(float), strictly_positive, error="period_min must be a number > 0"),
        'period_max': And(Use(float), strictly_positive, error="period_max must be a number > 0"),
        'bins_min': And(int, strictly_positive, error="bins_min must be an int > 0"),
        'bins_max': And(int, strictly_positive, error="bins_max must be an int > 0"),
        Optional('fpmin'): And(int, strictly_positive, error="fpmin must be an int > 0"),
        Optional('wtsp'): And(
            Use(float), lambda x: x > 1, error="wtsp must be a number > 1"),
        Optional('ducy_max'): And(
            float, lambda x: 0 < x < 1, error='ducy_max must be strictly between 0 and 1'),
        },

    'find_peaks': {
        Optional('smin'): And(
            Use(float), strictly_positive, error="smin must be a number > 0"),
        Optional('segwidth'): And(
            Use(float), strictly_positive, error="segwidth must be a number > 0"),
        Optional('nstd'): And(
            Use(float), strictly_positive, error="nstd must be a number > 0"),
        Optional('minseg'): And(
            int, strictly_positive, error="minseg must be an int > 0"),
        Optional('polydeg'): And(
            Use(float), strictly_positive, error="polydeg must be a number > 0"),
        Optional('clrad'): Or(
            And(Use(float), strictly_positive), None, error="clrad must be a number > 0"),
        },

    'candidates': {
        'bins': And(int, strictly_positive, error="candidates.bins must be an int > 0"),
        'subints': And(int, strictly_positive, error="candidates.subints must be an int > 0"),
        },
})


PIPELINE_CONFIG_SCHEMA = Schema({
    'processes': And(int, strictly_positive, error="processes must be an int > 0"),

    'data': {
        'format': Schema(
            lambda x: x in VALID_FORMATS, 
            error=f"format must be one of {VALID_FORMATS}"),
        'fmin': Or(
            And(Use(float), strictly_positive), None,
            error="fmin must be a number > 0 or null/blank"),
        'fmax': Or(
            And(Use(float), strictly_positive), None,
            error="fmax must be a number > 0 or null/blank"),
        'nchans': Or(
            And(int, strictly_positive), None,
            error="nchans must be a number > 0 or null/blank"),
        },

    'dmselect': {
        'min': Or(Use(float), None, error="Minimum DM must be a number or null/blank"),
        'max': Or(Use(float), None, error="Maximum DM must be a number or null/blank"),
        'dmsinb_max': Or(
            strictly_positive, None, 
            error="dmsinb_max must be a number > 0 or null/blank"),
        },

    'dereddening': {
        'rmed_width': Schema(strictly_positive, error="rmed_width must be a number > 0"),
        'rmed_minpts': Schema(strictly_positive, error="rmed_minpts must be a number > 0")
        },

    'ranges': [SEARCH_RANGE_SCHEMA],

    'clustering': {
        'radius': Schema(strictly_positive, error="clustering radius must be a number > 0"),
        },

    'harmonic_flagging': {
        'denom_max': And(int, strictly_positive, error="denom_max must be an int > 0"),
        'phase_distance_max': And(
            Use(float), strictly_positive, error="phase_distance_max must be a number > 0"),
        'dm_distance_max': And(
            Use(float), strictly_positive, error="dm_distance_max must be a number > 0"),
        'snr_distance_max': And(
            Use(float), strictly_positive, error="snr_distance_max must be a number > 0"),
        },

    'candidate_filters': {
        'dm_min': Or(Use(float), None, error='Candidate dm_min must be a float or null/blank'),
        'snr_min': Or(Use(float), None, error='Candidate snr_min must be a float or null/blank'),
        'remove_harmonics': Or(
            bool, None, error='remove_harmonics must be a boolean or null/blank'),
        'max_number': Or(
            And(int, strictly_positive), None, 
            error='Candidate max_number must be an int > 0 or null/blank'),
        },

    'plot_candidates': Schema(bool, error='plot_candidates must be a boolean'),
})


def validate_range(rg, tsamp_max):
    """ """
    # NOTE: In general, we leave the pipeline code to raise the exceptions,
    # except if it takes too long for it to detect them; for example, if the number of candidate
    # bins is too large, we don't want to wait until the candidate building stage to realize this.
    period_min = rg['ffa_search']['period_min']
    period_max = rg['ffa_search']['period_max']
    bins_min = rg['ffa_search']['bins_min']
    cand_bins = rg['candidates']['bins']

    if bins_min * tsamp_max > period_min:
        raise InvalidSearchRange(
            f"Search range {period_min:.3e} to {period_max:.3e} seconds: requested phase "
            "resolution is too high w.r.t. coarsest input time series "
            f"(tsamp = {tsamp_max:.3e} seconds). Use smaller bins_min or larger period_min.") 

    if cand_bins * tsamp_max > period_min:
        raise InvalidSearchRange(
            f"Search range {period_min:.3e} to {period_max:.3e} seconds: "
            f"cannot fold candidates with such high resolution ({cand_bins:d} bins). "
            f"The coarsest input time series ({tsamp_max:.3e} seconds) does not allow it")  


def validate_ranges_contiguity(ranges):
    """ """
    for a, b in zip(ranges[:-1], ranges[1:]):
        period_max_a = a['ffa_search']['period_max']
        period_min_b = b['ffa_search']['period_min']
        if not period_max_a == period_min_b:
            raise InvalidSearchRange(
                "Search ranges are not either non-contiguous, or not ordered by increasing trial "
                f"period (period_max ({period_max_a:.6e}) != next period_min ({period_min_b:.6e})")


def validate_ranges(ranges, tsamp_max):
    """ 
    Check that the search ranges are valid. Raise an exception if not.

    Parameters
    ----------
    ranges : list of dict
        Search ranges read from the pipeline configuration file
    tsamp_max : float
        Maximum sampling interval of the TimeSeries to process

    Raises
    ------
    InvalidSearchRange
    """
    for rg in ranges:
        validate_range(rg, tsamp_max)
    validate_ranges_contiguity(ranges)
    

def validate_pipeline_config(conf):
    """
    Validate pipeline configuration dictionary and raise an error if it is
    incorrect. This function only checks the format of the config and 
    the data types.

    Parameters
    ----------
    conf : dict
        Configuration dictionary loaded from the pipeline config file

    Returns
    -------
    validated : dict
        Validated configuration dictionary. Some data types may have been
        changed (e.g. into to float, or float to int when both are allowed
        for a config parameter).

    Raises
    ------
    InvalidPipelineConfig
    """
    try:
        validated = PIPELINE_CONFIG_SCHEMA.validate(conf)
    except Exception as ex:
        # Suppress long and confusing exception chain caused by schema library
        raise InvalidPipelineConfig(str(ex)) from None
    return validated
