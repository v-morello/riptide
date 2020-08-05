from schema import Schema, Use, Optional, And, Or


def strictly_positive(x):
    return x > 0


VALID_FORMATS = ('presto', 'sigproc')


# TODO: finish this
SEARCH_RANGE_SCHEMA = Schema({
    'name': str,

    'ffa_search': {
        'period_min': float,
        'period_max': float,
        'bins_min': int,
        'bins_max': int,
        Optional('fpmin'): int,
        Optional('wtsp'): float,
        Optional('ducy_max'): float,
        },

    'find_peaks': {
        Optional('smin'): float,
        Optional('segwidth'): float,
        Optional('nstd'): float,
        Optional('minseg'): float,
        Optional('polydeg'): float,
        Optional('clrad'): float,
        },

    'candidates': {
        'bins': int,
        'subints': int
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


def validate_config(conf):
    """
    Validate pipeline configuration dictionary and raise an error if it is
    incorrect.

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
    """
    valid = PIPELINE_CONFIG_SCHEMA.validate(conf)
    # TODO
    # Further checks to perform:
    # * Ranges are contiguous
    # * Candidate bins is not too large
    # In general, leave the pipeline code raise errors, except when it takes too long to wait
    # for the error to be raised (e.g. it is raised only at Candidate building stage)
    return valid


if __name__ == '__main__':
    import yaml

    with open('/home/vince/repositories/riptide/riptide/pipeline/config/example.yaml', 'r') as fobj:
        conf = yaml.safe_load(fobj)

    print(
        yaml.dump(validate_config(conf), indent=4)
    )