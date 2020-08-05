import logging
import numpy as np
from numpy import sin, radians
from riptide import TimeSeries, Metadata


log = logging.getLogger('riptide.pipeline.dmiter')


# This is the standard "rounded value" of the dispersion constant in use by pulsar astronomers
# (page 129 of Manchester and Taylor 1977)
KDM = 1.0 / 2.41e-4


def select_dms(trial_dms, dm_start, dm_end, fmin, fmax, nchans, wmin):
    """
    Trial DMs are selected such that the amount of pulse broadening caused by
    DM error is no greater than max(tsmear, wmin) where tsmear is the amount
    of pulse broadening caused by intra-channel smearing, and wmin the minimum
    pulse width being searched for.

    Parameters
    ----------
    trial_dms : ndarray or iterable
    dm_start : float
        Start DM
    dm_end : float
        End DM
    [...] TODO [...]

    Returns
    -------
    selected : ndarray
    """
    trial_dms = np.asarray(trial_dms)
    mask = (trial_dms >= dm_start) & (trial_dms <= dm_end)
    trial_dms = trial_dms[mask]
    trial_dms = np.sort(trial_dms)

    if not trial_dms.size:
        raise ValueError(f"No trial DMs between {dm_start:.4f} and {dm_end:.4f}")

    # tdisp = kdisp * dm
    kdisp = KDM * (fmin**-2 - fmax**-2)
    cw = (fmax - fmin) / nchans
    fmid = (fmax + fmin) / 2.0

    # tsmear = ksmear * dm
    ksmear = KDM * ((fmid-cw/2)**-2 - (fmid+cw/2)**-2)

    # Coverage radius (in DM space) of every trial DM
    # Within this radius, the total smearing time is <= wmin
    radii = np.maximum(wmin, ksmear * trial_dms) / kdisp

    def dm_gap(i, j): # assumes i <= j
        return (trial_dms[j] - radii[j]) - (trial_dms[i] + radii[i])

    def largest_in_range(i):
        for j in range(i, trial_dms.size):
            if dm_gap(i, j) > 0:
                return max(j - 1, i + 1)
        return max(j, i + 1)

    icur = 0
    selected = [trial_dms[icur]]

    while True:
        inext = largest_in_range(icur)

        if inext >= trial_dms.size:
            break

        if inext == icur + 1 and dm_gap(icur, inext) > 0:
            log.warning(
                f"The step from trial DM {trial_dms[icur]:.4f} should not exceed "
                f"{2 * radii[icur]:.4f}, "
                f"but the next available trial DM lies farther, at {trial_dms[inext]:.4f}")
        selected.append(trial_dms[inext])
        icur = inext
    return np.asarray(selected)



def get_band_params(meta, fmt='presto'):
    """
    Returns (fmin, fmax, nchans) given a metadata dictionary loaded from
    a specific file format.
    """
    if fmt == 'presto':
        fbot = meta['fbot']
        nchans = meta['nchan']
        ftop = fbot + nchans * meta['cbw']
        fmin = min(fbot, ftop)
        fmax = max(fbot, ftop)
    elif fmt == 'sigproc':
        raise ValueError("Cannot parse observing band parameters from data in sigproc format")
    else:
        raise ValueError(f"Unknown format: {fmt}")
    return fmin, fmax, nchans


def infer_band_params(metadata_list, fmt='presto'):
    """
    Read observing band parameters of all given Metadata objects, and check
    that they are all the same (otherwise, raise RuntimeError).
    Returns a tuple (fmin, fmax, nchans).
    """
    if not metadata_list:
        raise ValueError(
            "Cannot infer observing band parameters from empty metadata list. "
            "It appears no TimeSeries were passed as input."
            )
    params = [get_band_params(md, fmt=fmt) for md in metadata_list]
    if not all([params[0] == p for p in params]):
        raise RuntimeError(
            "Observing band parameters are NOT identical across all dedispersed time series")
    return params[0]


def get_galactic_coordnates(metadata_list):
    """
    Read galactic coordinates of all given Metadata objects, and check that
    they are all the same (otherwise, raise RuntimeError).
    Returns a float tuple (gl_deg, gb_deg).
    """
    def galc(md):
        coord = md['skycoord'].galactic
        return coord.l.deg, coord.b.deg

    ref = galc(metadata_list[0])
    if not all([galc(md) == ref for md in metadata_list]):
        raise RuntimeError("Coordinates are NOT identical across all dedispersed time series")
    return ref


class DMIterator(object):
    """
    Iterate through the minimum subset of DM trials to achieve DM space coverage.
    Observing band parameters (fmin, fmax, nchans) are inferred from the
    input files if possible, in which case any values passed as arguments are
    ignored. If the input files do not contain this information (example:
    SIGPROC dedispersed data), then fmin, fmax and nchans must all be specified,
    otherwise a ValueError is raised.

    Parameters
    ----------
    filenames : list of str
        List of TimeSeries file names
    dm_start : float or None
        Lowest DM to cover. If None, start at the smallest available DM trial.
    dm_end : float or None
        Highest DM to cover. If None, stop at the largest available DM trial.
    fmt : str
        Format of the TimeSeries objects
    wmin : float
        Minimum pulse width being searched for, in seconds
    fmin: float or None
        Minimum observing frequency, in MHz. Leave this to None unless the
        time series header / metadata does not contain this value.
    fmax: float or None
        Maximum observing frequency, in MHz. Leave this to None unless the
        time series header / metadata does not contain this value.
    nchans : int or None
        Number of observing channels. Leave this to None unless the
        time series header / metadata does not contain this value.
    """

    METADATA_LOADERS = {
        'sigproc': Metadata.from_sigproc,
        'presto': Metadata.from_presto_inf
    }

    # TODO: actually implement dmsinb_max
    def __init__(self, filenames, dm_start, dm_end, dmsinb_max=45.0, fmt='presto', wmin=1.0e-3,
                 fmin=None, fmax=None, nchans=None):
        mdloader = self.METADATA_LOADERS[fmt]
        self.metadata_list = [mdloader(fname) for fname in filenames]
        self.dm_start = float(dm_start) if dm_start is not None \
            else min(md['dm'] for md in self.metadata_list)
        self.dm_end = float(dm_end) if dm_end is not None \
            else max(md['dm'] for md in self.metadata_list)
        self.dmsinb_max = float(dmsinb_max) if dmsinb_max is not None else None
        self.fmt = fmt
        self.wmin = wmin

        # Apply DM sin |b| cap, if any
        gl_deg, gb_deg = get_galactic_coordnates(self.metadata_list)
        if self.dmsinb_max is not None:
            galactic_dm_cap = self.dmsinb_max / abs(sin(radians(gb_deg)))
            log.info(
                f"Applying DM|sin b| cap of {self.dmsinb_max:.4f}: "
                f"At b = {gb_deg:.2f} deg this means a max DM of {galactic_dm_cap:.4f}"
                )
            self.dm_end = min(self.dm_end, galactic_dm_cap)
        
        log.info(f"Selecting DM trials in the range {self.dm_start:.4f} to {self.dm_end:.4f}")

        # Try to infer band parameters from the data
        try:
            (self.fmin, self.fmax, self.nchans) = infer_band_params(self.metadata_list, fmt=fmt)
            log.info(
                "Inferred observing band parameters from input files: "
                f"fmin = {self.fmin:.3f}, fmax = {self.fmax:.3f}, nchans = {self.nchans:d}. "
                "Any manually specified values of fmin/fmax/nchans will be ignored."
                )
        except (ValueError, RuntimeError) as err:
            log.info(f"Could not infer observing band parameters from input files: {err!s}")
            log.info("Using manually specified band parameters instead")
            if any([param is None for param in (fmin, fmax, nchans)]):
                raise ValueError("You MUST specify: fmin, fmax, nchans")
            else:
                (self.fmin, self.fmax, self.nchans) = (fmin, fmax, nchans)
                log.info(
                    f"Using: fmin = {self.fmin:.3f}, "
                    f"fmax = {self.fmax:.3f}, nchans = {self.nchans:d}")

        self.metadata_dict = {meta['dm']: meta for meta in self.metadata_list}

        log.info(
            f"Selecting minimal trial DM subset to cover the DM range "
            f"{self.dm_start:.4f} to {self.dm_end:.4f}")
        self.selected_dms = select_dms(
            list(self.metadata_dict.keys()),
            self.dm_start, self.dm_end, self.fmin, self.fmax, self.nchans, self.wmin
            )

        log.info(
            f"Selected {len(self.selected_dms)} DM trials for processing: "
            f"{list(self.selected_dms)}")

    def iterate_filenames(self, chunksize=1):
        """
        Iterate through selected DM trial filenames in chunks of given size
        """
        chunk = []
        for dm in self.selected_dms:
            fname = self.metadata_dict[dm]['fname']
            chunk.append(fname)
            if len(chunk) == chunksize:
                yield chunk
                chunk = []
        if chunk: # yield any non-empty, incomplete last chunk
            yield chunk

    def get_filename(self, dm):
        return self.metadata_dict[dm]['fname']

    def tobs_median(self):
        return np.median([md['tobs'] for md in self.metadata_list])

    def tsamp_max(self):
        return max([md['tsamp'] for md in self.metadata_list])
