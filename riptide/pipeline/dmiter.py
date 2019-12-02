import logging

import numpy as np

from riptide import Metadata, TimeSeries


log = logging.getLogger('riptide.dmiter')


class DMIterator(object):
    """
    """
    METADATA_LOADERS = {
        'sigproc': Metadata.from_sigproc,
        'presto': Metadata.from_presto_inf
    }

    TIMESERIES_LOADERS = {
        'sigproc': TimeSeries.from_sigproc,
        'presto': TimeSeries.from_presto_inf
    }

    KDM = 1.0 / 2.41e-4

    def __init__(self, dm_min=0.0, dm_max=1000.0, dmsinb_max=None, wmin=1.0e-3, fmin=1181.0, fmax=1581.0, nchans=1024):
        """
        dmsinb_max : float or None, optional
            If a float is specified, further limit the maximum trial DM to: abs(dmsinb_max * sin(glat))
            If None, do not apply any latitude-dependent trial DM cap
        """
        self.dm_min = float(dm_min)
        self.dm_max = float(dm_max)
        self.dmsinb_max = float(dmsinb_max) if dmsinb_max else None
        self.wmin = float(wmin)

        self.fmin = float(fmin)
        self.fmax = float(fmax)
        self.nchans = int(nchans)

        self.fcenter = 0.5 * (self.fmin + self.fmax)
        self.cw = abs(self.fmin - self.fmax) / self.nchans

        self.mdloader = None # Metadata loader
        self.tsloader = None # TimeSeries loader
        self.mdlist = []     # list of Metadata objects read from input files

    def tobs_median(self):
        if not self.mdlist:
            raise ValueError("No TimeSeries appear to have been loaded")
        return np.median([meta['tobs'] for meta in self.mdlist])

    def dm_step(self, dm):
        """ 
        Get the ideal DM step at given DM. This is used to select the optimal
        subset of available DM trials.
        """
        fbot = self.fcenter - self.cw/2.0
        ftop = self.fcenter + self.cw/2.0
        tsmear = self.KDM * dm * (fbot**-2 - ftop**-2)
        weff = max(self.wmin, tsmear)
        step = 2 * weff / (self.KDM * (self.fmin**-2 - self.fmax**-2))
        return step

    def get_metadata_loader(self, fmt):
        choices = list(self.METADATA_LOADERS.keys())
        if not fmt in choices:
            raise ValueError("fmt must be one of {}".format(choices))
        return self.METADATA_LOADERS[fmt]

    def get_timeseries_loader(self, fmt):
        choices = list(self.TIMESERIES_LOADERS.keys())
        if not fmt in choices:
            raise ValueError("fmt must be one of {}".format(choices))
        return self.TIMESERIES_LOADERS[fmt]

    def get_dm_trial(self, dm):
        """
        Get the TimeSeries with the given DM, if it exists. Raises ValueError
        otherwise.
        """
        try:
            meta = next(m for m in self.mdlist if m['dm'] == dm)
        except StopIteration:
            raise ValueError(f"No DM trial with specified DM = {dm}")
        return self.tsloader(meta['fname'])

    def prepare(self, files, fmt='presto'):
        """
        Prepare the DM iterator to iterate over the given list of files,
        selecting only those within the right DM range, and ensuring that
        the selected trials are spaced as much as possible.

        Parameters
        ----------
        files: list
        fmt: str
        """
        self.mdloader = self.get_metadata_loader(fmt)
        self.tsloader = self.get_timeseries_loader(fmt)

        mdlist = []
        for fname in files:
            try:
                mdlist.append(self.mdloader(fname))
            except Exception as err:
                log.warning("Failed to load {!r}: {}".format(fname, err))
        
        self.mdlist = mdlist
        log.info("Read metadata of {} input files".format(len(mdlist)))

        # Remove some trials if they are spaced too closely
        self.prune()

        selected_dm_trials = [meta['dm'] for meta in self.mdlist]
        log.info("DM trials selected for processing: {}".format(len(self.mdlist)))
        log.info("List of selected DM trials: {}".format(selected_dm_trials))

    def prune(self):
        log.info("Pruning DM trial list to achieve least expensive DM space coverage")

        def within_dm_range(meta):
            dm = meta['dm']
            result = self.dm_min <= dm <= self.dm_max
            if self.dmsinb_max is not None:
                gb_rad = meta['skycoord'].galactic.b.rad
                gal_limit = abs(self.dmsinb_max / np.sin(gb_rad))
                below_gal_limit = dm <= gal_limit
                result = result and below_gal_limit
            return result

        # Remove trials out of range
        self.mdlist = [meta for meta in self.mdlist if within_dm_range(meta)]
        log.debug(f"Retained {len(self.mdlist)} DM trials within specified DM range (min = {self.dm_min:.3f}, max = {self.dm_max:.3f}, dmsinb_max = {self.dmsinb_max})")

        # Sort by increasing DM
        self.mdlist = sorted(self.mdlist, key=lambda meta: meta['dm'])

        def get_dm(i):
            """ Get DM of trial index i """
            n = len(self.mdlist)
            if i < 0:
                return get_dm(0)
            elif i >= n:
                return self.dm_max
            else:
                return self.mdlist[i]['dm']

        def get_dist(i, j):
            return abs(get_dm(i) - get_dm(j))

        if not self.mdlist:
            return

        # DM trial indices to keep: always keep first trial
        keep = [0]
        log.debug("DM = {:.3f}: Keep first trial".format(get_dm(0)))

        for icur in range(1, len(self.mdlist)):
            ilast = keep[-1]
            dmcur = get_dm(icur)
            step = self.dm_step(get_dm(icur)) # ideal DM step at current DM

            # Current DM trial index i can be pruned
            # if abs(dm_last - dm_{i+1}) <= step
            dist = get_dist(ilast, icur+1)
            if dist <= step:
                log.debug("DM = {:.3f}, step = {:.3f}: Discard".format(dmcur, step))
                continue
            else:
                log.debug("DM = {:.3f}, step = {:.3f}: Keep".format(dmcur, step))
                keep.append(icur)
        self.mdlist = [self.mdlist[ii] for ii in keep]
    
    def iterate(self, chunksize=1):
        chunk = []
        for meta in self.mdlist:   
            ts = self.tsloader(meta['fname'])
            chunk.append(ts)

            if len(chunk) == chunksize:
                yield chunk
                chunk = []

        # Don't forget to yield any non-empty, incomplete last chunk
        if chunk:
            yield chunk
