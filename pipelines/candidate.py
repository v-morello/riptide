##### Standard imports #####
import logging
import operator
import copy

##### Non-standard imports #####
import numpy as np
import h5py

##### Local imports #####
from .. import Metadata, SubIntegrations


def create_dpw_cube(detections, logger=None):
    """ Create a DM-Period-Width cube from a list of clustered Detections.

    Returns:
    --------
        dm_trials
        period_trials
        width_trials
        cube
    """
    if not logger:
        logger = logging.getLogger('Candidate')

    # Each detection contains a S/N versus Period and Width array
    # The width trials are guaranteed to be all the same, but the period
    # trials may differ.
    # We need to resample those to a common set of period trials
    pmin = max([det.period_trials[ 0] for det in detections])
    pmax = min([det.period_trials[-1] for det in detections])

    # Choose the number of common period trials 'ncpt'
    ncpt = max([len(det.period_trials) for det in detections])
    logger.info("Re-sampling PW planes with: pmin = {pmin:.9e}, pmax = {pmax:.9e}, num_periods = {ncpt:4d}".format(pmin=pmin, pmax=pmax, ncpt=ncpt))

    period_trials = np.linspace(pmin, pmax, ncpt)
    width_trials = detections[0].width_trials
    dm_trials = np.sort(np.asarray([det.dm for det in detections]))

    ndm = len(detections)
    nw = len(width_trials)
    cube = np.zeros(shape=(ndm, ncpt, nw))

    # Re-sample all P-W planes in increasing DM order
    # Create the final DM-Period-Width cube from the output
    for idm, det in enumerate(sorted(detections, key=operator.attrgetter('dm'))):
        pwp = det.snr_trials
        for iw, width in enumerate(width_trials):
            snr = pwp[:, iw]
            cube[idm, :, iw] = np.interp(period_trials, det.period_trials, snr)

    return dm_trials, period_trials, width_trials, cube



class Candidate(object):
    """ """
    def __init__(self, dm_trials, period_trials, width_trials, dpw_cube, subints=None, metadata=None):
        self.dm_trials = dm_trials
        self.period_trials = period_trials
        self.width_trials = width_trials
        self.dpw_cube = dpw_cube
        self.subints = subints
        self.metadata = metadata if metadata is not None else Metadata({})

    def __str__(self):
        name = type(self).__name__
        md = self.metadata
        return '{name:s} [P0 = {period:.9e}, W = {width:3d}, DM = {dm:.3f}, S/N = {snr:.2f}]'.format(
            name=name,
            period=md['best_period'],
            width=md['best_width'],
            dm=md['best_dm'],
            snr=md['best_snr']
            )

    def __repr__(self):
        return str(self)


    @classmethod
    def from_pipeline_output(cls, cluster, tseries, nbins=128, nsubs=64, logger=None):
        """
        Parameters:
        -----------
            cluster: list of Detection objects
                Cluster of all associated detections.
            tseries: TimeSeries
                TimeSeries from which originated the most signficant detection.
        """
        if not logger:
            logger = logging.getLogger('Candidate')

        # Detection with the highest S/N
        topdet = cluster.top_detection

        logmsg = "Creating Candidate from group of {:4d} detections. P0 = {:.9e}, DM = {:7.2f}, S/N = {:6.2f}".format(
            len(cluster),
            topdet.period,
            topdet.dm,
            topdet.snr
            )
        logger.info(logmsg)

        # Each detection contains a S/N versus Period and Width array
        # We need to resample those to a common set of period trials
        dm_trials, period_trials, width_trials, dpw_cube = create_dpw_cube(cluster, logger=logger)

        # NOTE: Must use deepcopy(), otherwise we just pass a reference and
        # candidates end up sharing the same Metadata object
        md = copy.deepcopy(topdet.metadata)
        md['best_period'] = topdet.period
        md['best_width'] = topdet.width
        md['best_dm'] = topdet.dm
        md['best_snr'] = topdet.snr

        ### Build subints
        # Check that nsubs and nbins do not exceed acceptable maximum
        nsubs = min(nsubs, int(md['tobs'] / topdet.period))
        nbins = min(nbins, int(topdet.period / md['tsamp']))
        logger.info("Creating SubIntegrations with nbins = {:d}, nsubs = {:d}".format(nbins, nsubs))
        try:
            subints = SubIntegrations.from_time_series(tseries, topdet.period, nbins=nbins, nsubs=nsubs)
        except Exception as ex:
            msg = "Failed to build SubIntegrations: {!s}".format(ex)
            logger.error(msg)
            subints = None
        return cls(dm_trials, period_trials, width_trials, dpw_cube, subints=subints, metadata=md)


    def save_hdf5(self, fname):
        with h5py.File(fname, 'w') as fobj:
            # Create a group to store metadata, as attribute of said group
            self.metadata._save_to_hdf5_file(fobj)

            # Save subints to their own group
            try:
                self.subints._save_to_hdf5_file(fobj)
            except:
                pass

            # Save DPW cube
            cube_group = fobj.create_group('dpw_cube')
            cube_group.create_dataset('dpw_cube', data=self.dpw_cube, dtype=np.float32)
            cube_group.create_dataset('dm_trials', data=self.dm_trials, dtype=np.float32)
            cube_group.create_dataset('period_trials', data=self.period_trials, dtype=np.float32)
            cube_group.create_dataset('width_trials', data=self.width_trials, dtype=np.float32)


    @classmethod
    def load_hdf5(cls, fname):
        with h5py.File(fname, 'r') as fobj:
            metadata = Metadata._from_hdf5_file(fobj)

            try:
                subints = SubIntegrations._from_hdf5_file(fobj)
            except:
                subints = None

            cube_group = fobj['dpw_cube']
            dm_trials = cube_group['dm_trials'].value
            period_trials = cube_group['period_trials'].value
            width_trials = cube_group['width_trials'].value
            dpw_cube = cube_group['dpw_cube'].value

        return cls(dm_trials, period_trials, width_trials, dpw_cube, subints=subints, metadata=metadata)
