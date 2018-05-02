##### Standard imports #####
import logging
import operator
import copy

##### Non-standard imports #####
import numpy as np
import h5py

##### Local imports #####
from .. import Metadata, SubIntegrations



class ResponseCurve(object):
    """ Stores a tuple of arrays representing e.g. a S/N versus DM curve,
    or any other trial parameter. """
    def __init__(self, trials, snr):
        self.trials = trials
        self.snr = snr


class Candidate(object):
    """ The final data product of a search. """
    def __init__(self, dm_curve, period_curve, width_curve, subints=None, metadata=None):
        self.dm_curve = dm_curve
        self.period_curve = period_curve
        self.width_curve = width_curve
        self.subints = subints
        self.metadata = metadata if metadata is not None else Metadata({})

    def __str__(self):
        name = type(self).__name__
        md = self.metadata
        return '{name:s} [P0 = {period:.9e}, W = {width:3d}, DM = {dm:7.3f}, S/N = {snr:6.2f}]'.format(
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
            cluster: DetectionCluster
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

        ### New simpler way to build DM and period curves
        # Just extract the period and width curves at best DM
        # topdet.snr_trials has shape (num_period_trials, num_width_trials)
        iperiod = abs(topdet.period_trials - topdet.period).argmin()
        iwidth = abs(topdet.width_trials - topdet.width).argmin()
        period_curve = ResponseCurve(topdet.period_trials, topdet.snr_trials[:, iwidth])
        width_curve = ResponseCurve(topdet.width_trials, topdet.snr_trials[iperiod, :])

        # DM curve: for each DM trial, extract the best S/N at best width only
        dm_trials = np.asarray([det.dm for det in cluster])
        snr_trials = np.asarray([det.snr_trials[:, iwidth].max() for det in cluster])
        order = dm_trials.argsort()
        dm_curve = ResponseCurve(dm_trials[order], snr_trials[order])

        # NOTE: Must use deepcopy(), otherwise we just pass a reference and
        # candidates end up sharing the same Metadata object
        md = copy.deepcopy(topdet.metadata)
        md['best_period'] = topdet.period
        md['best_width'] = topdet.width
        md['best_dm'] = topdet.dm
        md['best_snr'] = topdet.snr
        md['best_ducy'] = topdet.ducy

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
        return cls(dm_curve, period_curve, width_curve, subints=subints, metadata=md)


    def save_hdf5(self, fname):
        with h5py.File(fname, 'w') as fobj:
            # Create a group to store metadata, as attribute of said group
            self.metadata._save_to_hdf5_file(fobj)

            # Save subints to their own group
            try:
                self.subints._save_to_hdf5_file(fobj)
            except:
                pass

            # Save response curves
            cube_group = fobj.create_group('response_curves')
            cube_group.create_dataset('dm_curve_trials', data=self.dm_curve.trials, dtype=np.float32)
            cube_group.create_dataset('period_curve_trials', data=self.period_curve.trials, dtype=np.float32)
            cube_group.create_dataset('width_curve_trials', data=self.width_curve.trials, dtype=np.float32)
            cube_group.create_dataset('dm_curve_snr', data=self.dm_curve.snr, dtype=np.float32)
            cube_group.create_dataset('period_curve_snr', data=self.period_curve.snr, dtype=np.float32)
            cube_group.create_dataset('width_curve_snr', data=self.width_curve.snr, dtype=np.float32)



    @classmethod
    def load_hdf5(cls, fname):
        with h5py.File(fname, 'r') as fobj:
            metadata = Metadata._from_hdf5_file(fobj)

            try:
                subints = SubIntegrations._from_hdf5_file(fobj)
            except:
                subints = None

            curves_group = fobj['response_curves']

            dm_curve = ResponseCurve(
                curves_group['dm_curve_trials'].value,
                curves_group['dm_curve_snr'].value
                )

            period_curve = ResponseCurve(
                curves_group['period_curve_trials'].value,
                curves_group['period_curve_snr'].value
                )

            width_curve = ResponseCurve(
                curves_group['width_curve_trials'].value,
                curves_group['width_curve_snr'].value
                )

        return cls(dm_curve, period_curve, width_curve, subints=subints, metadata=metadata)
