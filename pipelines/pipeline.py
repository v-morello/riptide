##### Standard imports #####
import os
import glob
import logging
import operator

##### Non-standard imports #####
import numpy as np
import yaml

##### Local imports #####
from .. import TimeSeries, ffa_search, find_peaks
from ..clustering import cluster_1d
from .candidate import Candidate



def load_config(fname):
    """ Load YAML pipeline configuration file into a (nested) dictionary. """
    with open(fname, 'r') as fobj:
        config = yaml.safe_load(fobj)
    return config


def time_series_iterator(directory, pattern='*.inf'):
    """ Iterate through a set of PRESTO .inf files in given directory, whose
    file names match a given UNIX pattern. """
    logger = logging.getLogger('Pipeline')
    fnames = sorted(glob.glob(os.path.join(directory, pattern)))
    for fname in fnames:
        try:
            fpath = os.path.join(directory, fname)
            time_series = TimeSeries.from_presto_inf(fpath)
            yield fpath, time_series
        except:
            msg = "Failed to load file: {0:s}".format(fpath)
            logger.error(msg)
    raise StopIteration




class DetectionCluster(list):
    """ Cluster of Detection objects. It also keeps in memory the TimeSeries
    associated with the brightest detection, and only that one. """
    def __init__(self, detections):
        super(DetectionCluster, self).__init__(detections)
        ### Allow only the top detection to keep its 'time_series' attribute
        time_series = self.top_detection.time_series
        for det in self:
            if hasattr(det, 'time_series'):
                delattr(det, 'time_series')
        self.top_detection.time_series = time_series

    @property
    def top_detection(self):
        return max(self, key=operator.attrgetter('snr'))

    def __str__(self):
        name = type(self).__name__
        return "{name:s} [size = {size:4d}, P0 = {top.period:.9e}, DM = {top.dm:8.3f}, S/N = {top.snr:6.2f}]".format(
            top=self.top_detection,
            name=name,
            size=len(self)
            )

    def __repr__(self):
        return str(self)


class Pipeline(object):
    """ """
    def __init__(self, config_fname, input_dir, output_dir, pattern='*.inf'):
        self.config = load_config(config_fname)
        self.input_dir = os.path.realpath(input_dir)
        self.output_dir = os.path.realpath(output_dir)
        self.pattern = pattern

        # Configure logger
        self.logger = logging.getLogger('Pipeline')
        self.logger.setLevel(logging.DEBUG)

        # Without this check, creating a new Pipeline instance would add an
        # extra handler and we would get duplicate log messages
        if not self.logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter(fmt='%(asctime)s.%(msecs)03d: %(message)s',datefmt='%Y-%m-%d,%H:%M:%S')
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)


    def get_iterator(self):
        return time_series_iterator(self.input_dir, pattern=self.pattern)

    def _run_search(self):
        self.detections = []
        self.clusters = []

        # Track tobs of every time series, this will be useful later for
        # clustering
        tobs = []

        for fname, tseries in self.get_iterator():
            tobs.append(tseries.length)
            dm = tseries.metadata['dm']

            # Compute periodogram
            self.logger.debug("Computing Periodogram @ DM = {0:.3f}".format(dm))
            ts, plan, pgram = ffa_search(tseries, **self.config['search'])

            # Identify detections. Give them a temporary attribute time_series
            # Only detections that are the brightest of their cluster will
            # be allowed to keep it in memory
            self.logger.debug("Finding Peaks @ DM = {0:.3f}".format(dm))
            dets = find_peaks(pgram, **self.config['peak_detection'])
            for det in dets:
                # NOTE: we keep the de-redenned and normalised copy
                # of the time series, NOT the original one.
                det.time_series = ts

            self.detections = self.detections + dets
            self.logger.debug("New Detections: {0:d}".format(len(dets)))
            self.logger.debug("Total Detections: {0:d}".format(len(self.detections)))
            self._cluster_detections(tobs=np.median(tobs))
            nts = self.num_memorised_time_series()
            self.logger.debug("TimeSeries in RAM: {0:d}".format(nts))

    def _cluster_detections(self, tobs):
        self.logger.debug("Clustering Detections ...")
        periods = np.asarray([det.period for det in self.detections])
        dbi = tobs / periods

        cluster_indices = cluster_1d(
            dbi,
            self.config['peak_detection']['peak_clustering_radius']
            )

        self.clusters = [
            DetectionCluster([self.detections[ii] for ii in indices])
            for indices in cluster_indices
            ]
        self.logger.debug("Total Clusters: {0:d}".format(len(self.clusters)))

    def _build_candidates(self):
        nbins = self.config['candidates']['nbins']
        nsubs = self.config['candidates']['nsubs']
        self.candidates = [
            Candidate.from_pipeline_output(cluster, nbins=nbins, nsubs=nsubs)
            for cluster in self.clusters
            ]

    def num_memorised_time_series(self):
        return len(set(id(det.time_series) for det in self.detections if hasattr(det, 'time_series')))

    def _remove_harmonics(self):
        pass

    def run(self):
        self._run_search()
        self._remove_harmonics()
        if self.clusters:
            self._build_candidates()
