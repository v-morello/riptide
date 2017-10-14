### Standard imports ###
import os
import glob
import logging
import argparse
import itertools
import operator
from multiprocessing import Pool

### Non-standard imports ###
import yaml
import numpy as np

### Local imports ###
from riptide import TimeSeries, ffa_search, find_peaks
from riptide.pipelines import Candidate
from riptide.clustering import cluster_1d


def parse_yaml_config(fname):
    with open(fname, 'r') as fobj:
        config = yaml.load(fobj)
    return config

def get_logger(name, level=logging.INFO):
    logger = logging.getLogger(name)
    logger.setLevel(level)

    formatter = logging.Formatter(
        fmt="%(asctime)s.%(msecs)03d - %(name)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
        )
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger

def grouper(iterable, n):
    """ Iterate through iterable, yielding groups of n elements. The last
    group may have less than n elements. """
    # we want izip_longest() in python 2.x
    # or zip_longest() in python 3.x
    if hasattr(itertools, 'zip_longest'):
        zipper = itertools.zip_longest
    else:
        zipper = itertools.izip_longest
    args = [iter(iterable)] * n
    for group in zipper(*args, fillvalue=None):
        filtered_group = [val for val in group if val is not None]
        yield filtered_group


class DetectionCluster(list):
    """ Cluster of Detection objects. It also keeps in memory the TimeSeries
    associated with the brightest detection, and only that one. """
    def __init__(self, detections):
        super(DetectionCluster, self).__init__(detections)

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



class PulsarSearchWorker(object):
    """ Function-like object that takes a single TimeSeries as an argument,
    and outputs a list of detections. This is to circumvent a limitation of
    multiprocessing.Pool.map() which requires the mapped function to:
    - take just one argument.
    - be pickle-able, and therefore be defined at the top level of a module
    """
    def __init__(self, config):
        """
        Parameters:
            config : dict
                Configuration parameters of the parent PulsarSearch object.
        """
        self.config = config

    def __call__(self, tseries):
        ts, plan, pgram = ffa_search(tseries, **self.config['search'])
        dets = find_peaks(pgram, **self.config['detect'])
        return dets



class PulsarSearch(object):
    """ Gets fed time series and accumulates detections at various DM trials
    during processing. Once processing is over, build candidates out of the
    accumulated detections.

    The idea is to create one PulsarSearch object for different non-overlapping
    period search ranges, where we search longer periods with more phase bins.
    """
    def __init__(self, manager, config):
        """
        Parameters:
            manager: PipelineManager
                Parent manager.
            config : dict
                Configuration parameters of the search. It is read from a
                YAML config file, and updated by extra parameters from command
                line arguments.
        """
        self.manager = manager
        self.config = config
        self.configure_logger()
        self.detections = []
        self.clusters = []
        self.candidates = []

    @property
    def name(self):
        return self.config['name']

    @property
    def num_processes(self):
        return self.manager.config['num_processes']

    def configure_logger(self):
        logger_name = '.'.join(['PulsarSearch', self.name])
        self.logger = get_logger(logger_name)

    def process_time_series_batch(self, batch):
        """ Processes several time series in parallel using the muliprocessing
        module.

        Parameters:
        -----------
            batch: list
                list of TimeSeries to process.
        """
        pool = Pool(processes=self.num_processes)

        # Split the TimeSeries to search between processes
        # 'output' is a list of lists of Detections
        self.logger.info("Searching batch of {:d} TimeSeries using {:d} worker processes ...".format(len(batch), self.num_processes))
        output = pool.map(PulsarSearchWorker(self.config), batch)
        pool.close()
        pool.join()

        # Accumulate the new Detections
        new_detections = [det for sublist in output for det in sublist]
        self.logger.info("Search complete. New detections: {:d}".format(len(new_detections)))
        self.detections = self.detections + new_detections
        self.logger.info("Total detections stored: {:d}".format(len(self.detections)))

    def close(self):
        """ Close the search: turn the list of accumulated detections into
        Candidates and save them to disk. """
        self.cluster_detections()
        self.remove_harmonics()
        self.build_candidates()
        self.save_candidates()

    def cluster_detections(self):
        self.logger.debug("Clustering Detections ...")
        periods = np.asarray([det.period for det in self.detections])
        tobs = np.median([det.metadata['tobs'] for det in self.detections])
        dbi = tobs / periods

        clrad = self.config['detect']['peak_clustering_radius']
        cluster_indices = cluster_1d(dbi, clrad)

        self.clusters = [
            DetectionCluster([self.detections[ii] for ii in indices])
            for indices in cluster_indices
            ]
        self.logger.debug("Clustering complete. Total Clusters: {0:d}".format(len(self.clusters)))


    def remove_harmonics(self):
        pass

    def build_candidates(self):
        self.logger.debug("Building Candidates ...")
        self.candidates = []
        rmed_width = self.config['search']['rmed_width']
        rmed_minpts = self.config['search']['rmed_minpts']
        nbins = self.config['candidates']['nbins']
        nsubs = self.config['candidates']['nsubs']

        for cluster in self.clusters:
            # Re-load TimeSeries associated to the top detection, and run
            # the same pre-processing again.
            fname = cluster.top_detection.metadata['fname']
            tseries = self.manager.loader(fname)
            tseries.deredden(rmed_width, minpts=rmed_minpts, inplace=True)
            tseries.normalise(inplace=True)

            candidate = Candidate.from_pipeline_output(cluster, tseries, nbins=nbins, nsubs=nsubs, logger=self.logger)
            self.candidates.append(candidate)
        self.logger.debug("Done building candidates.")

    def save_candidates(self):
        pass



class PipelineManager(object):
    """ Responsible for the outermost DM loop and top-level pulsar search
    management. """
    def __init__(self, config):
        self.config = config
        self.configure_logger()
        self.configure_loader()
        self.configure_searches()

    def configure_logger(self):
        logger_name = 'PipelineManager'
        self.logger = get_logger(logger_name)

    def configure_searches(self):
        self.searches = []
        config_dir, config_name = os.path.split(self.config['config'])
        for search_config_fname in self.config['search_configs']:
            search_config_path = os.path.join(config_dir, search_config_fname)
            search = PulsarSearch(self, parse_yaml_config(search_config_path))
            self.searches.append(search)
            self.logger.info("Configured PulsarSearch '{:s}'".format(search.name))
        self.logger.info("Configured a total of {:d} searches.".format(len(self.searches)))

    def configure_loader(self):
        fmt = self.config['data_format'].strip().lower()
        if fmt == 'presto':
            self.loader = TimeSeries.from_presto_inf
            self.logger.info("Specified file format: {:s}".format(fmt))
        else:
            raise ValueError("Invalid data format '{s}'".format(fmt))

    def iter_batches(self):
        """ Iterate through input time series in batches. Yields a list of
        num_processes TimeSeries at each iteration. """
        num_processes = self.config['num_processes']
        paths = sorted(glob.glob(self.config['input_pattern']))
        num_dm_trials = len(paths)
        self.logger.info("Preparing to iterate through DM trials. Number of input files: {:d}". format(num_dm_trials))

        for batch in grouper(paths, num_processes):
            tsbatch = list(map(self.loader, batch))
            yield tsbatch

    def run(self):
        self.logger.info("Starting pipeline ...")
        for tsbatch in self.iter_batches():
            dms = sorted([ts.metadata['dm'] for ts in tsbatch])
            self.logger.info("Processing DM trials: {!s}".format(dms))
            for search in self.searches:
                search.process_time_series_batch(tsbatch)
        self.logger.info("All DM trials have been processed. Closing searches ...")

        for search in self.searches:
            search.close()
        self.logger.info("Searches closed. Pipeline run complete.")


################################################################################

def parse_arguments():
    """ Parse command line arguments with which the script was called. Returns
    an object containing them all.
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--config', type=str, required=True,
        help="YAML configuration file for the PipelineManager"
        )
    parser.add_argument(
        '--input_pattern', type=str, required=True,
        help="UNIX shell pattern matching the input files to be processed. IMPORTANT: Must be delimited by double quotes."
        )
    parser.add_argument(
        '--outdir', type=str, required=True,
        help="Output directory for data products"
        )
    args = parser.parse_args()
    return args



def main():
    args = parse_arguments()

    # IMPORTANT: get the full path of the manager config file. We need to know
    # in which directory it lies, because that's where we also look for the
    # config files for individual PulsarSearches.
    args.config = os.path.realpath(args.config)

    config = parse_yaml_config(args.config)
    config.update(vars(args))

    manager = PipelineManager(config)
    manager.run()


if __name__ == '__main__':
    main()
