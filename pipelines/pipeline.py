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
from riptide.reading import PrestoInf, SigprocHeader
from riptide.pipelines.harmonic_filtering import test_harmonic_relationship

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
    during processing. Once processing is over, build clusters out of the
    accumulated detections. The main PipelineManager object is responsible
    for removing harmonics and saving candidates.

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

    def cluster_detections(self):
        self.logger.info("Clustering Detections ...")
        periods = np.asarray([det.period for det in self.detections])
        tobs = np.median([det.metadata['tobs'] for det in self.detections])
        dbi = tobs / periods

        clrad = self.config['detect']['peak_clustering_radius']
        cluster_indices = cluster_1d(dbi, clrad)

        self.clusters = [
            DetectionCluster([self.detections[ii] for ii in indices])
            for indices in cluster_indices
            ]
        self.logger.info("Clustering complete. Total Clusters: {0:d}".format(len(self.clusters)))



class PipelineManager(object):
    """ Responsible for the outermost DM loop and top-level pulsar search
    management. """
    def __init__(self, config_path, override_keys={}):
        """
        Parameters:
        -----------
            config_path: str
                Path to the YAML config file of the PipelineManager
            override_keys: dict
                Updated values for some keys of the YAML config file,
                specified through the command line arguments of this script.
        """
        self.config_path = os.path.realpath(config_path)
        self.config = parse_yaml_config(self.config_path)
        self.config.update(override_keys)

        self.clusters = []
        self.candidates = []

        self.configure_logger()
        self.configure_loaders()
        self.configure_searches()

    def configure_logger(self):
        logger_name = 'PipelineManager'
        self.logger = get_logger(logger_name)

    def configure_searches(self):
        self.searches = []
        config_dir, config_name = os.path.split(self.config_path)
        for search_config_fname in self.config['search_configs']:
            search_config_path = os.path.join(config_dir, search_config_fname)
            search = PulsarSearch(self, parse_yaml_config(search_config_path))
            self.searches.append(search)
            self.logger.info("Configured PulsarSearch '{:s}'".format(search.name))
        self.logger.info("Configured a total of {:d} searches.".format(len(self.searches)))

    def configure_loaders(self):
        fmt = self.config['data_format'].strip().lower()
        if fmt == 'presto':
            self.loader = TimeSeries.from_presto_inf
            self.dm_getter = lambda fname: PrestoInf(fname).dm
        elif fmt == 'sigproc':
            self.loader = TimeSeries.from_sigproc
            self.dm_getter = lambda fname: SigprocHeader(fname)['refdm']
        else:
            raise ValueError("Invalid data format '{s}'".format(fmt))
        self.logger.info("Specified file format: {:s}".format(fmt))

    def select_dm_trials(self):
        """ Build a list of files to process """
        glob_pattern = self.config['glob']
        filenames = sorted(glob.glob(glob_pattern))
        self.logger.info("Found a total of {:d} file names corresponding to specified pattern \"{:s}\"".format(len(filenames), glob_pattern))
        self.logger.info("Retching DM trial values from headers. This may take a while ...")
        dm_trials = {
            self.dm_getter(fname): fname
            for fname in filenames
            }
        self.logger.info("DM trial values have been read.")

        # Helper iterator, used to select a sequence of DM trials according to
        # config parameters
        def iter_steps(sequence, vmin, vmax, step):
            # Ignore steps outside of bounds
            array = np.asarray(sorted(list(sequence)))
            mask = (array >= vmin) & (array <= vmax)
            array = array[mask]

            # Yield values separated by at least 'step'
            last = None
            rtol = 1e-7 # Deal with float rounding errors
            for value in array:
                if last is None or value - last >= step * (1-rtol):
                    last = value
                    yield value

        # Set max DM trial as a function of both the hard maximum limit and
        # dmsinb_max
        dm_min = self.config['dm_min']
        dm_max = self.config['dm_max']
        dm_step = self.config['dm_step']
        dmsinb_max = self.config['dmsinb_max']
        if filenames:
            # Get coordinates of the first file in the list
            tseries = self.loader(filenames[0])
            skycoord = tseries.metadata['skycoord']
            glat_radians = skycoord.galactic.b.rad
            self.logger.info("Read galactic latitude from \"{:s}\" b = {:.3f} deg".format(filenames[0], skycoord.galactic.b.deg))

            galactic_dm_max = dmsinb_max / (np.sin(abs(glat_radians)) + 1e-4)
            self.logger.info("Requested maximum value of DM x sin |b| ({:.3f}) corresponds to DM = {:.3f}".format(dmsinb_max, galactic_dm_max))

            # Update value of dm_max
            dm_max = min(dm_max, galactic_dm_max)

        self.logger.info("Selecting DM trials in range [{:.3f}, {:.3f}] with a minimum step of {:.3f}".format(dm_min, dm_max, dm_step))

        # NOTE: this is an iterator
        dm_trial_values = iter_steps(dm_trials.keys(), dm_min, dm_max, dm_step)
        self.dm_trial_paths = [dm_trials[value] for value in dm_trial_values]
        self.logger.info("Selected {:d} DM trials to process".format(len(self.dm_trial_paths)))

    def iter_batches(self):
        """ Iterate through input time series in batches. Yields a list of
        num_processes TimeSeries at each iteration. """
        num_processes = self.config['num_processes']
        paths = self.dm_trial_paths

        num_dm_trials = len(paths)
        self.logger.info("Preparing to iterate through DM trials. Number of input files: {:d}". format(num_dm_trials))

        for batch in grouper(paths, num_processes):
            tsbatch = list(map(self.loader, batch))
            yield tsbatch

    def fetch_clusters(self):
        """ Place all DetectionCluster objects from all the searches into a
        single list. Give each DetectionCluster a new attribute tracking
        which PulsarSearch it belongs to."""
        self.clusters = []
        for search in self.searches:
            for cl in search.clusters:
                cl.search = search
                self.clusters.append(cl)

    def remove_harmonics(self):
        self.logger.info("Removing harmonics ...")

        max_denominator = self.config['harmonic_filtering']['max_denominator']
        snr_tol = self.config['harmonic_filtering']['snr_tol']
        dft_bins_tol = self.config['harmonic_filtering']['dft_bins_tol']

        self.clusters = sorted(self.clusters, key=lambda cl: cl.top_detection.snr, reverse=True)
        for cl in self.clusters:
            cl.is_harmonic = False
        num_harmonics_flagged = 0
        for fd, hm in itertools.combinations(self.clusters, 2):
            # If fundamental was previously flagged, move on.
            if fd.is_harmonic:
                continue
            is_harmonic, fraction = test_harmonic_relationship(
                fd, hm,
                max_denominator=max_denominator,
                snr_tol=snr_tol,
                dft_bins_tol=dft_bins_tol)
            if is_harmonic:
                hm.is_harmonic = True
                num_harmonics_flagged += 1
                self.logger.info("{!s} was flagged as a harmonic of {!s} with ratio {:d}/{:d}".format(hm, fd, fraction.numerator, fraction.denominator))
        self.logger.info("Flagged {:d} harmonics".format(num_harmonics_flagged))
        self.clusters = [cl for cl in self.clusters if not cl.is_harmonic]
        self.logger.info("Retained {:d} final Candidates".format(len(self.clusters)))


    def _apply_candidate_filter(self, filter_name, func):
        num_clusters = len(self.clusters)
        valid_clusters = list(filter(func, self.clusters))
        num_invalid = num_clusters - len(valid_clusters)
        self.logger.info("Applied candidate filter \"{:s}\" on {:d} DetectionClusters: {:d} were removed".format(filter_name, num_clusters, num_invalid))
        self.clusters = valid_clusters

    def apply_candidate_filters(self):
        params = self.config['candidate_filters']
        dm_min = params['dm_min']
        snr_min = params['snr_min']
        max_number = params['max_number']

        if dm_min:
            self._apply_candidate_filter(
                "DM >= {:.2f}".format(dm_min),
                lambda cl: cl.top_detection.dm >= dm_min)

        if snr_min:
            self._apply_candidate_filter(
                "S/N >= {:.2f}".format(snr_min),
                lambda cl: cl.top_detection.snr >= snr_min)

        if max_number:
            self.logger.info("Keeping only the top {:d} brightest candidates".format(max_number))
            self.clusters = self.clusters[:max_number]


    def build_candidates(self):
        self.logger.info("Building Candidates ...")
        self.candidates = []

        for cluster in self.clusters:
            search = cluster.search

            # Get the original parameters of the PulsarSearch that Found
            # this cluster
            rmed_width = search.config['search']['rmed_width']
            rmed_minpts = search.config['search']['rmed_minpts']
            nbins = search.config['candidates']['nbins']
            nsubs = search.config['candidates']['nsubs']

            # Re-load TimeSeries associated to the top detection, and run
            # the same pre-processing again.
            fname = cluster.top_detection.metadata['fname']

            try:
                tseries = self.loader(fname)
                tseries.deredden(rmed_width, minpts=rmed_minpts, inplace=True)
                tseries.normalise(inplace=True)

                candidate = Candidate.from_pipeline_output(cluster, tseries, nbins=nbins, nsubs=nsubs, logger=self.logger)
                self.candidates.append(candidate)
            except Exception as error:
                self.logger.error("ERROR: Failed to build candidate from {!s}. Reason: {!s}".format(cluster, error))


        self.candidates = sorted(self.candidates, key=lambda cd: cd.metadata['best_snr'], reverse=True)
        self.logger.info("Done building candidates.")

    def save_candidates(self):
        outdir = self.config['outdir']
        self.logger.info("Saving {:d} candidates to output directory: {:s}".format(len(self.candidates), outdir))
        for index, cand in enumerate(self.candidates, start=1):
            outpath = os.path.join(outdir, "riptide_cand_{:04d}.h5".format(index))
            self.logger.info("Saving {!s} to file {:s}".format(cand, outpath))
            cand.save_hdf5(outpath)

    def run(self):
        self.logger.info("Starting pipeline ...")
        self.logger.info("Selecting DM trials ...")
        self.select_dm_trials()

        for tsbatch in self.iter_batches():
            dms = sorted([ts.metadata['dm'] for ts in tsbatch])
            self.logger.info("Processing DM trials: {!s}".format(dms))
            for search in self.searches:
                search.process_time_series_batch(tsbatch)

        self.logger.info("All DM trials have been processed. Clustering detections ...")
        for search in self.searches:
            search.cluster_detections()

        self.fetch_clusters()
        # NOTE: As of 22 Jan 2018 I am turning this off for safety
        # Tests on LOTAAS beams have shown that in the presence
        # of strong RFI, pulsars can be removed.
        #self.remove_harmonics()
        self.apply_candidate_filters()
        self.build_candidates()
        self.save_candidates()

        self.logger.info("Searches closed. Pipeline run complete.")


################################################################################

def parse_arguments():
    """ Parse command line arguments with which the script was called. Returns
    an object containing them all.
    """
    ### IMPORTANT NOTE: The 'config' argument specifies the path to the main
    ### YAML configuration file. All the other arguments override keys
    ### already present in the YAML config file.
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'config', type=str,
        help="YAML configuration file for the PipelineManager"
        )
    parser.add_argument(
        'glob', type=str,
        help="glob pattern matching the input files to be processed. IMPORTANT: Must be delimited by double quotes."
        )
    parser.add_argument(
        'outdir', type=str,
        help="Output directory for data products"
        )
    args = parser.parse_args()
    return args



def main():
    args = parse_arguments()
    # Get absolute paths right away, just to be safe
    args.config = os.path.realpath(args.config)
    args.outdir = os.path.realpath(args.outdir)
    manager = PipelineManager(args.config, override_keys=vars(args))
    manager.run()


if __name__ == '__main__':
    main()
