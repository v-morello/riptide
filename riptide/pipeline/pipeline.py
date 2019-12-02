import logging
import argparse
import os
import itertools
import multiprocessing
import traceback
from collections import defaultdict

import yaml
import json
import numpy as np
import pandas

import threadpoolctl

from riptide import __version__, ffa_search, find_peaks, Candidate
from riptide.clustering import cluster1d

from riptide.pipeline.dmiter import DMIterator
from riptide.pipeline.worker_pool import WorkerPool
from riptide.pipeline.peak_cluster import PeakCluster, clusters_to_dataframe
from riptide.pipeline.harmonic_testing import htest

from riptide.serialization import save_json
from riptide.timing import timing


log = logging.getLogger('riptide.pipeline')


class CandidateWriter(object):
    """
    func-like object to be used in conjunction with multiprocessing.Pool
    to write candidates with multiple processes
    """
    def __init__(self, outdir, plot=False):
        self.outdir = os.path.realpath(outdir)
        self.plot = plot

    def __call__(self, arg):
        """
        arg is a tuple (rank, candidate) with types
        (int, Candidate)
        """
        rank, cand = arg
        fname = os.path.join(self.outdir, f'candidate_{rank:04d}.json')
        log.debug(f"Saving to {fname}: {cand}")
        save_json(fname, cand)
        if self.plot:
            fname = os.path.join(self.outdir, f'candidate_{rank:04d}.png')
            log.debug(f"Saving plot to {fname}")
            cand.savefig(fname)


class Pipeline(object):
    """
    Top-level class that runs a multiple DM trial search.

    Parameters
    ----------
    conf: dict
        Configuration dictionary loaded from YAML file
    """
    def __init__(self, config):
        # TODO: validate config
        self.config = config
        self.dmiter = None
        self.worker_pool = None
        self.peaks = []
        self.clusters = []
        self.candidates = []

    def wmin(self):
        """ Minimum pulse width being searched for """
        search_ranges = self.config['ranges']
        min_widths = [
            kw['ffa_search']['period_min'] / kw['ffa_search']['bins_min']
            for kw in search_ranges
        ]
        return min(min_widths)

    def get_search_range(self, period):
        """
        Get the search range parameters (from configuration file) that the
        given candidate period falls into. This is used at Candidate building
        stage to retrieve how many phase bins and subints should be used
        when folding the data at this particular period.

        Parameters
        ----------
        period: float
            Candidate period in seconds

        Returns
        -------
        range: dict
            Parameters of the search range the given periods falls into
        """
        # TODO: The range parameters need to be validated on initialization
        # We need to be sure that they form a partition of a wider range
        # ie. do not overlap and leave no gaps in between
        # The code below can return wrong results if the ranges do not connect
        # perfectly with each other
        ranges = sorted(
            self.config['ranges'], 
            key=lambda r: r['ffa_search']['period_max']
            )
        
        pmin_global = min(rng['ffa_search']['period_min'] for rng in ranges)
        pmax_global = max(rng['ffa_search']['period_max'] for rng in ranges) 

        if period < pmin_global:
            msg = (
                f"Given period={period:.9f} is shorter than the minimum search period={pmin_global:.9f}."
                " This will not affect the processing but it should NOT be happening."
            )
            log.warning(msg)
            return dict(ranges[0])

        # This can rightfully happen on occasion
        # We actually search slightly higher than pmax_global in practice
        if period >= pmax_global:
            return dict(ranges[-1])

        for rng in ranges:
            pmin = rng['ffa_search']['period_min']
            pmax = rng['ffa_search']['period_max']
            if pmin <= period < pmax:
                return dict(rng)

    @timing
    def prepare(self, files):
        """
        Inspect input files and select a minimal set to process

        files: list
        """
        log.info("Preparing pipeline")
        conf = self.config

        # The DM iterator is in charge of:
        # - checking the time series
        # - selecting the DM trials in the right range with the right dynamic step
        # - yielding them in chunks of size = number of parallel processes
        self.dmiter = DMIterator(
            dm_min=conf['dmselect']['min'],
            dm_max=conf['dmselect']['max'],
            dmsinb_max=conf['dmselect']['dmsinb_max'],
            wmin=self.wmin(),
            fmin=conf['data']['fmin'],
            fmax=conf['data']['fmax'],
            nchans=conf['data']['nchans'],
        )
        self.worker_pool = WorkerPool(conf['dereddening'], conf['ranges'], processes=conf['processes'])

        log.debug("Input files: {}".format(len(files)))
        self.dmiter.prepare(files, fmt=self.config['data']['format'])
        log.info("Pipeline ready")
        
    @timing
    def search(self):
        """
        Search all selected files
        """
        log.info("Running search")
        peaks = []
        for chunk in self.dmiter.iterate(chunksize=self.config['processes']):
            peaks.extend(
                self.worker_pool.process_chunk(chunk)
            )        
        self.peaks = sorted(peaks, key=lambda p: p.period)
        log.info("Total peaks found: {}".format(len(peaks)))
        log.info("Search complete")

    @timing
    def cluster_peaks(self):
        if not self.peaks:
            log.info("No peaks found: skipping clustering")
            return

        log.info("Clustering peaks")
        conf = self.config
        tmed = self.dmiter.tobs_median()
        clrad = conf['clustering']['radius'] / tmed

        log.debug(f"Median Tobs = {tmed:.2f} s")
        log.debug(f"Frequency clustering radius = {clrad:.3e} Hz")

        # NOTE: peaks is assumed to be sorted by increasing period
        freqs = np.asarray([p.freq for p in self.peaks])
        cluster_ids = cluster1d(freqs, clrad, already_sorted=True)

        self.clusters = [
            PeakCluster((self.peaks[ii] for ii in ids))
            for ids in cluster_ids
        ]

        log.info(f"Total clusters found: {len(self.clusters)}")
        log.info("Clustering done")

    @timing
    def flag_harmonics(self):
        if not self.clusters:
            log.info("No clusters found: skipping harmonic flagging")
            return

        log.info("Flagging harmonics")
        tobs = self.dmiter.tobs_median()
        fmin = self.config['data']['fmin']
        fmax = self.config['data']['fmax']
        kwargs = self.config['harmonic_flagging']
        
        clusters_decreasing_snr = sorted(self.clusters, key=lambda c: c.centre.snr, reverse=True)

        # Assign ranks first
        for rank, cl in enumerate(clusters_decreasing_snr):
            cl.rank = rank

        # Flag harmonics
        for F, H in itertools.combinations(clusters_decreasing_snr, 2):
            if F.is_harmonic or H.is_harmonic:
                continue

            related, fraction = htest(F.centre, H.centre, tobs, fmin, fmax, **kwargs)

            if related:
                H.parent_fundamental = F
                H.hfrac = fraction
        
        harmonics = list(filter(lambda c: c.is_harmonic, self.clusters))
        log.info(f"Harmonics flagged: {len(harmonics)}")
        log.info(f"Fundamental clusters: {len(self.clusters) - len(harmonics)}")
        log.info("Harmonic flagging done")

    @timing
    def apply_candidate_filters(self):
        log.info("Applying candidate filters")
        params = self.config['candidate_filters']

        clusters_trimmed = self.clusters

        # DM cut
        dm_min = params['dm_min']
        if dm_min is not None:
            log.warning(f"Applying DM threshold of {dm_min}")
            clusters_trimmed = list(filter(lambda c: c.centre.dm >= dm_min, clusters_trimmed))

        # S/N cut
        snr_min = params['snr_min']
        if snr_min is not None:
            log.warning(f"Applying S/N threshold of {snr_min}")
            clusters_trimmed = list(filter(lambda c: c.centre.snr >= snr_min, clusters_trimmed))

        # Harmonic removal
        if params['remove_harmonics']:
            log.warning("Harmonic removal is enabled, clusters flagged as harmonics will NOT be output as candidates")
            clusters_trimmed = list(filter(lambda c: not c.is_harmonic, clusters_trimmed))

        # Cap on number of candidates to build
        nmax = params['max_number']
        if nmax:
            if len(clusters_trimmed) > nmax:
                nleft = len(clusters_trimmed)
                nexcess = nleft - nmax
                log.warning(
                    f"Number of clusters remaining ({nleft}) exceeds the maximum specified number of candidates ({nmax}). "
                    f"The faintest {nexcess} will not be saved as candidates"
                )
            clusters_trimmed = sorted(clusters_trimmed, key=lambda c: c.centre.snr, reverse=True)
            clusters_trimmed = clusters_trimmed[:nmax]

        self.clusters = clusters_trimmed
        log.info(f"Candidate filters applied. Clusters remaining: {len(self.clusters)}")


    @timing
    def build_candidates(self):
        log.info("Building candidates")
        #remove_harmonics = self.config['remove_harmonics']         
        clusters_decreasing_snr = sorted(self.clusters, key=lambda c: c.centre.snr, reverse=True)

        if not clusters_decreasing_snr:
            log.info("No clusters: no candidates to build")
            return

        # if remove_harmonics:
        #     log.warning("Harmonic removal is enabled, clusters flagged as harmonics will NOT be output as candidates")
        #     clusters_decreasing_snr = list(filter(lambda c: not c.is_harmonic, clusters_decreasing_snr))

        # Group candidates by DM, so that we can avoid re-loading some TimeSeries
        # multiple times
        grouped_clusters = defaultdict(list)
        for cl in clusters_decreasing_snr:
            grouped_clusters[cl.centre.dm].append(cl)

        log.debug(f"{len(clusters_decreasing_snr)} candidates to build from {len(grouped_clusters)} TimeSeries")

        for dm, clusters in grouped_clusters.items():
            ts = self.dmiter.get_dm_trial(dm)
            ts = ts.deredden(
                width=self.config['dereddening']['rmed_width'], 
                minpts=self.config['dereddening']['rmed_minpts'])
            ts = ts.normalise()

            for cl in clusters:
                try:
                    rng = self.get_search_range(cl.centre.period)
                    cand = Candidate.from_pipeline_output(
                        ts, cl, rng['candidates']['bins'], subints=rng['candidates']['subints'])
                    self.candidates.append(cand)

                # NOTE: this should never happen ideally
                # But just in case, avoid losing all results due to one candidate
                # failing to be built
                except Exception as err:
                    log.error(err)
                    log.error(traceback.format_exc())

        self.candidates = sorted(self.candidates, key=lambda c: c.params['snr'], reverse=True)
        log.info(f"Total candidates: {len(self.candidates)}")
        log.info("Done building candidates")

    @timing
    def save_products(self, outdir=os.getcwd()):
        """
        """
        log.info("Building products")

        if not self.peaks:
            log.info("No peaks found: no data products to save")
            return

        # CSV of Peak data
        df_peaks = pandas.DataFrame.from_dict([
            peak.summary_dict() for peak in self.peaks
        ])
        df_peaks_fname = os.path.join(outdir, 'peaks.csv')
        df_peaks.to_csv(df_peaks_fname, sep=',', index=False, float_format='%.9f')
        log.info("Saved Peak data to {!r}".format(df_peaks_fname))

        ### CSV of cluster data
        df_clusters = clusters_to_dataframe(self.clusters)
        df_clusters_fname = os.path.join(outdir, 'clusters.csv')
        df_clusters.to_csv(df_clusters_fname, sep=',', index=False, float_format='%.9f')
        log.info("Saved Cluster data to {!r}".format(df_peaks_fname))

        ### CSV of basic candidate parameters
        df_cands = pandas.DataFrame.from_dict([
            cand.params for cand in self.candidates
        ])
        df_cands_fname = os.path.join(outdir, 'candidates.csv')
        df_cands.to_csv(df_cands_fname, sep=',', index=False, float_format='%.9f')

        ### Candidates and candidate plots
        log.info("Writing candidate files")
        pool = multiprocessing.Pool(processes=self.config['processes'])
        writer = CandidateWriter(outdir, plot=self.config['plot_candidates'])
        arglist = [(rank, cand) for rank, cand in enumerate(self.candidates)]
        pool.map(writer, arglist)

        log.info("Data products written")       

    @timing
    def process(self, files, outdir):
        self.prepare(files)
        self.search()
        self.cluster_peaks()
        self.flag_harmonics()

        # NOTE: apply candidate filters *after* harmonic flagging
        # In case, for example, a bright zero DM candidate has some harmonics
        # that end up just above the DM threshold
        self.apply_candidate_filters()
        self.build_candidates()
        self.save_products(outdir=outdir)

    @classmethod
    def from_yaml_config(cls, fname):
        log.debug("Creating pipeline from config file: {}".format(fname))
        with open(fname, 'r') as fobj:
            conf = yaml.safe_load(fobj)
        log.debug("Pipeline configuration: {}".format(json.dumps(conf, indent=4)))
        return cls(conf)


###############################################################################


help_formatter = lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=16)


def parse_arguments():
    def outdir(path):
        """ Function that checks the outdir argument """
        if not os.path.isdir(path):
            msg = "Specified output directory {!r} does not exist".format(path)
            raise argparse.ArgumentTypeError(msg)
        return path

    parser = argparse.ArgumentParser(
        formatter_class=help_formatter,
        description=f"Search multiple DM trials with the riptide end-to-end FFA pipeline."
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        required=True,
        help="Pipeline configuration file",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=outdir,
        default=os.getcwd(),
        help="Output directory for the data products",
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default='DEBUG',
        help="Logging level for the riptide logger",
        choices=['DEBUG', 'INFO', 'WARNING']
    )
    parser.add_argument(
        "--log-timings",
        action='store_true',
        help="If this flag is specified, log the execution times of all major functions"
    )
    parser.add_argument(
        '--version', action='version', version=__version__
        )
    parser.add_argument("files", type=str, nargs="+", help="Input file(s) of the right format")
    args = parser.parse_args()
    return args


def run_program():
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s %(filename)18s:%(lineno)-4s %(levelname)-8s %(message)s"
        )
    args = parse_arguments()

    logging.getLogger('matplotlib').setLevel('WARNING') # otherwise it gets annoying
    logging.getLogger('riptide').setLevel(args.log_level)

    if args.log_timings:
        logging.getLogger('riptide.timing').setLevel('DEBUG')
    else:
        logging.getLogger('riptide.timing').setLevel('WARNING')

    pipeline = Pipeline.from_yaml_config(args.config)
    pipeline.process(args.files, args.outdir)

    # If you have seen the movie "The Martian" and always wanted to look like 
    # an actual scientist to your friends and family. Thank me later.
    log.info("CALCULATIONS CORRECT")


# NOTE: main() is the entry point of the console script
def main():
    ### Select non-interactive backend
    # matplotlib.use('Agg') would not work here, due to importing order 
    # the console_scripts entry point design means that 'riptide' is always imported first,
    # importing everything else in riptide's __init__.py, which ends up setting the backend
    # before the first line of this script is reached
    # Another alternative is to call the pipeline command with the MPLBACKEND=Agg prefix

    # NOTE: We don't do this at the top of the script, in case someone wants
    # to import the Pipeline class without switching backends
    import matplotlib.pyplot as plt
    plt.switch_backend('Agg')

    # NOTE (IMPORTANT): Force all numpy libraries to use a single thread/CPU
    # Each DM trial is assigned to a different process, and for optimal 
    # performance, each process should be limited to 1 CPU
    with threadpoolctl.threadpool_limits(limits=1):
        run_program()


if __name__ == '__main__':
    main()