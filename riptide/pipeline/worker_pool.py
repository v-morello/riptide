import logging
import multiprocessing

from riptide import TimeSeries, ffa_search, find_peaks


log = logging.getLogger('riptide.worker_pool')


class WorkerPool(object):
    """ 
    deredden_params : dict
    range_confs : list of dicts
        List of dicts from the 'ranges' section of the YAML config file
    loader : func
        Function that takes a file path as its only argument, and returns
        a TimeSeries object
    processes : int
        Number of parallel processes
    fmt : str
        TimeSeries file format
    """

    TIMESERIES_LOADERS = {
        'sigproc': TimeSeries.from_sigproc,
        'presto': TimeSeries.from_presto_inf
    }

    def __init__(self, deredden_params, range_confs, processes=1, fmt='presto'):
        self.deredden_params = deredden_params
        self.range_confs = range_confs
        self.loader = self.TIMESERIES_LOADERS[fmt]
        self.processes = int(processes)

    def process_fname_list(self, fnames):
        pool = multiprocessing.Pool(processes=self.processes)
        # results is a list of lists of Detections
        results = pool.map(self.process_fname, fnames)
        # NOTE: don't forget to close the pool to free up RAM
        # NOTE: and don't forget to join, otherwise the coverage module
        # does not properly report coverage for sub-processes spawned by
        # the pool
        pool.close()
        pool.join()
        return [det for dlist in results for det in dlist]

    def process_fname(self, fname):
        allpeaks = []
        ts = self.loader(fname)
        dm = ts.metadata['dm']
        log.debug("Searching DM = {:.3f}".format(dm))

        # Make pre-processing common to all ranges to save time
        ts = ts.deredden(
            self.deredden_params['rmed_width'], 
            minpts=self.deredden_params['rmed_minpts']
            )
        ts = ts.normalise()

        for conf in self.range_confs:
            kw_search = dict(conf['ffa_search'])
            kw_search.update({
                'deredden': False,
                'already_normalised': True
                })
            tsdr, pgram = ffa_search(ts, **kw_search)
            peaks, polycos = find_peaks(pgram, **conf['find_peaks'])
            allpeaks.extend(peaks)
            del tsdr, pgram, peaks, polycos # Free RAM ASAP
        log.debug(f"Done searching DM = {dm:.3f}, peaks found: {len(allpeaks)}")
        return allpeaks