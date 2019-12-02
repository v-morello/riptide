import logging
import multiprocessing

from riptide import ffa_search, find_peaks


log = logging.getLogger('riptide.worker_pool')


class WorkerPool(object):
    """ 
    deredden_params: dict
    range_confs: list of dicts
        List of dicts from the 'ranges' section of the YAML config file
    processes: int
        Number of parallel processes
    """
    def __init__(self, deredden_params, range_confs, processes):
        self.deredden_params = deredden_params
        self.range_confs = range_confs
        self.processes = int(processes)

    def process_chunk(self, chunk):
        pool = multiprocessing.Pool(processes=self.processes)
        # results is a list of lists of Detections
        results = pool.map(self.process_time_series, chunk)
        # NOTE: don't forget to close the pool to free up RAM
        pool.close()
        return [det for dlist in results for det in dlist]

    def process_time_series(self, ts):
        allpeaks = []
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
            __, __, pgram = ffa_search(ts, **kw_search)
            peaks, polycos = find_peaks(pgram, **conf['find_peaks'])
            allpeaks.extend(peaks)
        return allpeaks

