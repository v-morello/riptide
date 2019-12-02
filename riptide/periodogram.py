##### Non-standard imports #####
import numpy as np
import matplotlib.pyplot as plt

##### Local imports #####
from .metadata import Metadata
from .processing_plan import ProcessingPlan

class Periodogram(object):
    """ Stores the raw output of the FFA search of a time series. """
    def __init__(self, plan, periods, widths, snrs, metadata=None):
        self.plan = plan
        self.periods = periods
        self.widths = widths
        self.snrs = snrs.reshape(periods.size, widths.size)
        self.metadata = metadata if metadata is not None else Metadata({})

    @property
    def freqs(self):
        return 1.0 / self.periods

    @property
    def tobs(self):
        return self.metadata['tobs']

    @property
    def bins_avg(self):
        """ Average number of phase bins used in the search """
        return self.plan.bins_avg

    def to_dict(self):
        return {
            'plan': self.plan,
            'periods': self.periods,
            'widths': self.widths,
            'snrs': self.snrs,
            'metadata': self.metadata
        }

    @classmethod
    def from_dict(cls, items):
        return cls(items['plan'], items['periods'], items['widths'], items['snrs'], metadata=items['metadata'])

    def plot(self, iwidth=None):
        if iwidth is None:
            snr = self.snrs.max(axis=1)
        else:
            snr = self.snrs[:, iwidth]

        plt.plot(self.periods, snr, marker='o', markersize=2, alpha=0.5)
        plt.xlim(self.periods.min(), self.periods.max())
        plt.xlabel('Trial Period (s)', fontsize=16)
        plt.ylabel('S/N', fontsize=16)

        if iwidth is None:
            plt.title('Best S/N at any trial width', fontsize=18)
        else:
            width_bins = self.widths[iwidth]
            plt.title('S/N at trial width = %d' % width_bins, fontsize=18)

        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.grid(linestyle=':')
        plt.tight_layout()

    def display(self, iwidth=None, figsize=(20,5), dpi=100):
        plt.figure(figsize=figsize, dpi=dpi)
        self.plot(iwidth=iwidth)
        plt.show()

