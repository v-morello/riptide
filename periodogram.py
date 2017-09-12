import matplotlib.pyplot as plt

class Periodogram(object):
    """ Stores the raw output of the FFA search of a time series. """
    def __init__(self, periods, widths, snrs):
        self.periods = periods
        self.widths = widths
        self.snrs = snrs.reshape(periods.size, widths.size)
        
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

    def save_hdf5(self, fname):
        pass

    @classmethod
    def load_hdf5(cls, fname):
        pass
