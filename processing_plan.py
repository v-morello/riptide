### Non-standard imports
import pandas
from numpy import ceil

### Local module imports
from .ffautils import generate_width_trials



class ProcessingPlan(object):
    """ Defines the sequence of downsamplings and partial searches to perform on a time series. """
    def __init__(self, nsamp, tsamp, period_min=1.0, period_max=30.0, fpmin=8, bins_min=240, bins_max=260, ducy_max=0.20, wtsp=1.5):
        """
        Parameters:
        -----------
            nsamp: int
                Number of samples in the time series to be processed with this plan.
            tsamp: float
                Sampling time of the time series to be processed.
            period_min: float
                Minimum period to search.
        """
        self.tsamp = float(tsamp)
        self.nsamp = int(nsamp)
        self.fpmin = int(fpmin)
        self.tobs = self.tsamp * self.nsamp
        bins_min = int(bins_min)
        bins_max = int(bins_max)

        ### Check input
        if not bins_max > bins_min:
            raise ValueError('bins_max must be > bins_min')

        if self.tsamp * bins_min > period_min:
            raise ValueError('period_min must be larger than tsamp x bins_min')

        # Cap period_max so that we do not FFA transform less than 'fpmin' pulse periods
        period_max = min(period_max, self.tobs / self.fpmin)

        # average number of bins
        bins_avg = int(round(bins_min + bins_max - 1) / 2.0)
        
        # Width trials in number of bins
        self.widths = generate_width_trials(bins_avg, ducy_max=ducy_max, wtsp=wtsp)
        
        ### Initial downsampling factor
        ds_ini = period_min / bins_min / tsamp
        self.ds_ini = ds_ini
        self.nsamp_ds = int(self.nsamp / self.ds_ini)
        self.tsamp_ds = self.tsamp * self.ds_ini
        
        ds_growth = bins_max / float(bins_min)
        
        ### Compute plan steps
        steps = []
        ds = 1.0 # downsampling to apply after initial downsampling
        while True:
            current_period_min = bins_min * ds * ds_ini * tsamp
            current_period_max = bins_max * ds * ds_ini * tsamp
            current_step = (ds, ds * ds_ini * tsamp, bins_min, bins_max, current_period_min, current_period_max)
            steps.append(current_step)
            
            if current_period_max >= period_max:
                break

            ds *= ds_growth
                
        columns = ['dsfactor', 'tsamp', 'bins_min', 'bins_max', 'period_min', 'period_max']
        self.steps = pandas.DataFrame(steps, columns=columns)
        
        # Edit last step to stop the search as close as possible to period_max
        last_step = self.steps.iloc[-1]
        last_bins_max = int(ceil(period_max / last_step.tsamp))
        last_period_max = last_bins_max * last_step.tsamp
        
        ix = len(self.steps) - 1
        self.steps.loc[ix, 'bins_max'] = last_bins_max
        self.steps.loc[ix, 'period_max'] = last_period_max

    def __repr__(self):
        lines = [
            'Number of samples  : %d' % self.nsamp,
            'Sampling time      : %.6e' % self.tsamp,
            '',
            'Initial downsampling factor: %.6f' % self.ds_ini,
            'Num. samples after initial downsampling : %d' % self.nsamp_ds,
            'Sampling time after initial downsampling: %.6e' % self.tsamp_ds,
            '',
            'Width trials (bins): %s' % (', '.join(map(str, self.widths))),
            '',
            repr(self.steps)
            ]
        return '\n'.join(lines)
