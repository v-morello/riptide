### Non-standard imports
import pandas
from numpy import ceil

### Local module imports
from .ffautils import generate_width_trials


class ProcessingPlan(object):
    """ Defines the sequence of downsamplings and partial searches to perform
    on a time series. Meant to be stored as a class member of a Periodogram
    object. """
    _column_names = [
        'dsfactor', 'tsamp', 'bins_min', 'bins_max', 'period_min', 'period_max'
        ]

    def __init__(self, attrs, widths, steps):
        """ This should not be called directly. Use ProcessingPlan.create() instead. """
        self._widths = widths
        self._steps = steps
        self._attrs = attrs

    @property
    def widths(self):
        return self._widths

    @property
    def steps(self):
        return self._steps

    @property
    def tsamp(self):
        return self._attrs['tsamp']

    @property
    def nsamp(self):
        return self._attrs['nsamp']

    @property
    def bins_min(self):
        return int(self._steps.iloc[0].bins_min)

    @property
    def bins_max(self):
        return int(self._steps.iloc[0].bins_max)

    @property
    def bins_avg(self):
        """ Average number of phase bins used in the search. """
        return int(round(self.bins_min + self.bins_max - 1) / 2.0)

    @property
    def period_min(self):
        return self._steps.iloc[0].period_min

    @property
    def period_max(self):
        return self._steps.iloc[-1].period_max

    @staticmethod
    def create(nsamp, tsamp, period_min=1.0, period_max=30.0, fpmin=8, bins_min=240, bins_max=260, ducy_max=0.20, wtsp=1.5):
        """ Create a new ProcessingPlan instance

        Parameters
        ----------
        nsamp : int
            Number of samples in the time series to be processed with this plan.
        tsamp : float
            Sampling time of the time series to be processed.
        period_min : float
            Minimum period to search.
        """
        ### Check input properly
        tsamp = float(tsamp)
        nsamp = int(nsamp)
        fpmin = int(fpmin)
        tobs = tsamp * nsamp
        bins_min = int(bins_min)
        bins_max = int(bins_max)

        if not bins_max > bins_min:
            raise ValueError('bins_max must be > bins_min')

        if tsamp * bins_min > period_min:
            raise ValueError('period_min must be larger than tsamp x bins_min')

        # Cap period_max so that we do not FFA transform less than 'fpmin' pulse periods
        period_max = min(period_max, tobs / fpmin)

        # Width trials in number of bins
        bins_avg = int(round(bins_min + bins_max - 1) / 2.0)
        widths = generate_width_trials(bins_avg, ducy_max=ducy_max, wtsp=wtsp)

        # Downsampling factor for first plan step
        ds = period_min / bins_min / tsamp
        ds_growth = bins_max / float(bins_min)

        ### Compute plan steps
        steps = []
        while True:
            current_period_min = bins_min * ds * tsamp
            current_period_max = bins_max * ds * tsamp
            current_step = (ds, ds * tsamp, bins_min, bins_max, current_period_min, current_period_max)
            steps.append(current_step)

            if current_period_max >= period_max:
                break

            ds *= ds_growth

        steps = pandas.DataFrame(steps, columns=ProcessingPlan._column_names)

        # Edit last step to stop the search as close as possible to period_max
        last_step = steps.iloc[-1]
        last_bins_max = int(ceil(period_max / last_step.tsamp))
        last_period_max = last_bins_max * last_step.tsamp

        ix = len(steps) - 1
        steps.loc[ix, 'bins_max'] = last_bins_max
        steps.loc[ix, 'period_max'] = last_period_max

        # Pack basic attributes worth keeping as class members into a dictionary
        attrs = {
            'nsamp' : nsamp,
            'tsamp' : tsamp,
            }
        return ProcessingPlan(attrs, widths, steps)

    def __repr__(self):
        lines = [
            'Number of samples  : %d' % self.nsamp,
            'Sampling time      : %.6e' % self.tsamp,
            'Width trials (bins): %s' % (', '.join(map(str, self.widths))),
            '',
            repr(self.steps)
            ]
        return '\n'.join(lines)

    @classmethod
    def from_dict(cls, items):
        return cls(items['attrs'], items['widths'], items['steps'])

    def to_dict(self):
        return {
            'attrs': self._attrs,
            'widths': self._widths,
            'steps': self._steps
        }

