import pandas


class PeakCluster(list):
    """
    Basic list subclass to store a cluster of Peak objects.

    Parameters
    ----------
    peaks: iterable
        List or iterable of Peak objects
    rank: int or None, optional
        Rank within the search, 0 means brightest
        (default: None)
    parent_fundamental: PeakCluster or None, optional
        Parent fundamental PeakCluster, can be set later by the harmonic
        flagging procedure. None means that the cluster has no parent, and
        is thus a fundamental itself.
        (default: None)
    hfrac: fractions.Fraction or None, optional
        If there is a parent fundamental, this is the ratio between the
        cluster's frequency and its fundamental's frequency
        (default: None)
    """

    def __init__(self, peaks, rank=None, parent_fundamental=None, hfrac=None):
        super().__init__(peaks)
        self.rank = rank
        self.parent_fundamental = parent_fundamental
        self.hfrac = hfrac

    @property
    def is_harmonic(self):
        """Return whether this cluster is flagged as a harmonic."""
        return self.parent_fundamental is not None

    @property
    def centre(self):
        """Return the highest-S/N peak in the cluster."""
        return max(self, key=lambda peak: peak.snr)

    def summary_dataframe(self):
        """
        Return a DataFrame containing the parameters of member Peaks.

        The columns are the keys of the dictionary returned by the
        Peak.summary_dict() method.
        """
        return pandas.DataFrame.from_dict([peak.summary_dict() for peak in self])

    def summary_dict(self):
        """Return a dictionary summarizing the cluster."""
        return {
            **self.centre.summary_dict(),
            "npeaks": len(self),
            # NOTE: use defaults instead of None when there is no fundamental.
            # This works around a pandas.DataFrame limitation where columns
            # with missing values MUST be of type float, while these must be int.
            "rank": self.rank,
            "hfrac_num": self.hfrac.numerator if self.is_harmonic else 0,
            "hfrac_denom": self.hfrac.denominator if self.is_harmonic else 0,
            "fundamental_rank": (
                self.parent_fundamental.rank if self.is_harmonic else self.rank
            ),
        }

    def __str__(self):
        """Return a human-readable representation."""
        name = type(self).__name__
        return f"{name}(size={len(self)}, centre={self.centre})"

    def __repr__(self):
        """Return the developer representation."""
        return str(self)


def clusters_to_dataframe(clusters):
    """
    Convert PeakCluster objects to a DataFrame.

    Include a summary of their attributes and harmonic parameters, sorted by
    decreasing S/N.
    """
    clusters = sorted(clusters, key=lambda c: c.centre.snr, reverse=True)
    df = pandas.DataFrame.from_dict([cl.summary_dict() for cl in clusters])

    # Re-order columns
    columns = [
        "rank",
        "period",
        "dm",
        "snr",
        "ducy",
        "freq",
        "npeaks",
        "hfrac_num",
        "hfrac_denom",
        "fundamental_rank",
    ]
    df = df[columns]
    return df
