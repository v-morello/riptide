import numpy as np


def cluster1d(x, r, already_sorted=False):
    """
    Find clusters in one dimensional data using a simple friends of friends
    algorithm. Two points are in the same cluster if they are within no more
    than a distance 'r' from each other

    Parameters
    ----------
    x : ndarray
        Input data
    r : float
        Clustering radius
    already_sorted : bool, optional
        True if 'x' is already sorted in increasing or decreasing order.
        This can save some compute time.

    Returns
    -------
    clusters : list
        List of clusters. Each cluster is an array of indices of points in 'x'
    """
    if not len(x):
        return []

    if not already_sorted:
        indices = x.argsort()
        diff = np.diff(x[indices])
    else:
        indices = np.arange(len(x))
        diff = np.diff(x)

    # NOTE: diff is the sequence of consecutive differences
    # of x AFTER it has been sorted
    # Indices 
    ibreaks = np.where(abs(diff) > r)[0]

    # In this case, there is only one cluster
    if not len(ibreaks):
        return [indices]

    # Cluster bounds indices
    ibounds = np.concatenate(([0], ibreaks+1, [len(x)]))
    
    clusters = [
        indices[start:end]
        for start, end in zip(ibounds[:-1], ibounds[1:])
    ]
    return clusters