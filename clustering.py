import numpy as np

def iterslices(indices):
    # Special case where there is only one cluster
    if not len(indices):
        yield slice(None, None)
        raise StopIteration

    # Multiple clusters case
    else:
        yield slice(None, indices[0])
        for ii, jj in zip(indices[:-1], indices[1:]):
            yield slice(ii, jj)
        yield slice(indices[-1], None)
        raise StopIteration

def cluster_1d(x, radius):
    """ Perform clustering on 1D data. Elements of 'x' whose
    absolute difference is less than 'radius' are considered
    part of the same cluster. """
    if not len(x):
        return []

    # Sort input data, compute differences between
    # consecutive elements, spot those differences
    # that exceeed the clustering radius
    x = np.asarray(x)
    order = x.argsort()
    y = x[order]
    dy = np.diff(y)
    diffmask = dy > radius

    # Indices that mark the end of a cluster
    ibreaks = np.where(diffmask)[0] + 1

    clusters = [
        order[sl]
        for sl in iterslices(ibreaks)
        ]
    return clusters
