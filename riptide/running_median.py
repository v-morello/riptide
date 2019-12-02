import numpy
from bisect import insort, bisect_left


def running_median(data, width):
    if not width % 2:
        raise ValueError("width must be an odd number")
    l = width * [data[0]]
    mididx = (width - 1) // 2
    result = numpy.zeros_like(data)
    for idx, new_elem in enumerate(data):
        old_elem = data[max(0, idx - width)]
        del l[bisect_left(l, old_elem)]
        insort(l, new_elem)
        result[idx] = l[mididx]
    return result[2*mididx:]


def fast_running_median(data, width, min_points):
    """ Compute an approximate running median of data over large window sizes.

    Parameters
    ----------
    data : ndarray
        Input data
    width : int
        Required width of the running median window in number of samples
    min_points : int
        The running median is calculated of a time scrunched version of the
        input data to save time: minpts is the minimum number of
        scrunched samples that must fit in the running median window.
        Lower values make the running median calculation less accurate but
        faster, due to allowing a higher scrunching factor.
    """
    if width < 3:
        raise ValueError('width must be > 3')
    if min_points < 3:
        raise ValueError('min_points must be > 3')
    
    width = int(width)
    min_points = int(min_points)

    # num_points = effective width of running median on downsampled data
    dsfactor = max(int(width / float(min_points)), 1)
    num_points = int(numpy.ceil(width / float(dsfactor)))
    
    # Make sure num_points is an odd number, makes life easier when it comes
    # to padding the edges of the data
    num_points += int(not num_points % 2)

    ### Edge case: no downsampling needed
    if dsfactor <= 1:
        # num_points is our final window width even in this case
        y = numpy.pad(data, num_points // 2, 'median', stat_length=width)
        return running_median(y, num_points)
    
    ### Add 'width' elements to data on both edges.
    # Fill value at the beginning is the median of the first 'width' samples of data
    # Fill value at the end is the median of the last 'width' samples of data
    y = numpy.pad(data, width, 'median', stat_length=width)
    
    # x-coordinates of padded data in original referential
    x = numpy.arange(-width, data.size + width, 1, dtype=int)
    
    ### Downsample
    n = y.size - y.size % dsfactor
    yds = y[:n].reshape(-1, dsfactor).mean(axis=1)
    xds = x[:n].reshape(-1, dsfactor).mean(axis=1)
    
    ### Running Median
    rmed_ds = running_median(yds, num_points)
    x_rmed_ds = xds[numpy.arange(rmed_ds.size) + num_points//2]
    
    ### Upsample
    rmed = numpy.interp(numpy.arange(data.size), x_rmed_ds, rmed_ds)
    return rmed
