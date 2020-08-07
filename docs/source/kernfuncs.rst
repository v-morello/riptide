FFA kernel functions
====================

The python interface of ``riptide`` exposes some lower-level functions related to calculating the folding transform (interchangeably called FFA transform) of input data at some base integer period.
There is also a kernel function to calculate S/N. These are:  

* ``ffa2``: FFA transform of a two-dimensional input that represents a pulse stack. The ``m`` lines of the input represent pulses in chronological order, and the ``p`` columns represent the phase dimension  
* ``ffa1``: FFA transform of a one-dimensional input, that represents a time series. The function simply selects the largest number of entire pulses that fit in the data, reshapes them into a two-dimensional array, and calls ``ffa2()``  
* ``ffafreq``: Returns the trial folding frequencies corresponding to every line in the output of an FFA transform  
* ``ffaprd``: Same as ``ffafreq``, but returns trial periods instead
* ``boxcar_snr``: Compute the S/N ratio of pulse profile(s) by concolving them with a range of boxcar filters with different widths.

See :ref:`API Reference` for further details.