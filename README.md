## Overview

__riptide__ is an end-to-end pulsar searching package implementing the Fast Folding Algorithm (FFA). Its interface is entirely in python while the core algorithms are implemented in C. It is being developed with ease of use in mind, without sacrificing execution speed. The main features are:

- Reading time series data written by major pulsar software packages such as PRESTO and SIGPROC
- Saving / Loading data products to / from HDF5
- Carry metadata in a flexible format across pipeline stages and into HDF5 outputs.
- Make nice plots of everything

__riptide__ loosely stands for "seaRchIng for Pulsars in the TIme Domain". Mostly because there are no cool names containing the letters "FFA". Also a subtle reference to Australia where I did my Masters.


## Installation

Clone the repository:
```bash
git clone https://vmorello@bitbucket.org/vmorello/riptide.git
```

And then in the base directory of `riptide` run

```bash
make install
```

This simply runs ``pip install`` in [editable mode](https://pip.pypa.io/en/latest/reference/pip_install/#editable-installs), which means you can freely edit the code. It also installs any required dependencies with ``pip`` that are not present already.


## Docker image

The riptide Dockerfile is located in the 'docker' subdirectory. To build the image with the name 'riptide':

```bash
cd docker
docker build -t riptide .
```

Both python and ipython are installed within the docker image. To start a container:

```bash
docker run -it --rm riptide
```

Refer to your favourite docker cheat sheet for further information and advanced usage.


## Basic Usage: searching a single time series

The core functionality of riptide is to take a single time series as an input and calculate its *periodogram*: the S/N as a function of trial period and trial pulse width. This is what we will demonstrate in this section.

#### Load a time series, or create a fake one

Time series (or single DM trials) are encapsulated in the aptly named `TimeSeries` class. We can either create an artificial train of pulses (plus some white noise) for test purposes, or load some real data created with some popular pulsar software packages.

```python
from riptide import TimeSeries
# Generate an artificial train of pulses with a white noise background
# 600s of data, 256us sampling time, a period of pi seconds and a duty cycle of 2%
tseries_fake = TimeSeries.generate(length=600.0, tsamp=256e-6, period=3.14159, ducy=0.02)

# We can also load an existing time series created by PRESTO or SIGPROC
tseries_presto = TimeSeries.from_presto_inf("~/work/pulsars/presto_time_series/J1855+0307/J1855+0307_DM400.00.inf")
tseries_sigproc = TimeSeries.from_sigproc("~/work/pulsars/sigproc_time_series/J0636-4549.sigproc.tim")
```


#### Computing and manipulating periodograms

Periodograms are computed with the `ffa_search` function, which takes as input a `TimeSeries` and many keyword arguments. The search period range is specified via `period_min` and `period_max`. The duty cycle resolution of the search is set by `bins_min` and `bins_max`; due to how the search is performed under the hood, it is recommended to make `bins_max` approximately 10% higher than `bins_min`.

```python
from riptide import ffa_search

# Compute periodogram
ts, plan, pgram = ffa_search(tseries_sigproc, rmed_width=4.0, period_min=1.0, period_max=10.0, bins_min=240, bins_max=260)

# Plot S/N vs. trial period
pgram.display()
```

`ffa_search` returns three outputs:

1. The de-reddened and normalised copy of the input time series that was actually searched
2. The search plan that was followed, specifying how the data was iteratively downsampled across the entire search period ramge
3. The output `Periodogram`

Periodograms are actually two-dimensional: they represent a S/N as a function of both trial period and trial width, as shown below. The `display()` method only shows S/N for the best trial width.

```python
# Trial periods in seconds, trial widths in number of phase bins, output S/N
# pgram.display() only shows pgram.snrs.max(axis=1)
print(pgram.periods.size, pgram.widths.size, pgram.snrs.shape)
> 124778 10 (124778, 10)
```

#### Important detail: the metadata attribute

All `TimeSeries` objects and all derived data products (periodogram, pulsar search candidates, etc.) in riptide have a `metadata` dictionary carrying whatever information provided by the software package that created the input time series data. There is of course no header standardization for such data in pulsar astronomy, and the information contained in `metadata` will therefore vary across software packages and observatories. We do however attempt to guarantee some metadata uniformity in riptide, by always enforcing the presence of the following metadata keys and their associated data types:

* `dm`: `float`, dispersion measure of the input data
* `fname`: `str`, original file name
* `mjd`: `float`, epoch of observation
* `source_name`: `str`
* `skycoord`: `astropy.SkyCoord`, source coordinates
* `tobs`: `float`, integration time in seconds

If these required attributes were not provided by the original creator of the time series, they are set to the special value `None`. Here's for example the metadata for the SIGPROC time series loaded in code snippet above:

```python
print(tseries_sigproc.metadata)

{'az_start': 0.0,
 'barycentric': 0,
 'data_type': 2,
 'dm': 26.31,
 'fch1': 1581.8046875,
 'fname': '/home/vince/work/pulsars/sigproc_time_series/J0636-4549.sigproc.tim',
 'machine_id': 10,
 'mjd': 56771.1303125,
 'nbits': 32,
 'nchans': 1,
 'nifs': 1,
 'refdm': 26.31,
 'skycoord': <SkyCoord (ICRS): (ra, dec) in deg
    (99.17595833, -45.73472222)>,
 'source_name': 'G255.3-21.9_s',
 'src_dej': -454405.0,
 'src_raj': 63642.23,
 'telescope_id': 4,
 'tobs': 557.040896,
 'tsamp': 6.4e-05,
 'tstart': 56771.1303125,
 'za_start': 0.0}
```

## Large scale survey usage: searching multiple DM trials

Multiple DM trial searches are managed by the `Pipeline` class in the `riptide.pipelines` submodule. Documentation coming soon.
