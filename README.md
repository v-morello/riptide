### Overview

__riptide__ is an end-to-end pulsar searching package implementing the Fast Folding Algorithm (FFA). Its interface is entirely in python while the core algorithms are implemented in C. It is being developed with ease of use in mind, without sacrificing execution speed. The main features are / are planned to be:

- Reading time series data written by major pulsar software packages such as PRESTO and SIGPROC
- Saving / Loading data products to / from HDF5
- Carry metadata in a flexible format across pipeline stages and into HDF5 outputs.
- Make nice plots of everything

### The name

riptide loosely stands for "seaRchIng for Pulsars in the TIme Domain". Mostly because there are no cool names containing the letters "FFA". Also a subtle reference to Australia where I did my Masters.

### Dependencies

* numpy
* matplotlib
* pandas
* h5py
* astropy
* pyyaml

### Installation

Clone the repository, for example:
```bash
cd ~/software
git clone https://vmorello@bitbucket.org/vmorello/riptide.git
```

Build the C library libffa.so in the c_src/ subdirectory:

```bash
cd c_src
make all
```

And then make sure that the base directory of the riptide repository ("~/software/riptide" in this example) is in the PYTHONPATH.


### Docker image

The riptide Dockerfile is located in the 'docker' subdirectory. To build the image with the name 'riptide':

```bash
cd docker
docker build -t riptide .
```

Both python3 and ipython are installed within the docker image. To run it interactively:

```bash
docker run -it --rm riptide bash
```

Refer to your favourite docker cheat sheet for further information and advanced usage.


### Basic Usage

Coming soon.
