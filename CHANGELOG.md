# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 0.0.3 - 2019-11-30
### Added
- Full docstring for `ffa_search()`

### Changed
- Cleaner and faster FFA C kernels, which have been moved to separate files
- Slight optimisation to S/N calculation, where only the best value across all pulse phases is normalized, instead of normalizing at every phase trial
- S/N calculation separated into smaller functions
- Removed OpenMP multithreading from S/N calculation, it was slower in most usual cases. The benefits were visible only for very large input data. As a result, `ffa_search()` does not accept the 'threads' keyword anymore, and the 'threads' keyword has also been removed from the pipeline configuration files (in the 'search' section). Parallelism only happens at the TimeSeries level, i.e. one process per TimeSeries.

### Fixed
- Reinstated the `--fast-math` compilation option that had been accidentally removed in v0.0.2. The code is now much faster.


## 0.0.2 - 2019-11-05
### Added
- `riptide` is now properly packaged with a `setup.py` script. Updated installation instructions in `README.md`
- Updated Dockerfile. Build docker image with `make docker` command.

### Fixed
- When normalising TimeSeries, use a float64 accumulator when calculating mean and variance. This avoid saturation issues on data with high values, e.g. from 8-bit Parkes UWL observations.

### Changed
- Improved candidate plots: docstring for `Candidate` class, DM unit on plots, option to subtract the baseline value of the profile before displaying it, option to plot the profile as either a bar or line chart.

## 0.0.1 - 2018-10-25
### Added
- First stable release of riptide. This is the version that will be run on the LOTAAS survey.
