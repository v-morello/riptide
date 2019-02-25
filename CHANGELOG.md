# Changelog
All notable changes to this project will be documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## 0.0.2rc1 - In preparation
### Fixed
- When normalising TimeSeries, use a float64 accumulator when calculating mean and variance. This avoid saturation issues on data with high values, e.g. from 8-bit Parkes UWL observations.

## 0.0.1 - 2018-10-25
### Added
- First stable release of riptide. This is the version that will be run on the LOTAAS survey.