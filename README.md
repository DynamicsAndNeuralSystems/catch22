# _catch22_ - CAnonical Time-series CHaracteristics

[![DOI](https://zenodo.org/badge/146194807.svg)](https://zenodo.org/badge/latestdoi/146194807)

## About

_catch22_ is a collection of 22 time-series features coded in C.
The _catch22_ features are a high-performing subset of the over 7000 features in [_hctsa_](https://github.com/benfulcher/hctsa).

Features were selected based on their classification performance across a collection of 93 real-world time-series classification problems, as described in our open-access paper [&#x1F4D7;Lubba et al. (2019). _catch22_: CAnonical Time-series CHaracteristics](https://doi.org/10.1007/s10618-019-00647-x).
For the computational pipeline, see the [`op_importance`](https://github.com/chlubba/op_importance) repository.

For _catch22_-related information and resources, including a list of publications using _catch22_, see the [___catch22_ wiki__](https://github.com/chlubba/catch22/wiki).

## Installation: Python, R, Matlab, Julia, and compiled C

The fast, C-coded functions in this repository can be used in Python, Matlab, and R following the [detailed installation instructions on the wiki](https://github.com/chlubba/catch22/wiki/Installation-and-Testing).

Julia users can use [this Julia package](https://github.com/brendanjohnharris/Catch22.jl) to evaluate the _catch22_ feature set.

## Usage

- See language-specific usage information in the [wiki](https://github.com/chlubba/catch22/wiki/Installation-and-Testing).
- __Important Note:__ _catch22_ features only evaluate _dynamical_ properties of time series and do not respond to basic differences in the location (e.g., mean) or spread (e.g., variance).
  - If you think features of the raw distribution may be important for your application, we suggest you add them (in the simplest case, two additional features: the mean and standard deviation) to this feature set.
- Note that time series are _z_-scored internally which means e.g., constant time series will lead to `NaN` outputs.
