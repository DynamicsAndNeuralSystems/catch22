# _catch22_ - CAnonical Time-series CHaracteristics

[![DOI](https://zenodo.org/badge/146194807.svg)](https://zenodo.org/badge/latestdoi/146194807)

## About

This is a collection of 22 time series features contained in the [_hctsa_](https://github.com/benfulcher/hctsa) toolbox coded in C.
Features were selected by their classification performance across a collection of 93 real-world time-series classification problems (according to the [`op_importance`](https://github.com/chlubba/op_importance) repository).

For information on how this feature set was constructed see our open-access paper:

* C.H. Lubba, S.S. Sethi, P. Knaute, S.R. Schultz, B.D. Fulcher, N.S. Jones. [_catch22_: CAnonical Time-series CHaracteristics](https://doi.org/10.1007/s10618-019-00647-x). *Data Mining and Knowledge Discovery* **33**, 1821 (2019).

For information on the full set of over 7000 features, see the following open-access publications:

* B.D. Fulcher and N.S. Jones. [_hctsa_: A computational framework for automated time-series phenotyping using massive feature extraction](http://www.cell.com/cell-systems/fulltext/S2405-4712\(17\)30438-6). *Cell Systems* **5**, 527 (2017).
* B.D. Fulcher, M.A. Little, N.S. Jones [Highly comparative time-series analysis: the empirical structure of time series and their methods](http://rsif.royalsocietypublishing.org/content/10/83/20130048.full). *J. Roy. Soc. Interface* **10**, 83 (2013).

## Installation: Python, Matlab, R, Julia, and compiled C

The fast, C-coded functions in this repository can be used in Python, Matlab, and R following the [detailed installation instructions on the wiki](https://github.com/chlubba/catch22/wiki/Installation-and-Testing).

Julia users can use [this Julia package](https://github.com/brendanjohnharris/Catch22.jl) to evaluate the _catch22_ feature set.

## Using _catch22_

- For more information about _catch22_, including a list of publications using the feature set, see the [_catch22_ wiki](https://github.com/chlubba/catch22/wiki).
- __Important Note:__ _catch22_ features only evaluate _dynamical_ properties of time series and do not respond to basic differences in the location (e.g., mean) or spread (e.g., variance).
  - If you think features of the raw distribution may be important for your application, we suggest you add them (in the simplest case, two additional features: the mean and standard deviation) to this feature set.
- Note that time series are _z_-scored internally which means e.g., constant time series will lead to `NaN` outputs.
