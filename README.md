<p align="center"><img src="img/catch22_logo_square.png" alt="catch22 logo" height="220"/></p>

<h1 align="center"><em>catch22</em>: CAnonical Time-series CHaracteristics</h1>

<p align="center">
 	<a href="https://zenodo.org/badge/latestdoi/146194807"><img src="https://zenodo.org/badge/146194807.svg" height="20"/></a>
    <a href="https://www.gnu.org/licenses/gpl-3.0"><img src="https://img.shields.io/badge/License-GPLv3-blue.svg" height="20"/></a>
 	<a href="https://twitter.com/compTimeSeries"><img src="https://img.shields.io/twitter/url/https/twitter.com/compTimeSeries.svg?style=social&label=Follow%20%40compTimeSeries" height="20"/></a>
</p>

_catch22_ is a collection of 22 time-series features coded in C that can be run from Python, R, Matlab, and Julia, licensed under the [GNU GPL v3 license](http://www.gnu.org/licenses/gpl-3.0.html) (or later).
The _catch22_ features are a high-performing subset of the over 7000 features in [_hctsa_](https://github.com/benfulcher/hctsa).

The features were selected based on their classification performance across a collection of 93 real-world time-series classification problems, as described in our open-access paper, [&#x1F4D7; Lubba et al. (2019). _catch22_: CAnonical Time-series CHaracteristics](https://doi.org/10.1007/s10618-019-00647-x).

## [&#x1F4D9;&#x1F4D8;&#x1F4D7;___catch22_ documentation__](https://time-series-features.gitbook.io/catch22/)

There is [comprehensive documentation](https://time-series-features.gitbook.io/catch22/) for _catch22_, including:

- Installation instructions (across C, python, R, Julia, and Matlab)
- Information about the theory behind and behavior of each of the features,
- A list of publications that have used or extended _catch22_
- And more :yum:

## Installation and Usage in Python, R, Matlab, Julia, and compiled C

There are also native versions of this code for other programming languages:

- [Rcatch22](https://github.com/hendersontrent/Rcatch22) (R) `install.packages("Rcatch22")`
- [pycatch22](https://github.com/DynamicsAndNeuralSystems/pycatch22) (python) `pip install pycatch22`
- [Catch22.jl](https://github.com/brendanjohnharris/Catch22.jl) (Julia) `Pkg.add("Catch22")`

You can also use the C-compiled features directly, or in Matlab, following the [detailed installation instructions in the documentation](https://time-series-features.gitbook.io/catch22/).

## Acknowledgement :+1:

If you use this software, please read and cite this open-access article:

- &#x1F4D7; Lubba et al. [_catch22_: CAnonical Time-series CHaracteristics](https://doi.org/10.1007/s10618-019-00647-x), _Data Min Knowl Disc_ __33__, 1821 (2019).

## Performance Summary

Summary of the performance of the _catch22_ feature set across 93 classification problems, and a comparison to the [_hctsa_ feature set](https://github.com/benfulcher/hctsa) (cf. Fig. 4 from [our paper](https://doi.org/10.1007/s10618-019-00647-x)):

![](img/PerformanceComparisonFig4.png)

## Notes

- When presenting results using _catch22_, you must identify the version used to allow clear reproduction of your results. For example, `CO_f1ecac` was altered from an integer-valued output to a linearly interpolated real-valued output from v0.3.
- _catch22_ features only evaluate _dynamical_ properties of time series and do not respond to basic differences in the location (e.g., mean) or spread (e.g., variance).
- If the location and spread of the raw time-series distribution may be important for your application, you should apply the function argument `catch24 = true` (`TRUE` in R, `True` in Python) to your call to the _catch22_ function in the language of your choice. This will result in 24 features being calculated: the _catch22_ features in addition to mean and standard deviation.
- Time series are _z_-scored internally (for features other than mean and standard deviation), which means that, e.g., constant time series will lead to `NaN` outputs.
- Time-series data are taken as an ordered sequence of values (without time stamps). We assume an evenly sampled time series.
- See language-specific usage information in the [docs](https://time-series-features.gitbook.io/catch22/).
- The computational pipeline used to generate the _catch22_ feature set is in the [`op_importance`](https://github.com/chlubba/op_importance) repository.
