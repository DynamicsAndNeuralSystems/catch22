# _catch22_ - CAnonical Time-series CHaracteristics
This is a collection of 22 time series features contained in the [_hctsa_](https://github.com/benfulcher/hctsa) toolbox coded in C. Features were selected by their classification performance across a collection of 93 real-world time-series classification problems.

For information on how this feature set was constructed see our preprint:

* C.H. Lubba, S.S. Sethi, P. Knaute, S.R. Schultz, B.D. Fulcher, N.S. Jones. [_catch22_: CAnonical Time-series CHaracteristics](https://arxiv.org/abs/1901.10200). arXiv (2019).

For information on the full set of over 7000 features, see the following (open) publications:

* B.D. Fulcher and N.S. Jones. [_hctsa_: A computational framework for automated time-series phenotyping using massive feature extraction](http://www.cell.com/cell-systems/fulltext/S2405-4712\(17\)30438-6). *Cell Systems* **5**, 527 (2017).
* B.D. Fulcher, M.A. Little, N.S. Jones [Highly comparative time-series analysis: the empirical structure of time series and their methods](http://rsif.royalsocietypublishing.org/content/10/83/20130048.full). *J. Roy. Soc. Interface* **10**, 83 (2013).

# Using the _catch22_-features from Python, Matlab and R

The fast C-coded functions in this repository can be used in Python, Matlab, and R following the instructions below. Time series are z-scored internally which means e.g., constant time series will lead to NaN outputs. The wrappers are only tested on OS X so far and require Clang.

## Python

Installation of the Python wrapper differs slightly between Python 2 and 3.

### Installation Python 2

Go to the directory `wrap_Python` and run the following

```
python setup.py build
python setup.py install
```

or alternatively, using pip, go to main directory and run

```
pip install -e wrap_Python
```

### Installation Python 3

Only manual installation through `distutils`

```
python3 setup_P3.py build
python3 setup_P3.py install
```

### Test Python 2 and 3

To test that the _catch22_ wrapper was installed successfully and works run (NB: replace `python` with `python3` for Python 3):

```
$ python testing.py
```

The module is now available under the name `catch22`. Each feature function can be accessed individually and takes arrays as tuple or lists (not `Numpy`-arrays). E.g., for loaded data, `tsData` in Python:

```python
import catch22
catch22.CO_f1ecac(tsData)
```

All features are bundeled in the method `catch22_all` which also accepts numpy arrays and gives back a dictionary containing the entries `catch22_all['names']` for feature names and `catch22_all['values']` for feature outputs.

```python
from catch22 import catch22_all
catch22_all(tsData)
```

## R

This assumes your have `R` installed and the package `Rcpp` is available.

Copy all `.c`- and `.h`-files from `./C_Functions` to `./wrap_R/src`. Then go to the directory `./wrap_R` and run

```
R CMD build catch22
R CMD INSTALL catch22_x.y.tar.gz
```

To test if the installation was successful, navigate to `./wrap_R` in the console and run:

```
$ Rscript testing.R
```

The module is now available in `R` as `catch22`. Single functions can be accessed by their name, all functions are bundeled as `catch22_all` which can be called with a data vector `tsData` as an argument and gives back a data frame with the variables `name` for feature names and `values` for feature outputs:

```
library(catch22)
catch22_out = catch22_all(tsData);
print(catch22_out)
```

## Matlab

Go to the `wrap_Matlab` directory and call `mexAll` from within Matlab. Include the folder in your Matlab path to use the package.

To test, navigate to the `wrap_Matlab` directory from within Matlab and run:

```
testing
```

All feature can be called individually, e.g. `catch22_CO_f1ecac`. Alternatively, all features are bundeled in a function `catch22_all` which returns an array of feature outputs and, as a second output, a cell array of feature names. With loaded data `tsData`:

```
[vals, names] = catch22_all(data);
```

# Raw C

## Compilation

### OS X
```
gcc -o run_features main.c CO_AutoCorr.c DN_HistogramMode_10.c DN_HistogramMode_5.c DN_OutlierInclude.c FC_LocalSimple.c IN_AutoMutualInfoStats.c MD_hrv.c PD_PeriodicityWang.c SB_BinaryStats.c SB_CoarseGrain.c SB_MotifThree.c SB_TransitionMatrix.c SC_FluctAnal.c SP_Summaries.c butterworth.c fft.c helper_functions.c histcounts.c splinefit.c stats.c
```
### Ubuntu:
As for OS X but with `-lm` switch in from of every source-file name.

## Usage

### Single files

The compiled `run_features` program only takes one time series at a time. Usage is `./run_features <infile> <outfile>` in the terminal, where specifying `<outfile>` is optional, it prints to  `stdout` by default.

### Mutliple files

For multiple time series, put them – one file for each – into a folder `timeSeries` and call `./runAllTS.sh`. The output will be written into a folder `featureOutput`. Do change the permissions of `runAllTS.sh` to executable by calling `chmod 755 runAllTS.sh`.

### Output format

Each line of the output correponds to one feature; the three comma-separated entries per line correspond to feature value, feature name and feature execution time in milliseconds. E.g.
```
0.29910714285714, CO_Embed2_Basic_tau.incircle_1, 0.341000
0.57589285714286, CO_Embed2_Basic_tau.incircle_2, 0.296000
...
```

### Testing

Sample outputs for the time series `test.txt` and `test2.txt` are provided as `test_output.txt` and `test2_output.txt`. The first two entries per line should always be the same. The third one (execution time) will be different.
