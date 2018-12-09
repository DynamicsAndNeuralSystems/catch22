# _catch22_ - CAnonical Time-series CHaracteristics
This is a collection of 22 time series features contained in the [_hctsa_](https://github.com/benfulcher/hctsa) toolbox coded in C.

* B.D. Fulcher and N.S. Jones. [_hctsa_: A computational framework for automated time-series phenotyping using massive feature extraction](http://www.cell.com/cell-systems/fulltext/S2405-4712\(17\)30438-6). *Cell Systems* **5**, 527 (2017).
* B.D. Fulcher, M.A. Little, N.S. Jones [Highly comparative time-series analysis: the empirical structure of time series and their methods](http://rsif.royalsocietypublishing.org/content/10/83/20130048.full). *J. Roy. Soc. Interface* **10**, 83 (2013).

# Using the _catch22_-features from Python, Matlab and R

The fast C-coded functions in this repository can be used in Python, Matlab and R following the instructions below. Only tested on OS X so far.

## Python

Go to the directory `wrap_Python` and run

```
python setup.py build
python setup.py install
```

To test that this works run:

```
$ python testing.py
```

The module is now available under the name `catch22`. Each function contained in the module takes arrays as tuple or lists (not `Numpy`-arrays). In Python:

```
import catch22
```

## R

This assumes your have `R` installed and the package `Rcpp` is available.

Copy all `.c`- and `.h`-files from `./C_Functions` to `./wrap_R/src`. Then go to the directory `./wrap_R` and run

```
R CMD build catch22
R CMD INSTALL catch22_0.1.tar.gz
```

The module is now available under the name `catch22`. In `R`:

```
library(catch22)
```

To test, navigate to `./wrap_R` in the console and run:

```
$ Rscript testing.R
```

## Matlab

Go to the `./wrap_Matlab` directory and call `mexAll` from within Matlab. Include the folder in your Matlab path to use the package.

To test, navigate to `./wrap_Matlab` directory from within Matlab and run:

```
testing
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
