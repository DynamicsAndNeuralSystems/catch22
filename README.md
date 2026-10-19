<p align="center"><img src="img/catch22_logo_square.png" alt="catch22 logo" height="200"/></p>

<h1 align="center"><em>catch22</em>: CAnonical Time-series CHaracteristics</h1>

<p align="center">
  <a href="https://doi.org/10.1007/s10618-019-00647-x"><img src="https://img.shields.io/badge/paper-Data%20Min%20Knowl%20Disc%20(2019)-informational" alt="Paper"/></a>
  <a href="https://doi.org/10.5281/zenodo.4431165"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.4431165.svg" alt="DOI"/></a>
  <a href="https://www.gnu.org/licenses/gpl-3.0"><img src="https://img.shields.io/badge/license-GPLv3-blue.svg" alt="License: GPL v3"/></a>
  <a href="https://time-series-features.gitbook.io/catch22/"><img src="https://img.shields.io/badge/docs-GitBook-blue" alt="Documentation"/></a>
  <a href="https://github.com/DynamicsAndNeuralSystems/catch22/actions/workflows/ci.yml"><img src="https://github.com/DynamicsAndNeuralSystems/catch22/actions/workflows/ci.yml/badge.svg" alt="CI"/></a>
  <br>
  <a href="https://pepy.tech/projects/pycatch22"><img src="https://static.pepy.tech/personalized-badge/pycatch22?period=total&units=INTERNATIONAL_SYSTEM&left_color=GREY&right_color=BLUE&left_text=PyPI%20downloads" alt="pycatch22 PyPI downloads"/></a>
  <a href="https://anaconda.org/conda-forge/pycatch22"><img src="https://img.shields.io/conda/dn/conda-forge/pycatch22?label=conda%20downloads&color=blue" alt="pycatch22 conda downloads"/></a>
  <a href="https://cran.r-project.org/package=Rcatch22"><img src="https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fcranlogs.r-pkg.org%2Fdownloads%2Ftotal%2F2015-01-01%3Alast-day%2FRcatch22&query=%24%5B0%5D.downloads&label=CRAN%20downloads&color=blue&cacheSeconds=86400" alt="Rcatch22 CRAN downloads"/></a>
  <a href="https://juliapkgstats.com/pkg/Catch22"><img src="https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FCatch22&query=total_requests&suffix=%2Fmonth&label=Julia%20downloads&color=blue&cacheSeconds=86400" alt="Catch22.jl downloads per month"/></a>
</p>

_catch22_ is a set of 22 time-series features, implemented in C and available in Python, R, Julia and MATLAB.
The features are a high-performing, minimally redundant subset of the more than 7000 features in [_hctsa_](https://github.com/benfulcher/hctsa), selected for their classification performance across 93 real-world time-series classification problems ([Lubba et al., 2019](https://doi.org/10.1007/s10618-019-00647-x)).
They capture properties of a time series's dynamics, such as its autocorrelation, distribution of values, extreme events, scaling and predictability, and are fast to compute.

📖 **Full documentation:** [time-series-features.gitbook.io/catch22](https://time-series-features.gitbook.io/catch22/), covering installation, what each feature measures, and publications that use _catch22_.

## Installation

| Language | Package | Latest | Install |
|---|---|---|---|
| Python | [pycatch22](https://github.com/DynamicsAndNeuralSystems/pycatch22) | [![PyPI](https://img.shields.io/pypi/v/pycatch22?label=PyPI)](https://pypi.org/project/pycatch22/) [![conda-forge](https://img.shields.io/conda/vn/conda-forge/pycatch22?label=conda-forge)](https://anaconda.org/conda-forge/pycatch22) | `pip install pycatch22` or `conda install -c conda-forge pycatch22` |
| R | [Rcatch22](https://github.com/hendersontrent/Rcatch22) | [![CRAN](https://img.shields.io/cran/v/Rcatch22?label=CRAN)](https://cran.r-project.org/package=Rcatch22) | `install.packages("Rcatch22")` |
| Julia | [Catch22.jl](https://github.com/brendanjohnharris/Catch22.jl) | [![Julia](https://img.shields.io/github/v/tag/brendanjohnharris/Catch22.jl?label=Julia)](https://github.com/brendanjohnharris/Catch22.jl) | `] add Catch22` |
| MATLAB | [`wrap_Matlab/`](wrap_Matlab) in this repository | [![release](https://img.shields.io/github/v/release/DynamicsAndNeuralSystems/catch22?label=release)](https://github.com/DynamicsAndNeuralSystems/catch22/releases) | compile with `mexAll` ([instructions](https://time-series-features.gitbook.io/catch22/language-specific-docs/matlab)) |
| C | [`C/`](C) in this repository | [![release](https://img.shields.io/github/v/release/DynamicsAndNeuralSystems/catch22?label=release)](https://github.com/DynamicsAndNeuralSystems/catch22/releases) | `gcc -O2 C/*.c -o run_features -lm` |

## Quick start

In Python:

```python
import pycatch22

x = [0.1, 0.4, 0.3, 0.9, 0.5, 0.2, 0.7, 0.6, 0.8, 0.3]  # your time series
results = pycatch22.catch22_all(x, catch24=True)  # a dict of 'names' and 'values'
```

From the command line, with one value per line in `series.txt`:

```bash
./run_features series.txt features.txt  # asks whether to compute catch22 or catch24
```

Each line of the output is a feature value, its name, and its computation time in milliseconds.

See the [documentation](https://time-series-features.gitbook.io/catch22/) for usage in each language.

## The features

<details>
<summary>The 22 features (plus the 2 extra catch24 features)</summary>

| Short name | Feature name | What it measures |
|---|---|---|
| `mode_5` | `DN_HistogramMode_5` | Mode of the z-scored distribution (5-bin histogram) |
| `mode_10` | `DN_HistogramMode_10` | Mode of the z-scored distribution (10-bin histogram) |
| `outlier_timing_pos` | `DN_OutlierInclude_p_001_mdrmd` | Timing of extreme positive values |
| `outlier_timing_neg` | `DN_OutlierInclude_n_001_mdrmd` | Timing of extreme negative values |
| `acf_timescale` | `CO_f1ecac` | First 1/_e_ crossing of the autocorrelation function |
| `acf_first_min` | `CO_FirstMin_ac` | First minimum of the autocorrelation function |
| `low_freq_power` | `SP_Summaries_welch_rect_area_5_1` | Power in the lowest fifth of frequencies |
| `centroid_freq` | `SP_Summaries_welch_rect_centroid` | Centroid of the power spectrum |
| `forecast_error` | `FC_LocalSimple_mean3_stderr` | Error of a rolling 3-sample mean forecast |
| `whiten_timescale` | `FC_LocalSimple_mean1_tauresrat` | Change in autocorrelation timescale after differencing |
| `high_fluctuation` | `MD_hrv_classic_pnn40` | Proportion of large successive differences |
| `stretch_high` | `SB_BinaryStats_mean_longstretch1` | Longest stretch of above-mean values |
| `stretch_decreasing` | `SB_BinaryStats_diff_longstretch0` | Longest stretch of non-increasing steps |
| `entropy_pairs` | `SB_MotifThree_quantile_hh` | Entropy of successive pairs in a 3-letter symbolization |
| `ami2` | `CO_HistogramAMI_even_2_5` | Automutual information at lag 2 (5-bin histogram) |
| `trev` | `CO_trev_1_num` | Time-reversibility statistic |
| `ami_timescale` | `IN_AutoMutualInfoStats_40_gaussian_fmmi` | First minimum of the automutual information function |
| `transition_variance` | `SB_TransitionMatrix_3ac_sumdiagcov` | Variance in the transition matrix of a 3-letter symbolization |
| `periodicity` | `PD_PeriodicityWang_th0_01` | Periodicity measure of Wang et al. (2007) |
| `embedding_dist` | `CO_Embed2_Dist_tau_d_expfit_meandiff` | Exponential fit to distances in a 2-d embedding |
| `rs_range` | `SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1` | Rescaled-range fluctuation analysis |
| `dfa` | `SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1` | Detrended fluctuation analysis |
| `mean` | `DN_Mean` | Mean (catch24 only) |
| `std` | `DN_Spread_Std` | Standard deviation (catch24 only) |

The list is also in [`featureList.txt`](featureList.txt).
The [documentation](https://time-series-features.gitbook.io/catch22/information-about-catch22/feature-descriptions) describes each feature in detail.

</details>

## Usage notes

- **Location and spread.**
  Each time series is _z_-scored before the 22 features are computed, so they measure its dynamics and are insensitive to its mean and variance.
  If the mean and spread matter for your application, compute _catch24_ (e.g., `catch24=True` in Python, `catch24 = TRUE` in R), which adds the mean and standard deviation.
- **Sampling.**
  A time series is an ordered sequence of values without time stamps, assumed to be evenly sampled.

## Versions

Feature values can change between versions, so **report the version of _catch22_ you used** so that your results can be reproduced.

- **v0.6.0:** four features were corrected to match their MATLAB (_hctsa_) definitions.
  `PD_PeriodicityWang_th0_01` and `IN_AutoMutualInfoStats_40_gaussian_fmmi` return the lag itself (previously the lag minus one).
  `SB_BinaryStats_mean_longstretch1` and `SB_BinaryStats_diff_longstretch0` return the exact run length (previously one more for runs not at the start of the series), and zero differences count as non-increases in `SB_BinaryStats_diff_longstretch0`.
- **v0.3:** `CO_f1ecac` changed from an integer to a linearly interpolated real value.

See the [releases](https://github.com/DynamicsAndNeuralSystems/catch22/releases) for details.

## Citing _catch22_

If you use _catch22_, please cite the open-access paper:

> C. H. Lubba, S. S. Sethi, P. Knaute, S. R. Schultz, B. D. Fulcher, N. S. Jones.
> [_catch22_: CAnonical Time-series CHaracteristics](https://doi.org/10.1007/s10618-019-00647-x).
> _Data Mining and Knowledge Discovery_ **33**, 1821–1852 (2019).

<details>
<summary>BibTeX</summary>

```bibtex
@article{lubba2019catch22,
  title   = {catch22: {CA}nonical {T}ime-series {CH}aracteristics},
  author  = {Lubba, Carl H. and Sethi, Sarab S. and Knaute, Philip and Schultz, Simon R. and Fulcher, Ben D. and Jones, Nick S.},
  journal = {Data Mining and Knowledge Discovery},
  volume  = {33},
  number  = {6},
  pages   = {1821--1852},
  year    = {2019},
  doi     = {10.1007/s10618-019-00647-x}
}
```

</details>

## License

_catch22_ is licensed under the [GNU GPL v3](http://www.gnu.org/licenses/gpl-3.0.html) (or later).
