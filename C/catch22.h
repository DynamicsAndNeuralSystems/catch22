#ifndef CATCH22_H
#define CATCH22_H

// Single entry point for computing the catch22 (or catch24) feature set, shared
// by every language wrapper so that they all behave identically.

#define CATCH22_N_FEATURES 22
#define CATCH24_N_FEATURES 24

// Return codes of catch22_compute
#define CATCH22_OK 0
#define CATCH22_BAD_INPUT 1 // the series contains NaN or Inf: all outputs are NaN
#define CATCH22_BAD_INDEX 2 // catch22_compute_feature: no feature with that index

// Compute the features of y[0..n-1] into out[], in the order of
// catch22_feature_name(): 22 features, or 24 if catch24 is nonzero (adding the
// mean and standard deviation of the raw series).
//
// Features are computed on the z-scored series. A feature that cannot be
// estimated (e.g., for a constant or very short series) is NaN; if mask_nans is
// nonzero (recommended for statistical learning), it is instead set to a fixed
// value at the appropriate end of the feature's range (see catch22.c), so that
// 22 real numbers are always returned. NaN/Inf in the input is a data problem and
// is never masked. If masked is not NULL, masked[i] is set to 1 for each masked
// output (0 otherwise).
extern int catch22_compute(const double y[], const int n, const int catch24, const int mask_nans,
                           double out[], int masked[]);

// Compute only feature i (0 <= i < 24, in the order of catch22_feature_name()) of
// y[0..n-1] into *out, giving exactly the value that catch22_compute gives for it
// (with the same z-scoring, masking and NaN/Inf handling), but computing only that
// feature. If masked is not NULL, *masked is set to 1 if the value was masked (0
// otherwise). Returns CATCH22_BAD_INDEX (and NaN) if i is out of range.
extern int catch22_compute_feature(const double y[], const int n, const int i, const int mask_nans,
                                   double * out, int * masked);
// As catch22_compute, also reporting the computation time of each feature (in
// milliseconds) in time_ms[] if it is not NULL.
extern int catch22_compute_timed(const double y[], const int n, const int catch24, const int mask_nans,
                                 double out[], int masked[], double time_ms[]);

// The steps of catch22_compute_feature, for wrappers that compose the pipeline
// themselves (e.g., z-scoring a series once and then computing features one at a
// time). For 0 <= i < 22, catch22_compute_feature(y, n, i, mask_nans, &out, &masked) is
//     catch22_zscore(y, n, z);
//     catch22_compute_feature_raw(z, n, i, &out);
//     if (mask_nans && isnan(out)) out = catch22_mask_value(i, n);   // masked = 1
// (with CATCH22_BAD_INPUT from catch22_zscore returned as is); mean and std
// (i = 22, 23) are catch22_compute_feature_raw(y, n, i, &out), never masked.
// These, and catch22_compute, are the library's interface: the per-feature
// functions (DN_HistogramMode_5 and so on, in their own headers) are implementation
// details that expect finite, non-constant, z-scored input and check none of it, so
// wrappers should not call them directly.

// Z-score y[0..n-1] into z[0..n-1] as catch22_compute does. A constant series
// cannot be z-scored and is copied unchanged (catch22_compute_feature_raw then gives
// NaN for it). Where the sum of squared deviations overflows, or every squared
// deviation underflows to zero, z is constant or non-finite, and the features are
// NaN too. Returns CATCH22_BAD_INPUT (z all NaN) if y contains NaN/Inf or n < 1.
extern int catch22_zscore(const double y[], const int n, double z[]);

// Feature i of y[0..n-1] as given: no z-scoring and no masking. The catch22
// features (0 <= i < 22) are defined on the z-scored series (catch22_zscore); on
// other input they are computed as is, which is not the catch22 feature (e.g.,
// thresholds in standard deviations become thresholds in the series' units), and
// some return NaN for input far from z-scored. Features 0-21 are NaN for a constant
// series. Returns CATCH22_BAD_INPUT if y contains NaN/Inf or n < 1, and
// CATCH22_BAD_INDEX if i is out of range (both with NaN).
extern int catch22_compute_feature_raw(const double y[], const int n, const int i, double * out);

// The value that feature i (0 <= i < 22) takes when masked, for a series of
// length n >= 1 (see catch22.c). NaN for mean and std (never masked), or for an
// index or length out of range.
extern double catch22_mask_value(const int i, const int n);

// Number of features (22, or 24 with catch24), and their names (long names
// follow hctsa's feature names; short names are for readability).
extern int catch22_n_features(const int catch24);
extern const char * catch22_feature_name(const int i);
extern const char * catch22_short_name(const int i);
// The index of the feature with this long or short name (e.g., "CO_trev_1_num" or
// "trev"), or -1 if there is none, e.g., to compute a feature by name:
// catch22_compute_feature(y, n, catch22_feature_index("trev"), 1, &out, NULL)
// (which returns CATCH22_BAD_INDEX for an unknown name).
extern int catch22_feature_index(const char * name);
// The feature's long name in catch22 v0.x, for mapping earlier outputs (most
// definitions have since changed to match hctsa, so values differ; a few features
// were only renamed, e.g., DN_Mean is now DN_mean).
extern const char * catch22_legacy_name(const int i);
// The hctsa call that reproduces the feature exactly, as MATLAB code on x_z (the
// z-scored series) or x (the raw series), with the output field selected by a
// trailing .field, e.g. "SP_Summaries(x_z,'fft',[],[],false).wq_50". The long name
// is hctsa's name for it where hctsa registers it; otherwise the name hctsa's naming
// convention gives this call (MD_pNN_pnn40, CO_HistogramAMI_even_5bin_ami2,
// SB_TransitionMatrix_quantile_3_ac_sumdiagcov).
extern const char * catch22_hctsa_call(const int i);

#endif
