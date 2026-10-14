#ifndef CATCH22_H
#define CATCH22_H

// Single entry point for computing the catch22 (or catch24) feature set, shared
// by every language wrapper so that they all behave identically.

#define CATCH22_N_FEATURES 22
#define CATCH24_N_FEATURES 24

// Return codes of catch22_compute
#define CATCH22_OK 0
#define CATCH22_BAD_INPUT 1 // the series contains NaN or Inf: all outputs are NaN

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

// As catch22_compute, also reporting the computation time of each feature (in
// milliseconds) in time_ms[] if it is not NULL.
extern int catch22_compute_timed(const double y[], const int n, const int catch24, const int mask_nans,
                                 double out[], int masked[], double time_ms[]);

// Number of features (22, or 24 with catch24), and their names (long names
// match hctsa's feature names; short names are for readability).
extern int catch22_n_features(const int catch24);
extern const char * catch22_feature_name(const int i);
extern const char * catch22_short_name(const int i);
// The feature's long name in catch22 v0.x, for mapping earlier outputs (the
// definitions have since changed to match hctsa, so values differ).
extern const char * catch22_legacy_name(const int i);

#endif
