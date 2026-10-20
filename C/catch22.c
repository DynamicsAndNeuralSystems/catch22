// Single entry point for computing the catch22 (or catch24) feature set.

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "catch22.h"
#include "CO_AutoCorr.h"
#include "DN_HistogramMode.h"
#include "DN_Mean.h"
#include "DN_OutlierInclude.h"
#include "DN_Spread_Std.h"
#include "FC_LocalSimple.h"
#include "IN_AutoMutualInfoStats.h"
#include "MD_hrv.h"
#include "PD_PeriodicityWang.h"
#include "SB_BinaryStats.h"
#include "SB_MotifThree.h"
#include "SB_TransitionMatrix.h"
#include "SC_FluctAnal.h"
#include "SP_Summaries.h"
#include "stats.h"

// ---------------------------------------------------------------------------
// Mask values: what a feature reports, with masking on, when it cannot be
// estimated. A constant series is treated as the limit of an infinitely slow,
// perfectly predictable series, and a failure as the absence of the structure
// that the feature looks for, which places it at one end of the feature's range.
// ---------------------------------------------------------------------------

static double mask_zero(const int n) { (void) n; return 0.0; }
static double mask_one(const int n) { (void) n; return 1.0; }
// the autocorrelation never decays (the convention when the ACF never crosses)
static double mask_acf_timescale(const int n) { return n - 1; }
// no minimum of the ACF (the convention when no minimum is found)
static double mask_acf_first_min(const int n) { return n; }
// no minimum of the AMI within the lags analysed (the feature's own convention)
static double mask_ami_timescale(const int n) { return fmin(40, (n + 1) / 2); }
// perfectly predictable transitions: the largest value possible with three
// equiprobable symbols (the joint transition matrix of a permutation)
static double mask_transition_variance(const int n) { (void) n; return 1.0 / 9.0; }
// distances maximally far from exponential (all distances equal): a practical
// ceiling, far above the values of real series. The feature is at most
// 1/range(d) + 1/mean(d), with mean(d) >= sqrt(2)*sqrt((n-1)/n)/(n-2) for a
// z-scored series; the ceiling takes range(d) = sqrt(2)*0.001. (This was a strict
// bound when distances with SD < 0.001 counted as equal; hctsa now treats only
// exactly equal distances as such, so values just short of equality can exceed it.)
static double mask_embedding_dist(const int n)
{
    const double rangeTerm = 1.0 / (sqrt(2.0) * 0.001);
    const double meanTerm = (n > 2) ? (n - 2) / (sqrt(2.0) * sqrt((n - 1.0) / n)) : 0.0;
    return rangeTerm + meanTerm;
}

typedef struct {
    const char * name;      // long name (matches hctsa)
    const char * shortName;
    const char * legacyName; // long name in catch22 v0.x (mostly with an earlier definition)
    const char * hctsaCall;  // the hctsa call (on x_z, the z-scored series) that reproduces it
    double (*compute)(const double y[], const int size); // y finite and not constant (checked by the entry points)
    double (*mask)(const int n);
} catch22_feature;

static const catch22_feature FEATURES[CATCH22_N_FEATURES] = {
    {"DN_HistogramMode_5", "mode_5", "DN_HistogramMode_5",
        "DN_HistogramMode(x_z,5,true,false)",
        DN_HistogramMode_5, mask_zero},
    {"DN_HistogramMode_10", "mode_10", "DN_HistogramMode_10",
        "DN_HistogramMode(x_z,10,true,false)",
        DN_HistogramMode_10, mask_zero},
    {"DN_OutlierInclude_p_001_mdrmd", "outlier_timing_pos", "DN_OutlierInclude_p_001_mdrmd",
        "DN_OutlierInclude(x_z,'pos',0.01).mdrmd",
        DN_OutlierInclude_p_001_mdrmd, mask_zero},
    {"DN_OutlierInclude_n_001_mdrmd", "outlier_timing_neg", "DN_OutlierInclude_n_001_mdrmd",
        "DN_OutlierInclude(x_z,'neg',0.01).mdrmd",
        DN_OutlierInclude_n_001_mdrmd, mask_zero},
    {"CO_FirstCrossing_acf_1e_pointOfCrossing", "acf_timescale", "CO_f1ecac",
        "CO_FirstCrossing(x_z,'ac',1/exp(1),'both').pointOfCrossing",
        CO_FirstCrossing_acf_1e_pointOfCrossing, mask_acf_timescale},
    {"CO_FirstMin_ac", "acf_first_min", "CO_FirstMin_ac",
        "CO_FirstMin(x_z,'ac')",
        CO_FirstMin_ac, mask_acf_first_min},
    {"SP_Summaries_fft_area_5_1", "low_freq_power", "SP_Summaries_welch_rect_area_5_1",
        "SP_Summaries(x_z,'fft',[],[],false).area_5_1",
        SP_Summaries_fft_area_5_1, mask_one},
    {"SP_Summaries_fft_wq_50", "median_freq", "SP_Summaries_welch_rect_centroid",
        "SP_Summaries(x_z,'fft',[],[],false).wq_50",
        SP_Summaries_fft_wq_50, mask_zero},
    {"FC_LocalSimple_mean3_stde", "forecast_error", "FC_LocalSimple_mean3_stderr",
        "FC_LocalSimple(x_z,'mean',3).stde",
        FC_LocalSimple_mean3_stde, mask_zero},
    {"FC_LocalSimple_mean1_taurat", "whiten_timescale", "FC_LocalSimple_mean1_tauresrat",
        "FC_LocalSimple(x_z,'mean',1).taurat",
        FC_LocalSimple_mean1_taurat, mask_zero},
    {"MD_pNN_pnn40", "high_fluctuation", "MD_hrv_classic_pnn40",
        "MD_pNN(x_z).pnn40",
        MD_pNN_pnn40, mask_zero},
    {"SB_BinaryStats_mean_longstretch1", "stretch_high", "SB_BinaryStats_mean_longstretch1",
        "SB_BinaryStats(x_z,'mean').longstretch1",
        SB_BinaryStats_mean_longstretch1, mask_zero},
    {"SB_BinaryStats_diff_longstretch0", "stretch_decreasing", "SB_BinaryStats_diff_longstretch0",
        "SB_BinaryStats(x_z,'diff').longstretch0",
        SB_BinaryStats_diff_longstretch0, mask_zero},
    {"SB_MotifThree_quantile_hh", "entropy_pairs", "SB_MotifThree_quantile_hh",
        "SB_MotifThree(x_z,'quantile').hh",
        SB_MotifThree_quantile_hh, mask_zero},
    {"CO_HistogramAMI_even_5bin_ami2", "ami2", "CO_HistogramAMI_even_2_5",
        "CO_HistogramAMI(x_z,1:5,'even',5).ami2",
        CO_HistogramAMI_even_5bin_ami2, mask_zero},
    {"CO_trev_1_num", "trev", "CO_trev_1_num",
        "CO_trev(x_z,1).num",
        CO_trev_1_num, mask_zero},
    {"IN_AutoMutualInfoStats_40_gaussian_fmmi", "ami_timescale", "IN_AutoMutualInfoStats_40_gaussian_fmmi",
        "IN_AutoMutualInfoStats(x_z,40,'gaussian').fmmi",
        IN_AutoMutualInfoStats_40_gaussian_fmmi, mask_ami_timescale},
    {"SB_TransitionMatrix_quantile_3_ac_sumdiagcov", "transition_variance", "SB_TransitionMatrix_3ac_sumdiagcov",
        "SB_TransitionMatrix(x_z,'quantile',3,'ac').sumdiagcov",
        SB_TransitionMatrix_quantile_3_ac_sumdiagcov, mask_transition_variance},
    {"PD_PeriodicityWang_th0_01", "periodicity", "PD_PeriodicityWang_th0_01",
        "PD_PeriodicityWang(x_z).th0_01",
        PD_PeriodicityWang_th0_01, mask_zero},
    {"CO_Embed2_Dist_tau_d_expfit_meandiff", "embedding_dist", "CO_Embed2_Dist_tau_d_expfit_meandiff",
        "CO_Embed2_Dist(x_z,'tau').d_expfit_meandiff",
        CO_Embed2_Dist_tau_d_expfit_meandiff, mask_embedding_dist},
    {"SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1", "rs_range", "SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1",
        "SC_FluctAnal(x_z,2,'rsrangefit',50,1,[],true).prop_r1",
        SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1, mask_one},
    {"SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1", "dfa", "SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1",
        "SC_FluctAnal(x_z,2,'dfa',50,1,2,true).prop_r1",
        SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1, mask_one},
};

// catch24 adds the location and spread of the raw (not z-scored) series
static const char * CATCH24_NAMES[2] = {"DN_mean", "DN_Spread_std"};
static const char * CATCH24_SHORT_NAMES[2] = {"mean", "std"};
static const char * CATCH24_LEGACY_NAMES[2] = {"DN_Mean", "DN_Spread_Std"};
static const char * CATCH24_HCTSA_CALLS[2] = {"DN_Mean(x,'norm')", "DN_Spread(x,'std')"};

int catch22_n_features(const int catch24)
{
    return catch24 ? CATCH24_N_FEATURES : CATCH22_N_FEATURES;
}

const char * catch22_feature_name(const int i)
{
    if (i >= 0 && i < CATCH22_N_FEATURES) return FEATURES[i].name;
    if (i >= CATCH22_N_FEATURES && i < CATCH24_N_FEATURES) return CATCH24_NAMES[i - CATCH22_N_FEATURES];
    return NULL;
}

const char * catch22_short_name(const int i)
{
    if (i >= 0 && i < CATCH22_N_FEATURES) return FEATURES[i].shortName;
    if (i >= CATCH22_N_FEATURES && i < CATCH24_N_FEATURES) return CATCH24_SHORT_NAMES[i - CATCH22_N_FEATURES];
    return NULL;
}

int catch22_feature_index(const char * name)
{
    if (name == NULL) return -1;
    for (int i = 0; i < CATCH24_N_FEATURES; i++) {
        if (strcmp(name, catch22_feature_name(i)) == 0 || strcmp(name, catch22_short_name(i)) == 0) {
            return i;
        }
    }
    return -1;
}

const char * catch22_legacy_name(const int i)
{
    if (i >= 0 && i < CATCH22_N_FEATURES) return FEATURES[i].legacyName;
    if (i >= CATCH22_N_FEATURES && i < CATCH24_N_FEATURES) return CATCH24_LEGACY_NAMES[i - CATCH22_N_FEATURES];
    return NULL;
}

const char * catch22_hctsa_call(const int i)
{
    if (i >= 0 && i < CATCH22_N_FEATURES) return FEATURES[i].hctsaCall;
    if (i >= CATCH22_N_FEATURES && i < CATCH24_N_FEATURES) return CATCH24_HCTSA_CALLS[i - CATCH22_N_FEATURES];
    return NULL;
}

// NaN/Inf in the input (or no input) is a data problem, never masked. Also reports
// whether the series is constant (so cannot be z-scored).
static int valid_input(const double y[], const int n, int * isConstant)
{
    *isConstant = 1;
    for (int i = 0; i < n; i++) {
        if (!isfinite(y[i])) {
            return 0;
        }
        if (y[i] != y[0]) {
            *isConstant = 0;
        }
    }
    return n >= 1;
}

// Feature i of a series that passed valid_input, the one place features are
// computed: the catch22 features are NaN for a constant series, which they cannot be
// computed on
static double feature_value(const double y[], const int n, const int i, const int isConstant)
{
    if (i == CATCH22_N_FEATURES) return DN_mean(y, n);
    if (i == CATCH22_N_FEATURES + 1) return (n > 1) ? DN_Spread_std(y, n) : NAN;
    return isConstant ? NAN : FEATURES[i].compute(y, n);
}

int catch22_zscore(const double y[], const int n, double z[])
{
    int isConstant;
    if (!valid_input(y, n, &isConstant)) {
        for (int t = 0; t < n; t++) z[t] = NAN;
        return CATCH22_BAD_INPUT;
    }
    if (isConstant) {
        if (z != y) memcpy(z, y, n * sizeof * z); // in place: nothing to do
    }
    else {
        zscore_norm2(y, n, z);
    }
    return CATCH22_OK;
}

int catch22_compute_feature_raw(const double y[], const int n, const int i, double * out)
{
    *out = NAN;
    if (i < 0 || i >= CATCH24_N_FEATURES) {
        return CATCH22_BAD_INDEX;
    }
    int isConstant;
    if (!valid_input(y, n, &isConstant)) {
        return CATCH22_BAD_INPUT;
    }
    *out = feature_value(y, n, i, isConstant);
    return CATCH22_OK;
}

double catch22_mask_value(const int i, const int n)
{
    return (i >= 0 && i < CATCH22_N_FEATURES && n >= 1) ? FEATURES[i].mask(n) : NAN;
}

int catch22_compute(const double y[], const int n, const int catch24, const int mask_nans,
                    double out[], int masked[])
{
    return catch22_compute_timed(y, n, catch24, mask_nans, out, masked, NULL);
}

// The building blocks composed (catch22.h). Where z-scoring overflows or underflows,
// z is constant or non-finite and the feature NaN.
int catch22_compute_feature(const double y[], const int n, const int i, const int mask_nans,
                            double * out, int * masked)
{
    *out = NAN;
    if (masked) *masked = 0;
    if (i < 0 || i >= CATCH24_N_FEATURES) {
        return CATCH22_BAD_INDEX;
    }
    if (i >= CATCH22_N_FEATURES) { // mean and std: of the series itself, never masked
        return catch22_compute_feature_raw(y, n, i, out);
    }
    if (n < 1) {
        return CATCH22_BAD_INPUT;
    }

    double * z = malloc(n * sizeof * z);
    const int status = catch22_zscore(y, n, z);
    if (status == CATCH22_OK) {
        catch22_compute_feature_raw(z, n, i, out);
    }
    free(z);
    if (status == CATCH22_OK && mask_nans && isnan(*out)) {
        *out = catch22_mask_value(i, n);
        if (masked) *masked = 1;
    }
    return status;
}

// As catch22_compute_feature for each feature, z-scoring and checking the z-scores once
int catch22_compute_timed(const double y[], const int n, const int catch24, const int mask_nans,
                          double out[], int masked[], double time_ms[])
{
    const int nOut = catch22_n_features(catch24);
    for (int i = 0; i < nOut; i++) {
        out[i] = NAN;
        if (masked) masked[i] = 0;
        if (time_ms) time_ms[i] = 0.0;
    }
    if (n < 1) {
        return CATCH22_BAD_INPUT;
    }

    double * z = malloc(n * sizeof * z);
    if (catch22_zscore(y, n, z) != CATCH22_OK) {
        free(z);
        return CATCH22_BAD_INPUT;
    }
    if (catch24) {
        catch22_compute_feature_raw(y, n, CATCH22_N_FEATURES, &out[CATCH22_N_FEATURES]);
        catch22_compute_feature_raw(y, n, CATCH22_N_FEATURES + 1, &out[CATCH22_N_FEATURES + 1]);
    }
    int zConstant;
    const int zValid = valid_input(z, n, &zConstant);
    for (int i = 0; i < CATCH22_N_FEATURES; i++) {
        const clock_t begin = clock();
        out[i] = zValid ? feature_value(z, n, i, zConstant) : NAN;
        if (time_ms) time_ms[i] = (double)(clock() - begin) * 1000 / CLOCKS_PER_SEC;
    }
    free(z);

    if (mask_nans) {
        for (int i = 0; i < CATCH22_N_FEATURES; i++) {
            if (isnan(out[i])) {
                out[i] = catch22_mask_value(i, n);
                if (masked) masked[i] = 1;
            }
        }
    }
    return CATCH22_OK;
}
