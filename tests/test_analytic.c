// Analytic unit tests: feature values on constructed time series whose
// correct output follows from first principles (independent of any reference
// implementation).
//
// Build and run from the repository root:
//   gcc -O2 -IC tests/test_analytic.c $(ls C/*.c | grep -v main.c) -o test_analytic -lm
//   ./test_analytic

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "stats.h"
#include "CO_AutoCorr.h"
#include "FC_LocalSimple.h"
#include "IN_AutoMutualInfoStats.h"
#include "MD_hrv.h"
#include "PD_PeriodicityWang.h"
#include "SB_BinaryStats.h"
#include "SB_MotifThree.h"
#include "SP_Summaries.h"
#include "catch22.h"

#define PI 3.14159265358979323846

static int nFail = 0;
static int nPass = 0;

static void check(const char * what, const double got, const double expected, const double tol)
{
    if (fabs(got - expected) <= tol) {
        nPass++;
    }
    else {
        nFail++;
        printf("FAIL: %s: got %.10g, expected %.10g (tolerance %g)\n", what, got, expected, tol);
    }
}

// features are computed on the z-scored series, as in run_features
static double * zscored(const double y[], const int size)
{
    double * z = malloc(size * sizeof * z);
    catch22_zscore(y, size, z);
    return z;
}

static double * sine(const int size, const double period)
{
    double * y = malloc(size * sizeof * y);
    for (int i = 0; i < size; i++) {
        y[i] = sin(2 * PI * i / period);
    }
    return y;
}

// deterministic standard normal noise (LCG + Box-Muller), so tests are reproducible
static double * gaussian_noise(const int size)
{
    double * y = malloc(size * sizeof * y);
    unsigned long long state = 42;
    for (int i = 0; i < size; i += 2) {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        double u1 = ((state >> 11) + 1.0) / 9007199254740993.0;
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        double u2 = ((state >> 11) + 1.0) / 9007199254740993.0;
        double r = sqrt(-2 * log(u1));
        y[i] = r * cos(2 * PI * u2);
        if (i + 1 < size) {
            y[i+1] = r * sin(2 * PI * u2);
        }
    }
    return y;
}

static void test_sines(void)
{
    char what[100];
    const int periods[] = {8, 13, 20, 50};
    for (int k = 0; k < 4; k++) {
        const int P = periods[k];
        double * y = sine(2000, P);
        double * z = zscored(y, 2000);

        // the ACF of a sine peaks at its period and has its first minimum at half the period
        snprintf(what, sizeof what, "periodicity of a period-%d sine", P);
        check(what, PD_PeriodicityWang_th0_01(z, 2000), P, 0);
        snprintf(what, sizeof what, "acf_first_min of a period-%d sine", P);
        check(what, CO_FirstMin_ac(z, 2000), P / 2, 0);

        // ACF = cos(2 pi tau / P) first crosses 1/e between integer lags k and
        // k+1; acf_timescale interpolates linearly between them (the tolerance
        // allows for finite-sample error in the ACF, which grows with P)
        const double thresh = 1 / exp(1);
        const int k = (int) floor(P * acos(thresh) / (2 * PI));
        const double rk = cos(2 * PI * k / P), rk1 = cos(2 * PI * (k + 1) / P);
        snprintf(what, sizeof what, "acf_timescale of a period-%d sine", P);
        check(what, CO_FirstCrossing_acf_1e_pointOfCrossing(z, 2000), k + (rk - thresh) / (rk - rk1), 0.001 * P);

        // all power is at angular frequency 2 pi / P (within one frequency bin)
        snprintf(what, sizeof what, "median_freq of a period-%d sine", P);
        check(what, SP_Summaries_fft_wq_50(z, 2000), 2 * PI / P, 2 * PI / 2048);

        free(y);
        free(z);
    }

    // the Gaussian AMI is minimal where the ACF crosses zero: lag P/4 = 5
    double * y = sine(2000, 20);
    double * z = zscored(y, 2000);
    check("ami_timescale of a period-20 sine", IN_AutoMutualInfoStats_40_gaussian_fmmi(z, 2000), 5, 0);
    free(y);
    free(z);
}

static void test_stretches(void)
{
    // a ramp: no decreasing steps; its reverse decreases at every step; the top
    // half (including the final point) is above the mean
    double ramp[100], rampDown[100];
    for (int i = 0; i < 100; i++) {
        ramp[i] = i;
        rampDown[i] = 99 - i;
    }
    double * z = zscored(ramp, 100);
    check("stretch_decreasing of an increasing ramp", SB_BinaryStats_diff_longstretch0(z, 100), 0, 0);
    check("stretch_high of an increasing ramp (final point included)", SB_BinaryStats_mean_longstretch1(z, 100), 50, 0);
    check("forecast_error of a ramp", FC_LocalSimple_mean3_stde(z, 100), 0, 1e-12);
    free(z);
    z = zscored(rampDown, 100);
    check("stretch_decreasing of a decreasing ramp", SB_BinaryStats_diff_longstretch0(z, 100), 99, 0);
    free(z);

    // runs of 7 and 3 above the mean, not touching either end
    double runs[60] = {0};
    for (int i = 10; i < 17; i++) runs[i] = 1;
    for (int i = 30; i < 33; i++) runs[i] = 1;
    z = zscored(runs, 60);
    check("stretch_high with interior runs of 7 and 3", SB_BinaryStats_mean_longstretch1(z, 60), 7, 0);
    free(z);

    // a final run of 9 above the mean is the longest
    for (int i = 51; i < 60; i++) runs[i] = 1;
    z = zscored(runs, 60);
    check("stretch_high with a final run of 9", SB_BinaryStats_mean_longstretch1(z, 60), 9, 0);
    free(z);

    // rise, then exactly 5 decreasing steps, then rise again
    double steps[30];
    for (int i = 0; i < 10; i++) steps[i] = i;              // 0..9
    for (int i = 10; i < 15; i++) steps[i] = 18 - i;        // 8..4: 5 decreasing steps
    for (int i = 15; i < 30; i++) steps[i] = i - 10;        // 5..19
    z = zscored(steps, 30);
    check("stretch_decreasing with an interior run of 5 decreases", SB_BinaryStats_diff_longstretch0(z, 30), 5, 0);
    free(z);

    // zero differences count as non-increases (as in hctsa): 4 flat steps
    double flat[25];
    for (int i = 0; i < 10; i++) flat[i] = i;               // 0..9
    for (int i = 10; i < 14; i++) flat[i] = 9;              // 4 flat steps
    for (int i = 14; i < 25; i++) flat[i] = i - 4;          // 10..20
    z = zscored(flat, 25);
    check("stretch_decreasing with 4 flat steps", SB_BinaryStats_diff_longstretch0(z, 25), 4, 0);
    free(z);
    double tied[12] = {0, 5, 4, 3, 3, 2, 1, 0, 1, 2, 3, 4}; // 6 non-increasing steps, including a tie
    z = zscored(tied, 12);
    check("stretch_decreasing with a tie inside a decrease", SB_BinaryStats_diff_longstretch0(z, 12), 6, 0);
    free(z);
}

static void test_symbolic_and_increments(void)
{
    // alternating series: every step is far larger than 0.04 SD, and no run above the mean exceeds 1
    double alt[100];
    for (int i = 0; i < 100; i++) alt[i] = (i % 2 == 0) ? 1 : -1;
    double * z = zscored(alt, 100);
    check("high_fluctuation of an alternating series", MD_pNN_pnn40(z, 100), 1, 0);
    check("stretch_high of an alternating series", SB_BinaryStats_mean_longstretch1(z, 100), 1, 0);
    free(z);

    // repeating A,B,C: the only pairs are AB, BC and CA, with counts 100, 100 and 99
    // among the 299 successive pairs
    double abc[300];
    for (int i = 0; i < 300; i++) abc[i] = i % 3;
    z = zscored(abc, 300);
    double h = 0;
    const double counts[3] = {100, 100, 99};
    for (int k = 0; k < 3; k++) h -= counts[k] / 299 * log(counts[k] / 299);
    check("entropy_pairs of a repeating ABC sequence", SB_MotifThree_quantile_hh(z, 300), h, 1e-12);
    free(z);

    // sawtooth: slow rises with sharp drops give a negative trev; the reverse gives a positive one
    double saw[500], sawRev[500];
    for (int i = 0; i < 500; i++) {
        saw[i] = i % 10;
        sawRev[499 - i] = saw[i];
    }
    z = zscored(saw, 500);
    double trevSaw = CO_trev_1_num(z, 500);
    free(z);
    z = zscored(sawRev, 500);
    double trevRev = CO_trev_1_num(z, 500);
    free(z);
    check("trev is negative for slow rises with sharp drops", trevSaw < 0, 1, 0);
    check("trev changes sign under time reversal", trevSaw + trevRev, 0, 1e-12);
}

static void test_white_noise(void)
{
    const int N = 20000;
    double * y = gaussian_noise(N);
    double * z = zscored(y, N);

    // all nine symbol pairs equally likely: entropy log 9
    check("entropy_pairs of white noise", SB_MotifThree_quantile_hh(z, N), log(9), 0.005);
    // flat spectrum: a fifth of the power in the lowest fifth of frequencies
    check("low_freq_power of white noise", SP_Summaries_fft_area_5_1(z, N), 0.2, 0.02);
    // flat spectrum: half the power below pi/2
    check("median_freq of white noise", SP_Summaries_fft_wq_50(z, N), PI / 2, 0.05);
    // ACF drops to ~0 at lag 1, so the 1/e crossing is interpolated to 1 - 1/e
    check("acf_timescale of white noise", CO_FirstCrossing_acf_1e_pointOfCrossing(z, N), 1 - 1 / exp(1), 0.02);

    free(y);
    free(z);
}

static void test_catch22_compute(void)
{
    double out[CATCH24_N_FEATURES];
    int masked[CATCH24_N_FEATURES];

    // a constant series cannot be z-scored: masked by default, NaN without masking
    double constant[50];
    for (int i = 0; i < 50; i++) constant[i] = 3.0;
    check("constant series: return code", catch22_compute(constant, 50, 1, 1, out, masked), CATCH22_OK, 0);
    int nMasked = 0, nNaN = 0;
    for (int i = 0; i < CATCH22_N_FEATURES; i++) {
        nMasked += masked[i];
        nNaN += isnan(out[i]);
    }
    check("constant series: all 22 features masked", nMasked, 22, 0);
    check("constant series: no NaN when masked", nNaN, 0, 0);
    check("constant series: mean (catch24)", out[22], 3.0, 1e-12);
    check("constant series: std (catch24)", out[23], 0.0, 1e-12);
    check("constant series: rs_range masked to one regime", out[20], 1.0, 0);
    catch22_compute(constant, 50, 0, 0, out, masked);
    nNaN = 0;
    for (int i = 0; i < CATCH22_N_FEATURES; i++) nNaN += isnan(out[i]);
    check("constant series without masking: all NaN", nNaN, 22, 0);

    // NaN in the input is a data problem: never masked
    double withNaN[50];
    for (int i = 0; i < 50; i++) withNaN[i] = i;
    withNaN[10] = NAN;
    check("NaN input: return code", catch22_compute(withNaN, 50, 0, 1, out, masked), CATCH22_BAD_INPUT, 0);
    nNaN = 0;
    for (int i = 0; i < CATCH22_N_FEATURES; i++) nNaN += isnan(out[i]);
    check("NaN input: all features NaN even with masking", nNaN, 22, 0);

    // a short series gets 22 real numbers by default: dfa needs a series of length >= 79
    // (as in hctsa), so for length 60 it is masked (one scaling regime)
    double * y = gaussian_noise(60);
    catch22_compute(y, 60, 0, 1, out, masked);
    nNaN = 0;
    for (int i = 0; i < CATCH22_N_FEATURES; i++) nNaN += isnan(out[i]);
    check("length-60 series: no NaN", nNaN, 0, 0);
    check("length-60 series: dfa masked", masked[21], 1, 0);
    check("length-60 series: dfa masked to one regime", out[21], 1.0, 0);
    check("length-60 series: rs_range computed, not masked", masked[20], 0, 0);
    free(y);

    // embedding_dist on real series stays below its mask (the theoretical ceiling)
    y = gaussian_noise(1000);
    for (int i = 1; i < 1000; i++) y[i] += y[i-1]; // a random walk
    catch22_compute(y, 1000, 0, 1, out, masked);
    // the ceiling grows with n, so the masked value for a length-50 constant series is a
    // stricter test than the ceiling at n = 1000
    double ceiling[CATCH22_N_FEATURES];
    catch22_compute(constant, 50, 0, 1, ceiling, NULL);
    check("embedding_dist of a random walk is below the masked ceiling", out[19] < ceiling[19], 1, 0);
    free(y);
}

// catch22_compute_feature gives exactly the value (and mask flag) of catch22_compute
static void test_catch22_compute_feature(void)
{
    double constant[50], withNaN[50];
    for (int i = 0; i < 50; i++) {
        constant[i] = 3.0;
        withNaN[i] = i;
    }
    withNaN[10] = NAN;
    double * noise = gaussian_noise(1000);
    double * shortNoise = gaussian_noise(60);
    double * sinusoid = sine(500, 20);
    double oneSample[1] = {2.0};
    const double * series[] = {noise, shortNoise, sinusoid, constant, withNaN, oneSample};
    const int sizes[] = {1000, 60, 500, 50, 50, 1};
    const char * labels[] = {"noise", "length-60 noise", "sine", "constant", "NaN", "one sample"};

    for (int s = 0; s < 6; s++) {
        for (int mask_nans = 0; mask_nans <= 1; mask_nans++) {
            double all[CATCH24_N_FEATURES];
            int allMasked[CATCH24_N_FEATURES];
            const int allCode = catch22_compute(series[s], sizes[s], 1, mask_nans, all, allMasked);
            int nDiffer = 0;
            for (int i = 0; i < CATCH24_N_FEATURES; i++) {
                double one;
                int oneMasked;
                const int code = catch22_compute_feature(series[s], sizes[s], i, mask_nans, &one, &oneMasked);
                const int same = (one == all[i] || (isnan(one) && isnan(all[i])));
                if (!same || oneMasked != allMasked[i] || code != allCode) {
                    nDiffer++;
                    printf("  %s (mask_nans = %d), feature %d: %.17g vs %.17g\n", labels[s], mask_nans, i, one, all[i]);
                }
            }
            char what[100];
            snprintf(what, sizeof what, "compute_feature matches compute: %s, mask_nans = %d", labels[s], mask_nans);
            check(what, nDiffer, 0, 0);
        }
    }

    // by name: each long and short name finds its own feature (so no two features share a name)
    int nWrong = 0;
    for (int i = 0; i < CATCH24_N_FEATURES; i++) {
        nWrong += catch22_feature_index(catch22_feature_name(i)) != i;
        nWrong += catch22_feature_index(catch22_short_name(i)) != i;
    }
    check("feature_index: every long and short name finds its feature", nWrong, 0, 0);
    check("feature_index: unknown name", catch22_feature_index("not_a_feature"), -1, 0);
    check("feature_index: v0.x name (CO_f1ecac) is not a v1.0 name", catch22_feature_index("CO_f1ecac"), -1, 0);
    check("feature_index: NULL", catch22_feature_index(NULL), -1, 0);

    double out = 0.0;
    double all[CATCH24_N_FEATURES];
    catch22_compute(noise, 1000, 1, 1, all, NULL);
    catch22_compute_feature(noise, 1000, catch22_feature_index("trev"), 1, &out, NULL);
    check("compute_feature by name: trev", out, all[catch22_feature_index("CO_trev_1_num")], 0);
    check("compute_feature: unknown name gives CATCH22_BAD_INDEX",
          catch22_compute_feature(noise, 1000, catch22_feature_index("not_a_feature"), 1, &out, NULL), CATCH22_BAD_INDEX, 0);
    check("compute_feature: index 24 is out of range", catch22_compute_feature(noise, 1000, 24, 1, &out, NULL), CATCH22_BAD_INDEX, 0);
    check("compute_feature: out of range gives NaN", isnan(out), 1, 0);
    check("compute_feature: index -1 is out of range", catch22_compute_feature(noise, 1000, -1, 1, &out, NULL), CATCH22_BAD_INDEX, 0);
    free(noise);
    free(shortNoise);
    free(sinusoid);
}

// catch22_zscore, catch22_compute_feature_raw and catch22_mask_value compose to
// exactly catch22_compute_feature (value, mask flag and return code)
static void test_building_blocks(void)
{
    double constant[50], withNaN[50];
    for (int i = 0; i < 50; i++) {
        constant[i] = 3.0;
        withNaN[i] = i;
    }
    withNaN[10] = NAN;
    double * noise = gaussian_noise(1000);
    double * shortNoise = gaussian_noise(60);
    double * sinusoid = sine(500, 20);
    double oneSample[1] = {2.0};
    const double * series[] = {noise, shortNoise, sinusoid, constant, withNaN, oneSample};
    const int sizes[] = {1000, 60, 500, 50, 50, 1};
    const char * labels[] = {"noise", "length-60 noise", "sine", "constant", "NaN", "one sample"};

    for (int s = 0; s < 6; s++) {
        double * z = malloc(sizes[s] * sizeof * z);
        const int zCode = catch22_zscore(series[s], sizes[s], z);
        for (int mask_nans = 0; mask_nans <= 1; mask_nans++) {
            int nDiffer = 0;
            for (int i = 0; i < CATCH24_N_FEATURES; i++) {
                double one, composed;
                int oneMasked, composedMasked = 0;
                const int code = catch22_compute_feature(series[s], sizes[s], i, mask_nans, &one, &oneMasked);
                int composedCode = zCode;
                if (i >= CATCH22_N_FEATURES) {
                    composedCode = catch22_compute_feature_raw(series[s], sizes[s], i, &composed);
                }
                else {
                    composed = NAN;
                    if (zCode == CATCH22_OK) {
                        catch22_compute_feature_raw(z, sizes[s], i, &composed);
                        if (mask_nans && isnan(composed)) {
                            composed = catch22_mask_value(i, sizes[s]);
                            composedMasked = 1;
                        }
                    }
                }
                const int same = (one == composed || (isnan(one) && isnan(composed)));
                if (!same || oneMasked != composedMasked || code != composedCode) {
                    nDiffer++;
                    printf("  %s (mask_nans = %d), feature %d: %.17g vs %.17g\n", labels[s], mask_nans, i, composed, one);
                }
            }
            char what[100];
            snprintf(what, sizeof what, "building blocks compose to compute_feature: %s, mask_nans = %d", labels[s], mask_nans);
            check(what, nDiffer, 0, 0);
        }
        free(z);
    }

    // mask values: what a constant series (every feature masked) returns, for each length
    const int lengths[] = {1, 2, 50, 1000};
    int nWrong = 0;
    for (int k = 0; k < 4; k++) {
        double * y = malloc(lengths[k] * sizeof * y);
        for (int t = 0; t < lengths[k]; t++) y[t] = 1.5;
        for (int i = 0; i < CATCH22_N_FEATURES; i++) {
            double out;
            catch22_compute_feature(y, lengths[k], i, 1, &out, NULL);
            nWrong += out != catch22_mask_value(i, lengths[k]);
        }
        free(y);
    }
    check("mask_value: as masked on a constant series", nWrong, 0, 0);
    check("mask_value: ami_timescale for n = 50 is 25", catch22_mask_value(catch22_feature_index("ami_timescale"), 50), 25, 0);
    check("mask_value: mean is never masked", isnan(catch22_mask_value(CATCH22_N_FEATURES, 50)), 1, 0);
    check("mask_value: index out of range", isnan(catch22_mask_value(-1, 50)), 1, 0);
    check("mask_value: n = 0", isnan(catch22_mask_value(0, 0)), 1, 0);

    // raw: input checks, and constants
    double out;
    check("raw: NaN input gives CATCH22_BAD_INPUT", catch22_compute_feature_raw(withNaN, 50, 0, &out), CATCH22_BAD_INPUT, 0);
    check("raw: index 24 is out of range", catch22_compute_feature_raw(noise, 1000, 24, &out), CATCH22_BAD_INDEX, 0);
    check("raw: n = 0", catch22_compute_feature_raw(noise, 0, 0, &out), CATCH22_BAD_INPUT, 0);
    catch22_compute_feature_raw(constant, 50, catch22_feature_index("trev"), &out);
    check("raw: a catch22 feature of a constant series is NaN", isnan(out), 1, 0);
    catch22_compute_feature_raw(constant, 50, catch22_feature_index("mean"), &out);
    check("raw: mean of a constant series", out, 3.0, 0);
    double * zc = malloc(50 * sizeof * zc);
    check("zscore: a constant series is copied", catch22_zscore(constant, 50, zc) == CATCH22_OK && zc[7] == 3.0, 1, 0);
    check("zscore: NaN input gives CATCH22_BAD_INPUT", catch22_zscore(withNaN, 50, zc), CATCH22_BAD_INPUT, 0);
    free(zc);

    // where z-scoring over- or underflows, features are NaN (so masked), not computed
    // on the non-finite or constant z-scores
    double * scaled = malloc(1000 * sizeof * scaled);
    const double extremes[] = {1e200, 1e-200};
    for (int k = 0; k < 2; k++) {
        for (int t = 0; t < 1000; t++) scaled[t] = extremes[k] * noise[t];
        int nNotNaN = 0, nNotMasked = 0;
        for (int i = 0; i < CATCH22_N_FEATURES; i++) {
            double v;
            int vMasked;
            nNotNaN += catch22_compute_feature(scaled, 1000, i, 0, &v, NULL) != CATCH22_OK || !isnan(v);
            catch22_compute_feature(scaled, 1000, i, 1, &v, &vMasked);
            nNotMasked += !vMasked;
        }
        char what[100];
        snprintf(what, sizeof what, "z-scoring overflow: noise * %g gives NaN, masked", extremes[k]);
        check(what, nNotNaN + nNotMasked, 0, 0);
    }
    free(scaled);

    // raw on input far from z-scored: every feature returns without error (run under
    // the sanitizers in CI); outlier timing needs more thresholds than z-scored input can
    double * y = malloc(1000 * sizeof * y);
    const double scales[] = {1e-300, 1e10, 1e150, 1e300};
    int nNotOk = 0;
    for (int k = 0; k < 4; k++) {
        for (int offset = 0; offset <= 1; offset++) {
            for (int t = 0; t < 1000; t++) y[t] = offset * 1e15 + scales[k] * noise[t];
            for (int i = 0; i < CATCH24_N_FEATURES; i++) {
                nNotOk += catch22_compute_feature_raw(y, 1000, i, &out) != CATCH22_OK;
            }
        }
    }
    check("raw: far from z-scored input is computed without error", nNotOk, 0, 0);
    for (int t = 0; t < 1000; t++) y[t] = 1e10 * noise[t];
    catch22_compute_feature_raw(y, 1000, catch22_feature_index("outlier_timing_pos"), &out);
    check("raw: outlier_timing_pos is NaN beyond z-scored range", isnan(out), 1, 0);
    free(y);

    free(noise);
    free(shortNoise);
    free(sinusoid);
}

int main(void)
{
    test_sines();
    test_stretches();
    test_symbolic_and_increments();
    test_white_noise();
    test_catch22_compute();
    test_catch22_compute_feature();
    test_building_blocks();

    printf("%d passed, %d failed\n", nPass, nFail);
    return nFail > 0;
}
