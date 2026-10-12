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
    zscore_norm2(y, size, z);
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
        check(what, SP_Summaries_fft_wmax_50(z, 2000), 2 * PI / P, 2 * PI / 2048);

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

    // only strict decreases count: flat steps are not decreases, and they break a run
    double flat[25];
    for (int i = 0; i < 10; i++) flat[i] = i;               // 0..9
    for (int i = 10; i < 14; i++) flat[i] = 9;              // 4 flat steps
    for (int i = 14; i < 25; i++) flat[i] = i - 4;          // 10..20
    z = zscored(flat, 25);
    check("stretch_decreasing with 4 flat steps", SB_BinaryStats_diff_longstretch0(z, 25), 0, 0);
    double tied[12] = {0, 5, 4, 3, 3, 2, 1, 0, 1, 2, 3, 4}; // decreases of 2 and 3 steps, separated by a tie
    z = zscored(tied, 12);
    check("stretch_decreasing with a tie inside a decrease", SB_BinaryStats_diff_longstretch0(z, 12), 3, 0);
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
    check("median_freq of white noise", SP_Summaries_fft_wmax_50(z, N), PI / 2, 0.05);
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

int main(void)
{
    test_sines();
    test_stretches();
    test_symbolic_and_increments();
    test_white_noise();
    test_catch22_compute();

    printf("%d passed, %d failed\n", nPass, nFail);
    return nFail > 0;
}
