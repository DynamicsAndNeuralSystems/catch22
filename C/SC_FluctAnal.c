#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <stdio.h>
#include <time.h>
#include "stats.h"
#include "CO_AutoCorr.h"

// Sum of squared residuals of a least-squares straight-line fit of y on x
static double line_fit_sse(const double x[], const double y[], const int n)
{
    double m = 0.0, b = 0.0;
    linreg(n, x, y, &m, &b);
    double sse = 0.0;
    for (int i = 0; i < n; i++) {
        const double e = m * x[i] + b - y[i];
        sse += e * e;
    }
    return sse;
}

static int compare_ints(const void * a, const void * b)
{
    return (*(const int *)a > *(const int *)b) - (*(const int *)a < *(const int *)b);
}

// Fluctuation analysis with a two-regime fit, as hctsa's
// SC_FluctAnal(x,2,how,50,1,lag,true) output prop_r1: the proportion of the
// (log-spaced) timescales that fall in the first of two scaling regimes.
// isDFA selects "dfa" (detrended fluctuation analysis, linear detrending) over
// "rsrangefit" (range of linearly detrended windows).
static double SC_FluctAnal_2_50_1_logi_prop_r1(const double y[], const int size, const int lag, const int isDFA)
{
    const int tauStep = 50;

    // integrate every lag-th sample
    const int Ny = (size + lag - 1) / lag;
    if (Ny < 10) {
        return NAN;
    }
    double * yCS = malloc(Ny * sizeof * yCS);
    yCS[0] = y[0];
    for (int i = 1; i < Ny; i++) {
        yCS[i] = yCS[i-1] + y[i * lag];
    }

    // up to tauStep log-spaced timescales from 5 to floor(Ny/2), rounded and de-duplicated
    const double linLow = log(5.0);
    const double linHigh = log((double)(Ny / 2));
    int * tau = malloc(tauStep * sizeof * tau);
    for (int i = 0; i < tauStep; i++) {
        tau[i] = (int)round(exp(linLow + i * (linHigh - linLow) / (tauStep - 1)));
    }
    tau[tauStep - 1] = Ny / 2; // the endpoint exactly
    qsort(tau, tauStep, sizeof * tau, compare_ints);
    int nTau = 0;
    for (int i = 0; i < tauStep; i++) {
        if (nTau == 0 || tau[i] != tau[nTau - 1]) {
            tau[nTau++] = tau[i];
        }
    }
    if (nTau < 8) { // too few timescales to analyse
        free(yCS);
        free(tau);
        return NAN;
    }

    // fluctuation at each timescale, over the floor(Ny/tau) full windows
    double * xReg = malloc(tau[nTau - 1] * sizeof * xReg);
    for (int i = 0; i < tau[nTau - 1]; i++) {
        xReg[i] = i + 1;
    }
    double * buffer = malloc(tau[nTau - 1] * sizeof * buffer);
    double * logtt = malloc(nTau * sizeof * logtt);
    double * logFF = malloc(nTau * sizeof * logFF);
    for (int i = 0; i < nTau; i++) {
        const int t = tau[i];
        const int nWindows = Ny / t;
        double F = 0.0;
        for (int j = 0; j < nWindows; j++) {
            double m = 0.0, b = 0.0;
            linreg(t, xReg, yCS + j * t, &m, &b);
            for (int k = 0; k < t; k++) {
                buffer[k] = yCS[j * t + k] - (m * (k + 1) + b);
            }
            if (isDFA) {
                for (int k = 0; k < t; k++) {
                    F += buffer[k] * buffer[k];
                }
            }
            else {
                const double range = max_(buffer, t) - min_(buffer, t);
                F += range * range;
            }
        }
        F = isDFA ? sqrt(F / (nWindows * t)) : sqrt(F / nWindows);
        logtt[i] = log((double)t);
        logFF[i] = log(F);
    }

    // two-regime fit: choose the split minimising the total squared error of the two
    // line fits (which share the split point), each regime having at least minPoints
    // points (so at least 16 timescales are needed; otherwise NaN, masked by catch22_compute)
    const int minPoints = (int)fmax(8, round(0.25 * nTau));
    double out = NAN;
    if (nTau >= 2 * minPoints) {
        const int nSplits = nTau - 2 * minPoints + 1;
        double * err = malloc(nSplits * sizeof * err);
        double minErr = INFINITY;
        for (int s = 0; s < nSplits; s++) {
            const int i = minPoints + s; // number of points in the first regime
            const double sse = line_fit_sse(logtt, logFF, i) + line_fit_sse(logtt + i - 1, logFF + i - 1, nTau - i + 1);
            err[s] = sse / nTau;
            if (err[s] < minErr) {
                minErr = err[s];
            }
        }
        // the first split within a relative 1e-9 of the best, as hctsa, so that
        // splits tied up to rounding resolve the same way on every platform
        int breakPt = minPoints;
        for (int s = 0; s < nSplits; s++) {
            if (err[s] <= minErr * (1 + 1e-9)) {
                breakPt = minPoints + s;
                break;
            }
        }
        free(err);
        out = (double)breakPt / nTau;
    }

    free(yCS);
    free(tau);
    free(xReg);
    free(buffer);
    free(logtt);
    free(logFF);
    return out;
}

double SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1(const double y[], const int size){
    return SC_FluctAnal_2_50_1_logi_prop_r1(y, size, 2, 1);
}

double SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1(const double y[], const int size){
    return SC_FluctAnal_2_50_1_logi_prop_r1(y, size, 1, 0);
}