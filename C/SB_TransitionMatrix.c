//
//  SB_TransitionMatrix.c
//
//
//  Created by Carl Henning Lubba on 23/09/2018.
//

#include <math.h>
#include <stdlib.h>
#include "SB_TransitionMatrix.h"
#include "CO_AutoCorr.h"
#include "SB_CoarseGrain.h"
#include "stats.h"

// Zeroth-order modified Bessel function of the first kind, I0(x), from its
// power series sum_k ((x/2)^k / k!)^2. Converges quickly for the small
// arguments used here (|x| <= 5).
static double bessel_i0(const double x)
{
    const double halfX2 = (x / 2) * (x / 2);
    double term = 1.0;
    double sum = 1.0;
    for (int k = 1; k < 500; k++) {
        term *= halfX2 / ((double)k * k);
        sum += term;
        if (term < sum * 1e-17) {
            break;
        }
    }
    return sum;
}

// Anti-aliased downsampling of y (length size) by the integer factor q > 1.
// A windowed-sinc low-pass filter h(m), m = -K..K with K = 10q, is formed from
// the ideal response (1/q) sinc(m/q) times a Kaiser window (beta = 5), then
// scaled to unit DC gain. Output sample i is the zero-phase filtered value
// centred on input sample i*q, with zeros assumed beyond the ends of y:
//   out[i] = sum_j y[j] h(i*q - j),  i = 0..ceil(size/q)-1.
// Returns a newly allocated array (caller frees) and sets *nOut, or returns
// NULL on allocation failure.
static double * resample_down(const double y[], const int size, const int q, int * nOut)
{
    const double beta = 5.0;
    const int K = 10 * q;
    const int L = 2 * K + 1;
    const int M = (size + q - 1) / q; // ceil(size/q)

    double * h = malloc(L * sizeof(double));
    double * out = malloc(M * sizeof(double));
    if (h == NULL || out == NULL) {
        free(h);
        free(out);
        return NULL;
    }

    // windowed ideal low-pass filter, cutoff at the new Nyquist frequency
    const double i0Beta = bessel_i0(beta);
    double hSum = 0;
    for (int m = -K; m <= K; m++) {
        double u = (double)m / q;
        double g = (m == 0) ? 1.0 / q : sin(M_PI * u) / (M_PI * u) / q;
        double r = (double)m / K;
        double w = bessel_i0(beta * sqrt(1 - r * r)) / i0Beta;
        h[m + K] = g * w;
        hSum += h[m + K];
    }

    // unit DC gain
    for (int m = 0; m < L; m++) {
        h[m] /= hSum;
    }

    // delay-compensated filtering and decimation: only taps that land on
    // input samples inside the series contribute
    for (int i = 0; i < M; i++) {
        int centre = i * q;
        int jLo = centre - K < 0 ? 0 : centre - K;
        int jHi = centre + K > size - 1 ? size - 1 : centre + K;
        double acc = 0;
        for (int j = jLo; j <= jHi; j++) {
            acc += y[j] * h[centre - j + K];
        }
        out[i] = acc;
    }

    free(h);
    *nOut = M;
    return out;
}

double SB_TransitionMatrix_quantile_3_ac_sumdiagcov(const double y[], const int size)
{

    // NaN and const check
    int constant = 1;
    for(int i = 0; i < size; i++)
    {
        if(!isfinite(y[i]))
        {
            return NAN;
        }
        if(y[i] != y[0]){
            constant = 0;
        }
    }
    if (constant){
        return NAN;
    }

    const int numGroups = 3;

    // tau: first lag with a strictly negative autocorrelation (size-1 if
    // there is none), capped at floor(size/50)
    double * autocorrs = co_autocorrs(y, size);
    int tau = size - 1;
    for(int k = 1; k < size; k++){
        if(autocorrs[k] < 0){
            tau = k;
            break;
        }
    }
    free(autocorrs);
    if(tau > size/50){
        tau = size/50;
    }

    // downsample by tau with anti-aliasing (only when tau > 1)
    int nDown = size;
    double * yDown = NULL;
    if(tau > 1){
        yDown = resample_down(y, size, tau, &nDown);
        if(yDown == NULL){
            return NAN;
        }
    }
    else{
        yDown = malloc(size * sizeof(double));
        for(int i = 0; i < size; i++){
            yDown[i] = y[i];
        }
    }

    if(nDown < 2){
        free(yDown);
        return NAN;
    }

    // transfer to alphabet
    int * yCG = malloc(nDown * sizeof(int));
    sb_coarsegrain(yDown, nDown, "quantile", numGroups, yCG);

    // joint transition probabilities between consecutive symbols
    double T[3][3];
    for(int i = 0; i < numGroups; i++){
        for(int j = 0; j < numGroups; j++){
            T[i][j] = 0;
        }
    }
    for(int j = 0; j < nDown-1; j++){
        T[yCG[j]-1][yCG[j+1]-1] += 1;
    }
    for(int i = 0; i < numGroups; i++){
        for(int j = 0; j < numGroups; j++){
            T[i][j] /= (nDown-1);
        }
    }

    // trace of the covariance matrix between columns of T (rows are
    // observations): the sum of the three column variances
    double sumdiagcov = 0;
    for(int j = 0; j < numGroups; j++){
        double column[3];
        for(int i = 0; i < numGroups; i++){
            column[i] = T[i][j];
        }
        sumdiagcov += cov(column, column, numGroups);
    }

    free(yDown);
    free(yCG);

    return sumdiagcov;

}
