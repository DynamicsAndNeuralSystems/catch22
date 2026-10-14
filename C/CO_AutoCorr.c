#if __cplusplus
#   include <complex>
typedef std::complex< double > cplx;
#else
#   include <complex.h>
#if defined(__GNUC__) || defined(__GNUG__)
typedef double complex cplx;
#elif defined(_MSC_VER)
typedef _Dcomplex cplx;
#endif
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "stats.h"
#include "fft.h"
#include "histcounts.h"

#include "helper_functions.h"

#ifndef CMPLX
#define CMPLX(x, y) ((cplx)((double)(x) + _Imaginary_I * (double)(y)))
#endif
#define pow2(x) (1 << x)

int nextpow2(int n)
{
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;
    return n;
}

/*
static void apply_conj(cplx a[], int size, int normalize)
{   
    switch(normalize) {
        case(1):
            for (int i = 0; i < size; i++) {
                a[i] = conj(a[i]) / size;
            }
            break;
        default:
            for (int i = 0; i < size; i++) {
                a[i] = conj(a[i]);
            }
            break;
    }
}
 */

void dot_multiply(cplx a[], cplx b[], int size)
{
    for (int i = 0; i < size; i++) {
        a[i] = _Cmulcc(a[i], conj(b[i]));
    }
}

double * CO_AutoCorr(const double y[], const int size, const int tau[], const int tau_size)
{
    double m, nFFT;
    m = mean(y, size);
    nFFT = nextpow2(size) << 1;

    cplx * F = malloc(nFFT * sizeof *F);
    cplx * tw = malloc(nFFT * sizeof *tw);
    for (int i = 0; i < size; i++) {
        
        #if defined(__GNUC__) || defined(__GNUG__)
                F[i] = CMPLX(y[i] - m, 0.0);
        #elif defined(_MSC_VER)
                cplx tmp = { y[i] - m, 0.0 };
                F[i] = tmp;
        #endif
        
    }
    for (int i = size; i < nFFT; i++) {
        #if defined(__GNUC__) || defined(__GNUG__)
                F[i] = CMPLX(0.0, 0.0);
        #elif defined(_MSC_VER)
                cplx tmp = { 0.0, 0.0 };
                F[i] = tmp; // CMPLX(0.0, 0.0);
        #endif
        
    }
    // size = nFFT;

    twiddles(tw, nFFT);
    fft(F, nFFT, tw);
    dot_multiply(F, F, nFFT);
    fft(F, nFFT, tw);
    cplx divisor = F[0];
    for (int i = 0; i < nFFT; i++) {
        //F[i] = F[i] / divisor;
        F[i] = _Cdivcc(F[i], divisor);
    }
    
    double * out = malloc(tau_size * sizeof(out));
    for (int i = 0; i < tau_size; i++) {
        out[i] = creal(F[tau[i]]);
    }
    free(F);
    free(tw);
    return out;
}

double * co_autocorrs(const double y[], const int size)
{
    double m, nFFT;
    m = mean(y, size);
    nFFT = nextpow2(size) << 1;
    
    cplx * F = malloc(nFFT * 2 * sizeof *F);
    cplx * tw = malloc(nFFT * 2 * sizeof *tw);
    for (int i = 0; i < size; i++) {
        
        #if defined(__GNUC__) || defined(__GNUG__)
                F[i] = CMPLX(y[i] - m, 0.0);
        #elif defined(_MSC_VER)
                cplx tmp = { y[i] - m, 0.0 };
                F[i] = tmp;
        #endif
    }
    for (int i = size; i < nFFT; i++) {
        
        #if defined(__GNUC__) || defined(__GNUG__)
            F[i] = CMPLX(0.0, 0.0);
        #elif defined(_MSC_VER)
            cplx tmp = { 0.0, 0.0 };
            F[i] = tmp;
        #endif
    }
    //size = nFFT;
    
    twiddles(tw, nFFT);
    fft(F, nFFT, tw);
    dot_multiply(F, F, nFFT);
    fft(F, nFFT, tw);
    cplx divisor = F[0];
    for (int i = 0; i < nFFT; i++) {
        F[i] = _Cdivcc(F[i], divisor); // F[i] / divisor;
    }
    
    double * out = malloc(nFFT * 2 * sizeof(out));
    for (int i = 0; i < nFFT; i++) {
        out[i] = creal(F[i]);
    }
    free(F);
    free(tw);
    return out;
}

int co_firstzero(const double y[], const int size, const int maxtau)
{
    
    //double * autocorrs = malloc(size * sizeof * autocorrs);
    //autocorrs = co_autocorrs(y, size);
    
    double * autocorrs = co_autocorrs(y, size);
    
    int zerocrossind = 0;
    while(autocorrs[zerocrossind] > 0 && zerocrossind < maxtau)
    {
        zerocrossind += 1;
    }
    
    free(autocorrs);
    return zerocrossind;
    
}

// First crossing of a threshold by the autocorrelation function, linearly
// interpolated between the two lags either side of it (as hctsa's
// CO_FirstCrossing(y,'ac',threshold,'continuous')). Returns size-1 if the ACF
// never crosses, and NAN if the ACF is undefined.
double co_firstcrossing_continuous(const double y[], const int size, const double threshold)
{
    if (size < 2) {
        return NAN;
    }
    double * autocorrs = co_autocorrs(y, size);
    if (!isfinite(autocorrs[0])) {
        free(autocorrs);
        return NAN;
    }
    double out = size - 1;
    for (int k = 1; k < size; k++) {
        if (autocorrs[k] < threshold) {
            out = (k - 1) + (autocorrs[k-1] - threshold) / (autocorrs[k-1] - autocorrs[k]);
            break;
        }
    }
    free(autocorrs);
    return out;
}

double CO_FirstCrossing_acf_1e_pointOfCrossing(const double y[], const int size)
{
    // NaN check
    for (int i = 0; i < size; i++) {
        if (!isfinite(y[i])) {
            return NAN;
        }
    }
    // first 1/e crossing of the autocorrelation function
    return co_firstcrossing_continuous(y, size, 1.0/exp(1));
}

double CO_Embed2_Basic_tau_incircle(const double y[], const int size, const double radius, const int tau)
{
    int tauIntern = 0;
    
    if(tau < 0)
    {
        tauIntern = co_firstzero(y, size, size);
    }
    else{
        tauIntern = tau;
    }
    
    double insidecount = 0;
    for(int i = 0; i < size-tauIntern; i++)
    {
        if(y[i]*y[i] + y[i+tauIntern]*y[i+tauIntern] < radius)
        {
            insidecount += 1;
        }
    }
    
    return insidecount/(size-tauIntern);
}

double CO_Embed2_Dist_tau_d_expfit_meandiff(const double y[], const int size)
{
    // hctsa's CO_Embed2_Dist(y,'tau').d_expfit_meandiff: in the 2-d embedding at
    // tau (first zero of the ACF, at most N/10), the mean absolute difference
    // between the distribution ('pdf') of successive distances, binned by hctsa's
    // 'auto' rule (BF_HistEdges), and an exponential fitted to them, at the bin
    // centres. NaN if all distances are equal.
    for (int i = 0; i < size; i++) {
        if (isnan(y[i])) {
            return NAN;
        }
    }
    int tau = co_firstzero(y, size, size);
    if (tau > (double)size / 10) {
        tau = floor((double)size / 10);
    }
    const int n = size - tau - 1; // number of successive distances
    if (n < 1) {
        return NAN;
    }
    double * d = malloc(n * sizeof * d);
    for (int i = 0; i < n; i++) {
        const double a = y[i + 1] - y[i];
        const double b = y[i + tau + 1] - y[i + tau];
        d[i] = sqrt(a * a + b * b);
        if (isnan(d[i])) {
            free(d);
            return NAN;
        }
    }
    if (max_(d, n) == min_(d, n)) { // all distances equal: no distribution to fit
        free(d);
        return NAN;
    }
    const double l = mean(d, n); // exponential fit (maximum likelihood)

    const int numBins = hist_num_bins_auto(d, n);
    double * edges = malloc((numBins + 1) * sizeof * edges);
    int * counts = malloc(numBins * sizeof * counts);
    hist_edges(d, n, numBins, edges);
    hist_counts_edges(d, n, edges, numBins, counts);

    double sumAbs = 0;
    for (int i = 0; i < numBins; i++) {
        const double pdf = counts[i] / (n * (edges[i + 1] - edges[i])); // count/(N*width)
        const double centre = (edges[i] + edges[i + 1]) / 2;
        const double expf = (1 / l) * exp(-centre / l); // exponential pdf
        sumAbs += fabs(pdf - expf);
    }
    free(d);
    free(edges);
    free(counts);
    return sumAbs / numBins;
}

double CO_FirstMin_ac(const double y[], const int size)
{
    
    // NaN check
    for(int i = 0; i < size; i++)
    {
        if(!isfinite(y[i]))
        {
            return NAN;
        }
    }
    if (size < 3) {
        return NAN;
    }
    
    double * autocorrs = co_autocorrs(y, size);
    
    int minInd = size;
    for(int i = 1; i < size-1; i++)
    {
        if(autocorrs[i] < autocorrs[i-1] && autocorrs[i] < autocorrs[i+1])
        {
            minInd = i;
            break;
        }
    }
    
    free(autocorrs);
    
    return minInd;
    
}

double CO_trev_1_num(const double y[], const int size)
{
    
    // NaN check
    for(int i = 0; i < size; i++)
    {
        if(isnan(y[i]))
        {
            return NAN;
        }
    }
    
    int tau = 1;
    
    // need at least one difference
    if (size < 2) {
        return NAN;
    }
    
    double * diffTemp = malloc((size-1) * sizeof * diffTemp);
    
    for(int i = 0; i < size-tau; i++)
    {
        diffTemp[i] = pow(y[i+1] - y[i],3);
    }
    
    double out;
    
    out = mean(diffTemp, size-tau);
    
    free(diffTemp);
    
    return out;
}

#define tau 2
#define numBins 5

double CO_HistogramAMI_even_5bin_ami2(const double y[], const int size)
{
    // hctsa's CO_HistogramAMI(y,2,'even',5): automutual information at lag 2 from
    // 5 equal-width bins spanning the series (BF_HistEdges edges, histcounts2
    // counting), with the Panzeri-Treves bias correction
    for (int i = 0; i < size; i++) {
        if (isnan(y[i])) {
            return NAN;
        }
    }
    if (size <= tau) { // need at least one pair of values tau apart
        return NAN;
    }
    double edges[numBins + 1];
    const int nb = hist_edges(y, size, numBins, edges);
    const int nSamples = size - tau;

    double pij[numBins][numBins] = {{0}};
    for (int t = 0; t < nSamples; t++) {
        const int b1 = hist_bin_index(y[t], edges, nb);
        const int b2 = hist_bin_index(y[t + tau], edges, nb);
        pij[b1][b2] += 1;
    }
    double total = 0;
    for (int j = 0; j < nb; j++) { // column-major, as MATLAB
        for (int i = 0; i < nb; i++) {
            total += pij[i][j];
        }
    }
    for (int i = 0; i < nb; i++) {
        for (int j = 0; j < nb; j++) {
            pij[i][j] /= total;
        }
    }
    double pi[numBins] = {0}, pj[numBins] = {0}; // marginals over y1 (rows) and y2 (columns)
    for (int i = 0; i < nb; i++) {
        for (int j = 0; j < nb; j++) {
            pi[i] += pij[i][j];
        }
    }
    for (int j = 0; j < nb; j++) {
        for (int i = 0; i < nb; i++) {
            pj[j] += pij[i][j];
        }
    }

    double ami = 0;
    int Mxy = 0, Mx = 0, My = 0;
    for (int j = 0; j < nb; j++) { // column-major, as MATLAB's sum over pij(r)
        for (int i = 0; i < nb; i++) {
            if (pij[i][j] > 0) {
                ami += pij[i][j] * log(pij[i][j] / pi[i] / pj[j]);
                Mxy++;
            }
        }
    }
    for (int i = 0; i < nb; i++) {
        if (pi[i] > 0) Mx++;
        if (pj[i] > 0) My++;
    }
    ami -= (Mxy - Mx - My + 1) / (2.0 * nSamples);
    return ami;
}













