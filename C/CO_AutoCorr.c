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
    
    // NaN check
    for(int i = 0; i < size; i++)
    {
        if(isnan(y[i]))
        {
            return NAN;
        }
    }
    
    int tau = co_firstzero(y, size, size);
    
    //printf("co_firstzero ran\n");
    
    if (tau > (double)size/10){
        tau = floor((double)size/10);
    }
    //printf("tau = %i\n", tau);
    
    double * d = malloc((size-tau) * sizeof(double));
    for(int i = 0; i < size-tau-1; i++)
    {
        
        d[i] = sqrt((y[i+1]-y[i])*(y[i+1]-y[i]) + (y[i+tau]-y[i+tau+1])*(y[i+tau]-y[i+tau+1]));
        
        //printf("d[%i]: %1.3f\n", i, d[i]);
        if (isnan(d[i])){
            free(d);
            return NAN;
        }
        
        /*
        if(i<100)
            printf("%i, y[i]=%1.3f, y[i+1]=%1.3f, y[i+tau]=%1.3f, y[i+tau+1]=%1.3f, d[i]: %1.3f\n", i, y[i], y[i+1], y[i+tau], y[i+tau+1], d[i]);
         */
    }
    
    //printf("embedding finished\n");
    
    // mean for exponential fit
    double l = mean(d, size-tau-1);
    
    // count histogram bin contents
    /*
     int * histCounts;
    double * binEdges;
    int nBins = histcounts(d, size-tau-1, -1, &histCounts, &binEdges);
     */
    
    int nBins = num_bins_auto(d, size-tau-1);
    if (nBins == 0){ // distances (almost) constant: no distribution to fit
        free(d);
        return NAN;
    }
    int * histCounts = malloc(nBins * sizeof(double));
    double * binEdges = malloc((nBins + 1) * sizeof(double));
    histcounts_preallocated(d, size-tau-1, nBins, histCounts, binEdges);
    
    //printf("histcount ran\n");
    
    // normalise to a probability density
    double * histCountsNorm = malloc(nBins * sizeof(double));
    for(int i = 0; i < nBins; i++){
        //printf("histCounts %i: %i\n", i, histCounts[i]);
        // probability density (count / (n * bin width)), so it is comparable with
        // the exponential density below (as hctsa's 'pdf' normalisation)
        histCountsNorm[i] = (double)histCounts[i]/(double)(size-tau-1)/(binEdges[i+1] - binEdges[i]);
        //printf("histCounts norm %i: %1.3f\n", i, histCountsNorm[i]);
    }
    
    /*
    for(int i = 0; i < nBins; i++){
        printf("histCounts[%i] = %i\n", i, histCounts[i]);
    }
    for(int i = 0; i < nBins; i++){
        printf("histCountsNorm[%i] = %1.3f\n", i, histCountsNorm[i]);
    }
    for(int i = 0; i < nBins+1; i++){
        printf("binEdges[%i] = %1.3f\n", i, binEdges[i]);
    }
    */
     
    
    //printf("histcounts normed\n");
    
    double * d_expfit_diff = malloc(nBins * sizeof(double));
    for(int i = 0; i < nBins; i++){
        double expf = exp(-(binEdges[i] + binEdges[i+1])*0.5/l)/l;
        if (expf < 0){
            expf = 0;
        }
        d_expfit_diff[i] = fabs(histCountsNorm[i]-expf);
        //printf("d_expfit_diff %i: %1.3f\n", i, d_expfit_diff[i]);
    }
    
    double out = mean(d_expfit_diff, nBins);
    
    //printf("out = %1.6f\n", out);
    //printf("reached free statements\n");
    
    // arrays created dynamically in function histcounts
    free(d);
    free(d_expfit_diff);
    free(binEdges);
    free(histCountsNorm);
    free(histCounts);
    
    return out;
    
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
    
    // NaN check
    for(int i = 0; i < size; i++)
    {
        if(isnan(y[i]))
        {
            return NAN;
        }
    }
    
    //const int tau = 2;
    //const int numBins = 5;
    
    // need at least one pair of values tau apart
    if (size <= tau) {
        return NAN;
    }
    
    double * y1 = malloc((size-tau) * sizeof(double));
    double * y2 = malloc((size-tau) * sizeof(double));
    
    for(int i = 0; i < size-tau; i++){
        y1[i] = y[i];
        y2[i] = y[i+tau];
    }
    
    // set bin edges
    const double maxValue = max_(y, size);
    const double minValue = min_(y, size);
    
    double binStep = (maxValue - minValue + 0.2)/5;
    //double binEdges[numBins+1] = {0};
	double binEdges[5+1] = {0};
    for(int i = 0; i < numBins+1; i++){
        binEdges[i] = minValue + binStep*i - 0.1;
        // printf("binEdges[%i] = %1.3f\n", i, binEdges[i]);
    }
    
    
    // count histogram bin contents
    int * bins1;
    bins1 = histbinassign(y1, size-tau, binEdges, numBins+1);
    
    int * bins2;
    bins2 = histbinassign(y2, size-tau, binEdges, numBins+1);
    
    /*
    // debug
    for(int i = 0; i < size-tau; i++){
        printf("bins1[%i] = %i, bins2[%i] = %i\n", i, bins1[i], i, bins2[i]);
    }
    */
    
    // joint
    double * bins12 = malloc((size-tau) * sizeof(double));
    //double binEdges12[(numBins + 1) * (numBins + 1)] = {0};
	double binEdges12[(5 + 1) * (5 + 1)] = {0};    

    for(int i = 0; i < size-tau; i++){
        bins12[i] = (bins1[i]-1)*(numBins+1) + bins2[i];
        // printf("bins12[%i] = %1.3f\n", i, bins12[i]);
    }
    
    for(int i = 0; i < (numBins+1)*(numBins+1); i++){
        binEdges12[i] = i+1;
        // printf("binEdges12[%i] = %1.3f\n", i, binEdges12[i]);
    }
    
    // fancy solution for joint histogram here
    int * jointHistLinear;
    jointHistLinear = histcount_edges(bins12, size-tau, binEdges12, (numBins + 1) * (numBins + 1));
    
    /*
    // debug
    for(int i = 0; i < (numBins+1)*(numBins+1); i++){
        printf("jointHistLinear[%i] = %i\n", i, jointHistLinear[i]);
    }
    */
    
    // transfer to 2D histogram (no last bin)
    double pij[numBins][numBins];
    int sumBins = 0;
    for(int i = 0; i < numBins; i++){
        for(int j = 0; j < numBins; j++){
            pij[j][i] = jointHistLinear[i*(numBins+1)+j];
            
            // printf("pij[%i][%i]=%1.3f\n", i, j, pij[i][j]);
            
            sumBins += pij[j][i];
        }
    }
    
    // normalise
    for(int i = 0; i < numBins; i++){
        for(int j = 0; j < numBins; j++){
            pij[j][i] /= sumBins;
        }
    }

    // marginals
    //double pi[numBins] = {0};
	double pi[5] = {0};
    //double pj[numBins] = {0};
	double pj[5] = {0};
    for(int i = 0; i < numBins; i++){
        for(int j = 0; j < numBins; j++){
            pi[i] += pij[i][j];
            pj[j] += pij[i][j];
            // printf("pij[%i][%i]=%1.3f, pi[%i]=%1.3f, pj[%i]=%1.3f\n", i, j, pij[i][j], i, pi[i], j, pj[j]);
        }
    }
    
    /*
    // debug
    for(int i = 0; i < numBins; i++){
        printf("pi[%i]=%1.3f, pj[%i]=%1.3f\n", i, pi[i], i, pj[i]);
    }
    */
    
    // mutual information
    double ami = 0;
    for(int i = 0; i < numBins; i++){
        for(int j = 0; j < numBins; j++){
            if(pij[i][j] > 0){
                //printf("pij[%i][%i]=%1.3f, pi[%i]=%1.3f, pj[%i]=%1.3f, logarg=, %1.3f, log(...)=%1.3f\n",
                //       i, j, pij[i][j], i, pi[i], j, pj[j], pij[i][j]/(pi[i]*pj[j]), log(pij[i][j]/(pi[i]*pj[j])));
                ami += pij[i][j] * log(pij[i][j]/(pj[j]*pi[i]));
            }
        }
    }

    // Panzeri-Treves (Miller-Madow) bias correction, as hctsa's CO_HistogramAMI:
    // subtract (Mxy - Mx - My + 1)/(2n) for Mxy occupied joint bins, Mx and My
    // occupied marginal bins, and n pairs
    int Mxy = 0, Mx = 0, My = 0;
    for(int i = 0; i < numBins; i++){
        if(pi[i] > 0) Mx++;
        if(pj[i] > 0) My++;
        for(int j = 0; j < numBins; j++){
            if(pij[i][j] > 0) Mxy++;
        }
    }
    ami -= (Mxy - Mx - My + 1) / (2.0 * (size - tau));
    
    free(bins1);
    free(bins2);
    free(jointHistLinear);
    
    free(y1);
    free(y2);
    free(bins12);
    
    return ami;
}













