//
//  histcounts.c
//  C_polished
//
//  Created by Carl Henning Lubba on 19/09/2018.
//  Copyright © 2018 Carl Henning Lubba. All rights reserved.
//

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include "stats.h"
#include "histcounts.h"

int num_bins_auto(const double y[], const int size){
    
    double maxVal = max_(y, size);
    double minVal = min_(y, size);
    
    if (stddev(y, size) < 0.001){
        return 0;
    }
    
    return ceil((maxVal-minVal)/(3.5*stddev(y, size)/pow(size, 1/3.)));
    
}

int histcounts_preallocated(const double y[], const int size, int nBins, int * binCounts, double * binEdges)
{
    
    int i = 0;
    
    // check min and max of input array
    double minVal = DBL_MAX, maxVal=-DBL_MAX;
    for(int i = 0; i < size; i++)
    {
        // printf("histcountInput %i: %1.3f\n", i, y[i]);
        
        if (y[i] < minVal)
        {
            minVal = y[i];
        }
        if (y[i] > maxVal)
        {
            maxVal = y[i];
        }
    }
    
    // and derive bin width from it
    double binStep = (maxVal - minVal)/nBins;
    
    // variable to store counted occurances in
    for(i = 0; i < nBins; i++)
    {
        binCounts[i] = 0;
    }
    
    for(i = 0; i < size; i++)
    {
        
        int binInd = (y[i]-minVal)/binStep;
        if(binInd < 0)
            binInd = 0;
        if(binInd >= nBins)
            binInd = nBins-1;
        //printf("histcounts, i=%i, binInd=%i, nBins=%i\n", i, binInd, nBins);
        binCounts[binInd] += 1;
        
    }
    
    for(i = 0; i < nBins+1; i++)
    {
        binEdges[i] = i * binStep + minVal;
    }
    
    /*
     // debug
     for(i=0;i<nBins;i++)
     {
     printf("%i: count %i, edge %1.3f\n", i, binCounts[i], binEdges[i]);
     }
     */
    
    return 0;
    
}

int histcounts(const double y[], const int size, int nBins, int ** binCounts, double ** binEdges)
{

    int i = 0;
    
    // check min and max of input array
    double minVal = DBL_MAX, maxVal=-DBL_MAX;
    for(int i = 0; i < size; i++)
    {
        // printf("histcountInput %i: %1.3f\n", i, y[i]);
        
        if (y[i] < minVal)
        {
            minVal = y[i];
        }
        if (y[i] > maxVal)
        {
            maxVal = y[i];
        }
    }
    
    // if no number of bins given, choose spaces automatically
    if (nBins <= 0){
        nBins = ceil((maxVal-minVal)/(3.5*stddev(y, size)/pow(size, 1/3.)));
    }
    
    // and derive bin width from it
    double binStep = (maxVal - minVal)/nBins;
    
    // variable to store counted occurances in
    *binCounts = malloc(nBins * sizeof(int));
    for(i = 0; i < nBins; i++)
    {
        (*binCounts)[i] = 0;
    }
    
    for(i = 0; i < size; i++)
    {
        
        int binInd = (y[i]-minVal)/binStep;
        if(binInd < 0)
            binInd = 0;
        if(binInd >= nBins)
            binInd = nBins-1;
        (*binCounts)[binInd] += 1;
        
    }
    
    *binEdges = malloc((nBins+1) * sizeof(double));
    for(i = 0; i < nBins+1; i++)
    {
        (*binEdges)[i] = i * binStep + minVal;
    }
   
    /*
    // debug
    for(i=0;i<nBins;i++)
    {
        printf("%i: count %i, edge %1.3f\n", i, binCounts[i], binEdges[i]);
    }
    */
    
    return nBins;
    
}

int * histbinassign(const double y[], const int size, const double binEdges[], const int nEdges)
{
    
    
    // variable to store counted occurances in
    int * binIdentity = malloc(size * sizeof(int));
    for(int i = 0; i < size; i++)
    {
        // if not in any bin -> 0
        binIdentity[i] = 0;
        
        // go through bin edges
        for(int j = 0; j < nEdges; j++){
            if(y[i] < binEdges[j]){
                binIdentity[i] = j;
                break;
            }
        }
    }
    
    return binIdentity;
    
}

int * histcount_edges(const double y[], const int size, const double binEdges[], const int nEdges)
{
    
    
    int * histcounts = malloc(nEdges * sizeof(int));
    for(int i = 0; i < nEdges; i++){
        histcounts[i] = 0;
    }
    
    for(int i = 0; i < size; i++)
    {
        // go through bin edges
        for(int j = 0; j < nEdges; j++){
            if(y[i] <= binEdges[j]){
                histcounts[j] += 1;
                break;
            }
        }
    }
    
    return histcounts;
    
}

// ------------------------------------------------------------------------------
// Equal-width histograms as in hctsa's BF_HistEdges + histcounts
// ------------------------------------------------------------------------------

static int compare_doubles(const void * a, const void * b)
{
    const double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

// Interquartile range with MATLAB's documented percentile convention: the sorted
// values are the 100*(i - 0.5)/n percentiles, with linear interpolation between
// them and the extreme values beyond them.
double iqr_matlab(const double y[], const int size)
{
    double * s = malloc(size * sizeof * s);
    memcpy(s, y, size * sizeof * s);
    qsort(s, size, sizeof * s, compare_doubles);
    double q[2];
    const double p[2] = {0.25, 0.75};
    for (int k = 0; k < 2; k++) {
        const double r = size * p[k] - 0.5; // 0-based position
        if (r <= 0) {
            q[k] = s[0];
        }
        else if (r >= size - 1) {
            q[k] = s[size - 1];
        }
        else {
            const int i = (int)floor(r);
            const double frac = r - i;
            q[k] = (frac == 0) ? s[i] : s[i] + frac * (s[i + 1] - s[i]);
        }
    }
    free(s);
    return q[1] - q[0];
}

// hctsa's 'auto' rule: the larger of the Sturges and Freedman-Diaconis numbers
// of bins (Sturges when the IQR is 0), and between 1 and size.
int hist_num_bins_auto(const double y[], const int size)
{
    const double range = max_(y, size) - min_(y, size);
    const double maxAbs = fmax(fabs(max_(y, size)), fabs(min_(y, size)));
    const double numSturges = ceil(log2((double)size) + 1);
    double fdWidth = 0;
    const double q = iqr_matlab(y, size);
    if (range > 0 && q > 1e-12 * maxAbs) { // (an IQR at the level of rounding counts as 0)
        fdWidth = 2 * q * pow((double)size, -1.0 / 3.0);
    }
    // (in double: a rounding-level IQR gives a Freedman-Diaconis count far
    // beyond the range of an int, before it is capped at size)
    const double numFD = (fdWidth > 0) ? ceil(range / fdWidth) : numSturges;
    double numBins = fmax(numSturges, numFD);
    numBins = fmax(1, fmin(numBins, size));
    return (int)numBins;
}

// numBins + 1 equal-width bin edges spanning [min(y), max(y)], as hctsa's
// BF_HistEdges: interior edges lowered, and end edges widened, by 1e-6 of a bin
// width (so values on an ideal edge fall in the upper bin). Data constant to
// within rounding (range at most 1e-12 of the largest magnitude) get one bin of
// unit width. Returns the number of bins (1 if constant).
int hist_edges(const double y[], const int size, const int numBins, double edges[])
{
    const double lo = min_(y, size), hi = max_(y, size);
    if (hi - lo <= 1e-12 * fmax(fabs(lo), fabs(hi))) {
        edges[0] = lo - 0.5;
        edges[1] = lo + 0.5;
        return 1;
    }
    const double binWidth = (hi - lo) / numBins;
    const double tol = 1e-6 * binWidth;
    for (int i = 0; i <= numBins; i++) {
        edges[i] = (lo + i * binWidth) - tol;
    }
    edges[0] = lo - tol;
    edges[numBins] = hi + tol;
    return numBins;
}

// Index of the bin containing x, with MATLAB histcounts' convention: bin i is
// [edges[i], edges[i+1]), except the last, which also includes its right edge.
// -1 if x lies outside [edges[0], edges[numBins]].
int hist_bin_index(const double x, const double edges[], const int numBins)
{
    if (x < edges[0] || x > edges[numBins]) {
        return -1;
    }
    if (x == edges[numBins]) {
        return numBins - 1;
    }
    int lo = 0, hi = numBins; // edges[lo] <= x < edges[hi]
    while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;
        if (x >= edges[mid]) {
            lo = mid;
        }
        else {
            hi = mid;
        }
    }
    return lo;
}

// Counts in each of numBins bins (histcounts' convention, as hist_bin_index)
void hist_counts_edges(const double y[], const int size, const double edges[], const int numBins, int counts[])
{
    for (int i = 0; i < numBins; i++) {
        counts[i] = 0;
    }
    for (int i = 0; i < size; i++) {
        const int b = hist_bin_index(y[i], edges, numBins);
        if (b >= 0) {
            counts[b]++;
        }
    }
}
