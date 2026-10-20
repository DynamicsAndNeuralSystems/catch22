#include <math.h>
#include <stdlib.h>
#include "stats.h"
#include "histcounts.h"
#include "DN_HistogramMode.h"

// hctsa's DN_HistogramMode(y,nBins,true): the mean centre of the most populated of
// nBins equal-width bins (BF_HistEdges edges, histcounts counting), using edges[nBins + 1]
// and counts[nBins] as workspace
static double histogram_mode(const double y[], const int size, const int nBins, double edges[], int counts[])
{
    for (int i = 0; i < size; i++) {
        if (isnan(y[i])) {
            return NAN;
        }
    }
    const int numBins = hist_edges(y, size, nBins, edges);
    hist_counts_edges(y, size, edges, numBins, counts);

    int maxCount = -1;
    for (int i = 0; i < numBins; i++) {
        if (counts[i] > maxCount) {
            maxCount = counts[i];
        }
    }
    double sumCentres = 0;
    int numMaxs = 0;
    for (int i = 0; i < numBins; i++) {
        if (counts[i] == maxCount) {
            sumCentres += (edges[i] + edges[i + 1]) / 2;
            numMaxs++;
        }
    }
    return sumCentres / numMaxs;
}

double DN_HistogramMode_5(const double y[], const int size)
{
    double edges[5 + 1];
    int counts[5];
    return histogram_mode(y, size, 5, edges, counts);
}

double DN_HistogramMode_10(const double y[], const int size)
{
    double edges[10 + 1];
    int counts[10];
    return histogram_mode(y, size, 10, edges, counts);
}
