#include <math.h>
#include <stdlib.h>
#include "stats.h"
#include "histcounts.h"

// hctsa's DN_HistogramMode(y,10,true): the mean centre of the most populated of
// 10 equal-width bins (BF_HistEdges edges, histcounts counting)
double DN_HistogramMode_10(const double y[], const int size)
{
    for (int i = 0; i < size; i++) {
        if (isnan(y[i])) {
            return NAN;
        }
    }
    const int nBins = 10;
    double edges[10 + 1];
    int counts[10];
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
