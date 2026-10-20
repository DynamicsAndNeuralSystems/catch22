#include <math.h>
#include <stdlib.h>
#include "SB_CoarseGrain.h"
#include "helper_functions.h"

// hctsa's SB_MotifThree(y,'quantile').hh: the entropy of the pairs of successive
// letters of y's 3-letter quantile (equiprobable) alphabet
double SB_MotifThree_quantile_hh(const double y[], const int size)
{
    // NaN check
    for(int i = 0; i < size; i++)
    {
        if(isnan(y[i]))
        {
            return NAN;
        }
    }

    const int alphabet_size = 3;
    int * yt = calloc(size, sizeof * yt); // letters 1..3 (0 if unassigned, so not counted)
    sb_coarsegrain(y, size, alphabet_size, yt);

    // counts of each pair of successive letters
    int counts[3][3] = {{0}};
    for (int t = 0; t < size - 1; t++) {
        if (yt[t] > 0 && yt[t + 1] > 0) {
            counts[yt[t] - 1][yt[t + 1] - 1]++;
        }
    }
    free(yt);

    double hh = 0.0;
    int nonEmpty = 0;
    for (int i = 0; i < alphabet_size; i++) {
        double p[3];
        for (int j = 0; j < alphabet_size; j++) {
            p[j] = (double)counts[i][j] / ((double)(size) - (double)(1.0));
            if (p[j] > 0) {
                nonEmpty++;
            }
        }
        hh += f_entropy(p, alphabet_size);
    }
    // Miller-Madow bias correction (as hctsa's SB_MotifThree), with degrees of
    // freedom constrained by the (fixed, equiprobable) symbol marginals:
    // df = (non-empty cells) - 1 - wordLength*(alphabetSize - 1)
    const int df = nonEmpty - 1 - 2 * (alphabet_size - 1);
    if (df > 0) {
        hh += df / (2.0 * (size - 1));
    }
    return hh;
}
