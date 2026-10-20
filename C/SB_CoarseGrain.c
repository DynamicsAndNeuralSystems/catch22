#include <math.h>
#include <stdlib.h>
#include "stats.h"
#include "helper_functions.h"

// Labels 1..num_groups of y's quantile (equiprobable) coarse-graining
void sb_coarsegrain(const double y[], const int size, const int num_groups, int labels[])
{
    int i, j;
    double * th = malloc((num_groups + 1) * sizeof *th);
    double * ls = malloc((num_groups + 1) * sizeof *ls);
    linspace(0, 1, num_groups + 1, ls);
    for (i = 0; i < num_groups + 1; i++) {
        th[i] = quantile(y, size, ls[i]);
    }
    th[0] -= 1;
    for (i = 0; i < num_groups; i++) {
        for (j = 0; j < size; j++) {
            if (y[j] > th[i] && y[j] <= th[i + 1]) {
                labels[j] = i + 1;
            }
        }
    }
    
    free(th);
    free(ls);
}
