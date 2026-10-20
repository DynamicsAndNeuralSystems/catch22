#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "stats.h"
#include "CO_AutoCorr.h"

/* Tile width for the loop-interchanged accumulators (keep tiles in L1). */
#define FC_TILE 512

/* ------------------------------------------------------------------ */
/* helpers                                                            */
/* ------------------------------------------------------------------ */


#define FC_RESID_FIXED(LEN)                                             \
    do {                                                                \
        for (int i = 0; i < n; i++) {                                   \
            double yest = 0.0;                                          \
            for (int j = 0; j < (LEN); j++) yest += y[i + j];           \
            res[i] = y[i + (LEN)] - yest / (double)(LEN);               \
        }                                                               \
    } while (0)

static void fc_mean_residuals(const double * restrict y, const int size,
                              const int L, double * restrict res)
{
    const int n = size - L;

    switch (L) {
    case 1: FC_RESID_FIXED(1); return;
    case 2: FC_RESID_FIXED(2); return;
    case 3: FC_RESID_FIXED(3); return;
    case 4: FC_RESID_FIXED(4); return;
    default: break;
    }

    const double dL = (double)L;
    for (int i0 = 0; i0 < n; i0 += FC_TILE) {
        const int i1 = (i0 + FC_TILE < n) ? (i0 + FC_TILE) : n;

        for (int i = i0; i < i1; i++) res[i] = 0.0;

        for (int j = 0; j < L; j++) {
            const double * restrict yj = y + j;
            for (int i = i0; i < i1; i++) res[i] += yj[i];
        }

        for (int i = i0; i < i1; i++) res[i] = y[i + L] - res[i] / dL;
    }
}

/* ------------------------------------------------------------------ */
/* mean-forecast features                                             */
/* ------------------------------------------------------------------ */

double FC_LocalSimple_mean_taurat(const double y[], const int size, const int train_length)
{
    if (size <= train_length) return NAN;

    const int n = size - train_length;
    double * res = malloc(n * sizeof *res);

    fc_mean_residuals(y, size, train_length, res);

    // ratio of (interpolated) first zero crossings of the residuals' and the
    // series' autocorrelation functions, as hctsa's MF_ResidualAnalysis taurat
    double resAC1stZ = co_firstcrossing_continuous(res, n, 0.0);
    double yAC1stZ   = co_firstcrossing_continuous(y, size, 0.0);
    double output    = (yAC1stZ == 0 || !isfinite(yAC1stZ)) ? NAN : resAC1stZ / yAC1stZ;

    free(res);
    return output;
}

double FC_LocalSimple_mean_stderr(const double y[], const int size, const int train_length)
{
    if (size <= train_length) return NAN;

    const int n = size - train_length;
    double * res = malloc(n * sizeof *res);

    fc_mean_residuals(y, size, train_length, res);

    double output = stddev(res, n);

    free(res);
    return output;
}

double FC_LocalSimple_mean3_stde(const double y[], const int size)
{
    return FC_LocalSimple_mean_stderr(y, size, 3);
}

double FC_LocalSimple_mean1_taurat(const double y[], const int size)
{
    return FC_LocalSimple_mean_taurat(y, size, 1);
}
