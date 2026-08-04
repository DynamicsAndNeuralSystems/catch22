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

static int fc_has_nan(const double * restrict y, const int size)
{
    for (int i = 0; i < size; i++) {
        if (isnan(y[i])) return 1;
    }
    return 0;
}

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

double fc_local_simple(const double y[], const int size, const int train_length)
{
    (void)train_length;

    double m = 0.0;
    for (int i = 1; i < size; i++) {
        m += fabs(y[i] - y[i - 1]);
    }
    m /= (size - 1);
    return m;
}

/* ------------------------------------------------------------------ */
/* mean-forecast features                                             */
/* ------------------------------------------------------------------ */

double FC_LocalSimple_mean_tauresrat(const double y[], const int size, const int train_length)
{
    if (fc_has_nan(y, size)) return NAN;
    if (size <= train_length) return NAN;

    const int n = size - train_length;
    double * res = malloc(n * sizeof *res);

    fc_mean_residuals(y, size, train_length, res);

    double resAC1stZ = co_firstzero(res, n, n);
    double yAC1stZ   = co_firstzero(y, size, size);
    double output    = resAC1stZ / yAC1stZ;

    free(res);
    return output;
}

double FC_LocalSimple_mean_stderr(const double y[], const int size, const int train_length)
{
    if (fc_has_nan(y, size)) return NAN;
    if (size <= train_length) return NAN;

    const int n = size - train_length;
    double * res = malloc(n * sizeof *res);

    fc_mean_residuals(y, size, train_length, res);

    double output = stddev(res, n);

    free(res);
    return output;
}

double FC_LocalSimple_mean3_stderr(const double y[], const int size)
{
    return FC_LocalSimple_mean_stderr(y, size, 3);
}

double FC_LocalSimple_mean1_tauresrat(const double y[], const int size)
{
    return FC_LocalSimple_mean_tauresrat(y, size, 1);
}

/* NB: no NaN guard here. */
double FC_LocalSimple_mean_taures(const double y[], const int size, const int train_length)
{
    if (size <= train_length) return NAN;

    const int n = size - train_length;
    double * res = malloc(n * sizeof *res);

    fc_mean_residuals(y, size, train_length, res);

    int output = co_firstzero(res, n, n);

    free(res);
    return output;
}

double FC_LocalSimple_lfit_taures(const double y[], const int size)
{
    const int L = co_firstzero(y, size, size);

    if (size <= L) return NAN;

    const int n = size - L;
    double * res = malloc(n * sizeof *res);

    double sumx = 0.0, sumx2 = 0.0;
    for (int j = 0; j < L; j++) {
        const double xj = (double)(j + 1);
        sumx  += xj;
        sumx2 += xj * xj;
    }

    const double denom = ((double)L * sumx2 - sumx * sumx);
    const double xpred = (double)(L + 1);

    if (denom == 0.0) {
        /* linreg's singular branch: m = b = 0 for every window. */
        for (int i = 0; i < n; i++) res[i] = y[i + L] - (0.0 * xpred + 0.0);
    } else {
        double * sumy  = malloc(FC_TILE * sizeof *sumy);
        double * sumxy = malloc(FC_TILE * sizeof *sumxy);

        for (int i0 = 0; i0 < n; i0 += FC_TILE) {
            const int i1  = (i0 + FC_TILE < n) ? (i0 + FC_TILE) : n;
            const int cnt = i1 - i0;

            for (int k = 0; k < cnt; k++) { sumy[k] = 0.0; sumxy[k] = 0.0; }

            for (int j = 0; j < L; j++) {
                const double xj = (double)(j + 1);
                const double * restrict yj = y + i0 + j;
                for (int k = 0; k < cnt; k++) {
                    const double v = yj[k];
                    sumxy[k] += xj * v;
                    sumy[k]  += v;
                }
            }

            for (int k = 0; k < cnt; k++) {
                const double m = ((double)L * sumxy[k] - sumx * sumy[k]) / denom;
                const double b = (sumy[k] * sumx2 - sumx * sumxy[k]) / denom;
                res[i0 + k] = y[i0 + k + L] - (m * xpred + b);
            }
        }

        free(sumxy);
        free(sumy);
    }

    int output = co_firstzero(res, n, n);

    free(res);
    return output;
}