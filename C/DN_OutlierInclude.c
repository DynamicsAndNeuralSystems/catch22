#include <math.h>
#include <stdlib.h>
#include "stats.h"

static void bit_add(int *bit, int n, int i) { for (; i <= n; i += i & -i) bit[i]++; }
/* smallest index with prefix sum >= k */
static int bit_select(const int *bit, int n, int logn, int k) {
    int pos = 0;
    for (int pw = logn; pw > 0; pw >>= 1) {
        if (pos + pw <= n && bit[pos + pw] < k) { pos += pw; k -= bit[pos]; }
    }
    return pos + 1;
}

double DN_OutlierInclude_np_001_mdrmd(const double y[], const int size, const int sign)
{
    for (int i = 0; i < size; i++) if (isnan(y[i])) return NAN;

    double inc = 0.01;
    int tot = 0;
    double *yWork = malloc(size * sizeof(double));

    int constantFlag = 1;
    for (int i = 0; i < size; i++) {
        if (y[i] != y[0]) constantFlag = 0;
        yWork[i] = sign * y[i];
        if (yWork[i] >= 0) tot += 1;
    }
    if (constantFlag) { free(yWork); return 0; }

    double maxVal = max_(yWork, size);
    if (maxVal < inc) { free(yWork); return 0; }

    int nThresh = maxVal / inc + 1;

    /* bucket each position by the HIGHEST threshold index it satisfies.
       derived from the identical `>= j*inc` comparison, never from a bare floor(). */
    int *head = malloc(nThresh * sizeof(int));
    int *nxt  = malloc(size * sizeof(int));
    int *cnt  = calloc(nThresh, sizeof(int));
    for (int j = 0; j < nThresh; j++) head[j] = -1;

    for (int i = 0; i < size; i++) {
        double v = yWork[i];
        double q = v / inc;
        int j;
        if (q >= (double)(nThresh - 1)) j = nThresh - 1;
        else if (q < 0) j = -1;
        else j = (int)q;
        while (j + 1 < nThresh && v >= (j + 1) * inc) j++;   /* correct up   */
        while (j >= 0        && v <  j * inc)        j--;    /* correct down */
        if (j >= 0) { nxt[i] = head[j]; head[j] = i; cnt[j]++; }
    }

    /* counts at each threshold = suffix sums (monotone non-increasing) */
    int *C = malloc(nThresh * sizeof(int));
    int run = 0;
    for (int j = nThresh - 1; j >= 0; j--) { run += cnt[j]; C[j] = run; }

    /* mj and fbi need only the counts.
       msDti1[j] = mean of consecutive gaps = (max-min)/(count-1) by telescoping,
       and is NaN iff count-1 == 0. */
    int mj = 0, fbi = nThresh - 1;
    for (int j = 0; j < nThresh; j++) {
        if ((C[j] - 1) * 100.0 / tot > 2) mj = j;
    }
    for (int j = nThresh - 1; j >= 0; j--) {
        if (C[j] - 1 == 0) fbi = j;   /* 0/0 -> NaN */
    }
    int trimLimit = mj < fbi ? mj : fbi;

    /* sweep thresholds downward, inserting positions as they qualify */
    int *bit = calloc(size + 1, sizeof(int));
    int logn = 1; while ((logn << 1) <= size) logn <<= 1;

    double *msDti4 = malloc((trimLimit + 1) * sizeof(double));
    int m = 0;
    double denom = (double)size / 2;

    for (int j = nThresh - 1; j >= 0; j--) {
        for (int i = head[j]; i != -1; i = nxt[i]) { bit_add(bit, size, i + 1); m++; }
        if (j > trimLimit) continue;                 /* median not needed up here */
        double med;
        if (m % 2 == 1) {
            med = (double)bit_select(bit, size, logn, m / 2 + 1);
        } else {
            double lo = (double)bit_select(bit, size, logn, m / 2);
            double hi = (double)bit_select(bit, size, logn, m / 2 + 1);
            med = (hi + lo) / (double)2.0;
        }
        msDti4[j] = med / denom - 1;
    }

    double outputScalar = median(msDti4, trimLimit + 1);

    free(yWork); free(head); free(nxt); free(cnt); free(C); free(bit); free(msDti4);
    return outputScalar;
}

double DN_OutlierInclude_p_001_mdrmd(const double y[], const int size)
{
    return DN_OutlierInclude_np_001_mdrmd(y, size, 1.0);
}

double DN_OutlierInclude_n_001_mdrmd(const double y[], const int size)
{
    return DN_OutlierInclude_np_001_mdrmd(y, size, -1.0);
}