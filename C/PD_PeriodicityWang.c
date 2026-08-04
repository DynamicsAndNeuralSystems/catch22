//
//  PD_Periodicity.c
//
//  Original by Carl Henning Lubba, 28/09/2018.
//

#include <stdlib.h>
#include <math.h>

#include "PD_PeriodicityWang.h"
#include "splinefit.h"
#include "stats.h"

#ifndef ACF_B
#define ACF_B 4
#endif

/* Fill acf[tau-1] for tau in [tLo, tHi]. */
static void PD_acf_block(const double * restrict y, const int n, const int tLo,
                         const int tHi, double * restrict acf)
{
    int t = tLo;

    for (; t + (ACF_B-1) <= tHi; t += ACF_B) {
        const int tmax = t + (ACF_B-1);

        /* i < iEnd is in range for every tau in the block */
        int iEnd = n - tmax;
        if (iEnd < 0) iEnd = 0;

        double a[ACF_B];
        for (int k = 0; k < ACF_B; k++) a[k] = 0.0;

        for (int i = 0; i < iEnd; i++) {
            const double v = y[i];
            const double * restrict w = y + i + t;
            for (int k = 0; k < ACF_B; k++) a[k] += v * w[k];
        }

        /* ragged ends: shorter taus have a few more terms, appended in order */
        for (int k = 0; k < ACF_B; k++) {
            const int tau = t + k;
            const int m   = n - tau;
            for (int i = iEnd; i < m; i++) a[k] += y[i] * y[i+tau];
            acf[tau-1] = a[k] / m;
        }
    }

    for (; t <= tHi; t++) {
        const int m = n - t;
        double acc = 0.0;
        for (int i = 0; i < m; i++) acc += y[i] * y[i+t];
        acf[t-1] = acc / m;
    }
}

int PD_PeriodicityWang_th0_01(const double * y, const int size){

    // NaN check
    for(int i = 0; i < size; i++)
    {
        if(isnan(y[i]))
        {
            return 0;
        }
    }

    const double th = 0.01;

    double * ySpline = malloc(size * sizeof(double));

    // fit a spline with 3 nodes to the data
    splinefit(y, size, ySpline);

    // subtract spline from data to remove trend
    double * ySub = malloc(size * sizeof(double));
    for(int i = 0; i < size; i++){
        ySub[i] = y[i] - ySpline[i];
    }
    free(ySpline);

    const int acmax = (int)ceil((double)size/3);

    double * acf     = malloc(acmax * sizeof(double));
    double * troughs = malloc(acmax * sizeof(double));

    int nTroughs = 0;
    int jT   = -1;   // last trough before the peak under test; monotone
    int out  = 0;
    int have = 0;    // acf[0 .. have-1] computed, i.e. lags 1 .. have
    int nextI = 1;   // next acf index to classify
    int done = 0;

    // set to acmax to disable the early exit and always compute the full ACF
    const int CHUNK = 256;

    while(!done && have < acmax){

        int upto = have + CHUNK;
        if(upto > acmax) upto = acmax;

        PD_acf_block(ySub, size, have+1, upto, acf);
        have = upto;

        // classify everything that now has both neighbours available
        for(int i = nextI; i <= have-2; i++){

            const double slopeIn  = acf[i] - acf[i-1];
            const double slopeOut = acf[i+1] - acf[i];

            if(slopeIn < 0 && slopeOut > 0)
            {
                troughs[nTroughs] = i;
                nTroughs += 1;
            }
            else if(slopeIn > 0 && slopeOut < 0)
            {
                // a peak: test it now, in ascending order
                const int iPeak = i;
                const double thePeak = acf[iPeak];

                // find trough before this peak
                while(jT+1 < nTroughs && troughs[jT+1] < iPeak) jT++;
                if(jT == -1)
                    continue;

                const int iTrough = troughs[jT];
                const double theTrough = acf[iTrough];

                // (b) difference between peak and trough is at least 0.01
                if(thePeak - theTrough < th)
                    continue;

                // (c) peak corresponds to positive correlation
                if(thePeak < 0)
                    continue;

                // use this frequency that first fulfils all conditions.
                out = iPeak;
                done = 1;
                break;
            }
        }

        if(have-1 > nextI) nextI = have-1;
    }

    free(ySub);
    free(acf);
    free(troughs);

    return out;

}