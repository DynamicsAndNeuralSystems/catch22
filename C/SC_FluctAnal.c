#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <stdio.h>
#include <time.h>
#include "stats.h"
#include "CO_AutoCorr.h"

double SC_FluctAnal_2_50_1_logi_prop_r1(const double y[], const int size, const int lag, const char how[])
{
    // NaN check
    for (int i = 0; i < size; i++)
    {
        if (isnan(y[i]))
        {
            return NAN;
        }
    }

    // generate log spaced tau vector
    double linLow  = log(5);
    double linHigh = log(size/2);

    int nTauSteps = 50;
    double tauStep = (linHigh - linLow) / (nTauSteps-1);

    int tau[50];
    for (int i = 0; i < nTauSteps; i++)
    {
        tau[i] = round(exp(linLow + i*tauStep));
    }

    // check for uniqueness, use ascending order
    int nTau = nTauSteps;
    for (int i = 0; i < nTauSteps-1; i++)
    {
        while (tau[i] == tau[i+1] && i < nTau-1)
        {
            for (int j = i+1; j < nTauSteps-1; j++)
            {
                tau[j] = tau[j+1];
            }
            // lost one
            nTau -= 1;
        }
    }

    // fewer than 12 points -> leave.
    if (nTau < 12) {
        return 0;
    }

    // method dispatch resolved once instead of twice per window
    const int isDFA = (strcmp(how, "dfa") == 0);
    if (!isDFA && strcmp(how, "rsrangefit") != 0) {
        return 0.0;   // original returned 0.0 from inside the window loop
    }

    const int sizeCS = size/lag;
    double *yCS = malloc(sizeCS * sizeof *yCS);
    double *F   = malloc(nTau * sizeof *F);

    // transform input vector to cumsum
    yCS[0] = y[0];
    for (int i = 0; i < sizeCS-1; i++)
    {
        yCS[i+1] = yCS[i] + y[(i+1)*lag];
    }

    for (int i = 0; i < nTau; i++)
    {
        const int t       = tau[i];
        const int nBuffer = sizeCS/t;

        // x-moments of the regression depend only on tau -> hoisted out of the
        // window loop, same sequential summation order as linreg's inner loop
        double sumx = 0.0, sumx2 = 0.0;
        for (int k = 0; k < t; k++)
        {
            const double xv = (double)(k+1);
            sumx  += xv;
            sumx2 += xv*xv;
        }
        const double denom    = ((double)t * sumx2 - sumx*sumx);
        const int    singular = (denom == 0);

        double Fi = 0.0;

        for (int j = 0; j < nBuffer; j++)
        {
            const double *w = yCS + (size_t)j*t;

            double sumxy = 0.0, sumy = 0.0;
            for (int k = 0; k < t; k++)
            {
                const double xv = (double)(k+1);
                sumxy += xv * w[k];
                sumy  += w[k];
            }

            double m = 0.0, b = 0.0;
            if (!singular)
            {
                m = ((double)t * sumxy - sumx * sumy) / denom;
                b = (sumy * sumx2 - sumx * sumxy) / denom;
            }

            // detrend fused with the reduction: no scratch buffer
            if (isDFA)
            {
                for (int k = 0; k < t; k++)
                {
                    const double r = w[k] - (m * (k+1) + b);
                    Fi += r*r;
                }
            }
            else
            {
                double r  = w[0] - (m * (0+1) + b);
                double mx = r, mn = r;
                for (int k = 1; k < t; k++)
                {
                    r = w[k] - (m * (k+1) + b);
                    if (r > mx) mx = r;
                    if (r < mn) mn = r;
                }
                const double d = mx - mn;
                Fi += d*d;
            }
        }

        F[i] = isDFA ? sqrt(Fi/(nBuffer*t)) : sqrt(Fi/nBuffer);
    }

    double *logtt = malloc(nTau * sizeof *logtt);
    double *logFF = malloc(nTau * sizeof *logFF);
    int ntt = nTau;

    for (int i = 0; i < nTau; i++)
    {
        logtt[i] = log(tau[i]);
        logFF[i] = log(F[i]);
    }

    int minPoints = 6;
    int nsserr = (ntt - 2*minPoints + 1);
    double *sserr  = malloc(nsserr * sizeof *sserr);
    double *buffer = malloc((ntt - minPoints + 1) * sizeof *buffer);
    for (int i = minPoints; i < ntt - minPoints + 1; i++)
    {
        double m1 = 0.0, b1 = 0.0;
        double m2 = 0.0, b2 = 0.0;

        sserr[i - minPoints] = 0.0;

        linreg(i, logtt, logFF, &m1, &b1);
        linreg(ntt-i+1, logtt+i-1, logFF+i-1, &m2, &b2);

        for (int j = 0; j < i; j++)
        {
            buffer[j] = logtt[j] * m1 + b1 - logFF[j];
        }

        sserr[i - minPoints] += norm_(buffer, i);

        for (int j = 0; j < ntt-i+1; j++)
        {
            buffer[j] = logtt[j+i-1] * m2 + b2 - logFF[j+i-1];
        }

        sserr[i - minPoints] += norm_(buffer, ntt-i+1);
    }

    double firstMinInd = 0.0;
    double minimum = min_(sserr, nsserr);
    for (int i = 0; i < nsserr; i++)
    {
        if (sserr[i] == minimum)
        {
            firstMinInd = i + minPoints - 1;
            break;
        }
    }

    free(yCS);
    free(F);
    free(logtt);
    free(logFF);
    free(sserr);
    free(buffer);

    return (firstMinInd+1)/ntt;
}

double SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1(const double y[], const int size){
    return SC_FluctAnal_2_50_1_logi_prop_r1(y, size, 2, "dfa");
}

double SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1(const double y[], const int size){
    return SC_FluctAnal_2_50_1_logi_prop_r1(y, size, 1, "rsrangefit");
}