//
//  SP_Summaries.c
//  
//
//  Created by Carl Henning Lubba on 23/09/2018.
//

#include "SP_Summaries.h"
#include "CO_AutoCorr.h"

// Single-sided power spectrum of the (z-scored) series, matching hctsa's
// SP_Summaries(y,'fft'): zero-padded to the next power of two, with the DC bin
// dropped. On return, S[k] and w[k] (k = 0, ..., nBins-1) hold the power density
// in angular frequency and the angular frequency of FFT bins 1, ..., NFFT/2.
// Returns nBins, or 0 if the spectrum cannot be computed.
static int sp_fft_spectrum(const double y[], const int size, double ** S, double ** w)
{
    const double PI = 3.14159265358979323846;
    const int NFFT = nextpow2(size);
    const int nBins = NFFT / 2;
    if (nBins < 1) {
        return 0;
    }

    cplx * F = malloc(NFFT * sizeof *F);
    cplx * tw = malloc(NFFT * sizeof *tw);
    for (int i = 0; i < NFFT; i++) {
        #if defined(__GNUC__) || defined(__GNUG__)
            F[i] = (i < size ? y[i] : 0.0) + 0.0 * I;
        #elif defined(_MSC_VER)
            cplx tmp = { (i < size ? y[i] : 0.0), 0.0 };
            F[i] = tmp;
        #endif
    }
    twiddles(tw, NFFT);
    fft(F, NFFT, tw);

    *S = malloc(nBins * sizeof **S);
    *w = malloc(nBins * sizeof **w);
    for (int k = 0; k < nBins; k++) {
        const double a = cabs(F[k + 1]);
        (*S)[k] = 2 * a * a / size / (2 * PI); // single-sided density, in angular frequency
        (*w)[k] = 2 * PI * (k + 1) / NFFT;
    }
    free(F);
    free(tw);
    return nBins;
}

static double sp_fft_summary(const double y[], const int size, const char what[])
{
    for (int i = 0; i < size; i++) {
        if (!isfinite(y[i])) {
            return NAN;
        }
    }
    double * S;
    double * w;
    const int nBins = sp_fft_spectrum(y, size, &S, &w);
    if (nBins < 2) {
        if (nBins > 0) {
            free(S);
            free(w);
        }
        return NAN;
    }
    const double dw = w[1] - w[0];

    double output = NAN;
    if (strcmp(what, "wq_50") == 0) {
        // the frequency below which half of the power lies (the median frequency)
        double total = 0;
        for (int k = 0; k < nBins; k++) {
            total += S[k];
        }
        double cs = 0;
        for (int k = 0; k < nBins; k++) {
            cs += S[k];
            if (cs >= 0.5 * total) {
                output = w[k];
                break;
            }
        }
    }
    else if (strcmp(what, "area_5_1") == 0) {
        // power in the lowest fifth of the frequency bins
        double area = 0;
        for (int k = 0; k < nBins / 5; k++) {
            area += S[k];
        }
        output = area * dw;
    }
    free(S);
    free(w);
    return output;
}

double SP_Summaries_fft_area_5_1(const double y[], const int size)
{
    return sp_fft_summary(y, size, "area_5_1");
}

double SP_Summaries_fft_wq_50(const double y[], const int size)
{
    return sp_fft_summary(y, size, "wq_50");
}
