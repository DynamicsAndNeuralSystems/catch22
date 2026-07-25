//
//  IN_AutoMutualInfoStats.c
//  C_polished
//
//  Created by Carl Henning Lubba on 22/09/2018.
//  Copyright © 2018 Carl Henning Lubba. All rights reserved.
//

#include <math.h>

#include "IN_AutoMutualInfoStats.h"
#include "CO_AutoCorr.h"
#include "stats.h"

double IN_AutoMutualInfoStats_40_gaussian_fmmi(const double y[], const int size)
{
    // NaN check
    for (int i = 0; i < size; i++) {
        if (isnan(y[i])) {
            return NAN;
        }
    }

    // Maximum time delay
    int tau = 40;

    // Don't go above half the signal length
    int max_tau = (size + 1) / 2;
    if (tau > max_tau) {
        tau = max_tau;
    }

    // If there are fewer than 3 lags, no local minimum can exist
    if (tau < 3) {
        return tau;
    }

    // Compute first two AMI values
    double ac = autocorr_lag(y, size, 1);
    double ami_prev = -0.5 * log(1.0 - ac * ac);

    ac = autocorr_lag(y, size, 2);
    double ami_curr = -0.5 * log(1.0 - ac * ac);

    double fmmi = tau;

    // Stream through remaining AMI values
    for (int i = 1; i < tau - 1; i++) {
        ac = autocorr_lag(y, size, i + 2);
        double ami_next = -0.5 * log(1.0 - ac * ac);

        if (ami_curr < ami_prev && ami_curr < ami_next) {
            fmmi = i;
            break;
        }

        ami_prev = ami_curr;
        ami_curr = ami_next;
    }

    return fmmi;
}
