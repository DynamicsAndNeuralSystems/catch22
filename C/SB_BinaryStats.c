//
//  SB_BinaryStats.c
//  C_polished
//
//  Created by Carl Henning Lubba on 22/09/2018.
//  Copyright © 2018 Carl Henning Lubba. All rights reserved.
//

#include "SB_BinaryStats.h"
#include "stats.h"

// Length of the longest run of consecutive entries of yBin equal to val
static int longest_run(const int yBin[], const int size, const int val)
{
    int maxRun = 0;
    int run = 0;
    for (int i = 0; i < size; i++) {
        if (yBin[i] == val) {
            run++;
            if (run > maxRun) {
                maxRun = run;
            }
        }
        else {
            run = 0;
        }
    }
    return maxRun;
}

double SB_BinaryStats_diff_longstretch0(const double y[], const int size){
    
    // NaN check
    for(int i = 0; i < size; i++)
    {
        if(isnan(y[i]))
        {
            return NAN;
        }
    }
    
    // binarize: 1 for a stepwise increase, 0 otherwise (as in hctsa's BF_Binarize)
    int * yBin = malloc((size-1) * sizeof(int));
    for(int i = 0; i < size-1; i++){
        yBin[i] = (y[i+1] - y[i] > 0) ? 1 : 0;
    }
    
    // longest stretch of non-increasing steps
    int maxstretch0 = longest_run(yBin, size-1, 0);
    
    free(yBin);
    
    return maxstretch0;
}

double SB_BinaryStats_mean_longstretch1(const double y[], const int size){
    
    // NaN check
    for(int i = 0; i < size; i++)
    {
        if(isnan(y[i]))
        {
            return NAN;
        }
    }
    
    // binarize: 1 for values above the mean, 0 otherwise
    int * yBin = malloc(size * sizeof(int));
    double yMean = mean(y, size);
    for(int i = 0; i < size; i++){
        yBin[i] = (y[i] - yMean > 0) ? 1 : 0;
    }
    
    // longest stretch of above-mean values
    int maxstretch1 = longest_run(yBin, size, 1);
    
    free(yBin);
    
    return maxstretch1;
}
