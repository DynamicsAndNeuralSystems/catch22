//
//  MD_hrv.c
//  C_polished
//
//  Created by Carl Henning Lubba on 22/09/2018.
//  Copyright © 2018 Carl Henning Lubba. All rights reserved.
//

#include "MD_hrv.h"
#include "stats.h"

double MD_pNN_pnn40(const double y[], const int size){

    const int pNNx = 40;
    
    // need at least one difference
    if (size < 2) {
        return NAN;
    }
    
    // compute diff
    double * Dy = malloc((size-1) * sizeof(double));
    diff(y, size, Dy);
    
    double pnn40 = 0;
    for(int i = 0; i < size-1; i++){
        if(fabs(Dy[i])*1000 > pNNx){
            pnn40 += 1;
        }
    }
    
    free(Dy);
    
    return pnn40/(size-1);
}

