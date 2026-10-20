#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "SB_CoarseGrain.h"
#include "helper_functions.h"

double SB_MotifThree_quantile_hh(const double y[], const int size)
{
    // NaN check
    for(int i = 0; i < size; i++)
    {
        if(isnan(y[i]))
        {
            return NAN;
        }
    }
    
    int tmp_idx, r_idx;
    int dynamic_idx;
    int alphabet_size = 3;
    int array_size;
    int * yt = malloc(size * sizeof(yt)); // alphabetized array
    double hh; // output
    double * out = malloc(124 * sizeof(out)); // output array
    
    // transfer to alphabet
    sb_coarsegrain(y, size, "quantile", 3, yt);
    
    // words of length 1
    array_size = alphabet_size;
    int ** r1 = malloc(array_size * sizeof(*r1));
    int * sizes_r1 = malloc(array_size * sizeof(sizes_r1));
    double * out1 = malloc(array_size * sizeof(out1));
    for (int i = 0; i < alphabet_size; i++) {
        r1[i] = malloc(size * sizeof(r1[i])); // probably can be rewritten
        // using selfresizing array for memory efficiency. Time complexity
        // should be comparable due to ammotization.
        r_idx = 0;
        sizes_r1[i] = 0;
        for (int j = 0; j < size; j++) {
            if (yt[j] == i + 1) {
                r1[i][r_idx++] = j;
                sizes_r1[i]++;
            }
        }
    }
    
    // words of length 2
    array_size *= alphabet_size;
    // removing last item if it is == max possible idx since later we are taking idx + 1
    // from yt
    for (int i = 0; i < alphabet_size; i++) {
        if (sizes_r1[i] != 0 && r1[i][sizes_r1[i] - 1] == size - 1) {
            //int * tmp_ar = malloc((sizes_r1[i] - 1) * sizeof(tmp_ar));
            int* tmp_ar = malloc(sizes_r1[i] * sizeof(tmp_ar));
            subset(r1[i], tmp_ar, 0, sizes_r1[i]);
            memcpy(r1[i], tmp_ar, (sizes_r1[i] - 1) * sizeof(tmp_ar));
            sizes_r1[i]--;
            free(tmp_ar);
        }
    }
    
    /*
    int *** r2 = malloc(array_size * sizeof(**r2));
    int ** sizes_r2 = malloc(array_size * sizeof(*sizes_r2));
    double ** out2 = malloc(array_size * sizeof(*out2));
    */
    int*** r2 = malloc(alphabet_size * sizeof(**r2));
    int** sizes_r2 = malloc(alphabet_size * sizeof(*sizes_r2));
    double** out2 = malloc(alphabet_size * sizeof(*out2));
    

    // allocate separately
    for (int i = 0; i < alphabet_size; i++) {
        r2[i] = malloc(alphabet_size * sizeof(*r2[i]));
        sizes_r2[i] = malloc(alphabet_size * sizeof(*sizes_r2[i]));
        //out2[i] = malloc(alphabet_size * sizeof(out2[i]));
        out2[i] = malloc(alphabet_size * sizeof(**out2));
        for (int j = 0; j < alphabet_size; j++) {
            r2[i][j] = malloc(size * sizeof(*r2[i][j]));
        }
    }

    // fill separately
    for (int i = 0; i < alphabet_size; i++) {
    // for (int i = 0; i < array_size; i++) {
        //r2[i] = malloc(alphabet_size * sizeof(r2[i]));
        //sizes_r2[i] = malloc(alphabet_size * sizeof(sizes_r2[i]));
        //out2[i] = malloc(alphabet_size * sizeof(out2[i]));
        for (int j = 0; j < alphabet_size; j++) {
            //r2[i][j] = malloc(size * sizeof(r2[i][j]));
            sizes_r2[i][j] = 0;
            dynamic_idx = 0; //workaround as you can't just add elements to array
            // like in python (list.append()) for example, so since for some k there will be no adding,
            // you need to keep track of the idx at which elements will be inserted
            for (int k = 0; k < sizes_r1[i]; k++) {
                tmp_idx = yt[r1[i][k] + 1];
                if (tmp_idx == (j + 1)) {
                    r2[i][j][dynamic_idx++] = r1[i][k];
                    sizes_r2[i][j]++;
                    // printf("dynamic_idx=%i, size = %i\n", dynamic_idx, size);
                }
            }
            double tmp = (double)sizes_r2[i][j] / ((double)(size) - (double)(1.0));
            out2[i][j] =  tmp;
        }
    }

    hh = 0.0;
    int nonEmpty = 0;
    for (int i = 0; i < alphabet_size; i++) {
        hh += f_entropy(out2[i], alphabet_size);
        for (int j = 0; j < alphabet_size; j++) {
            if (out2[i][j] > 0) {
                nonEmpty++;
            }
        }
    }
    // Miller-Madow bias correction (as hctsa's SB_MotifThree), with degrees of
    // freedom constrained by the (fixed, equiprobable) symbol marginals:
    // df = (non-empty cells) - 1 - wordLength*(alphabetSize - 1)
    const int df = nonEmpty - 1 - 2 * (alphabet_size - 1);
    if (df > 0) {
        hh += df / (2.0 * (size - 1));
    }

    free(yt);
    free(out);
    free(out1);

    free(sizes_r1);

    // free nested array
    for (int i = 0; i < alphabet_size; i++) {
        free(r1[i]);
    }
    free(r1);
    // free(sizes_r1);
    
    for (int i = 0; i < alphabet_size; i++) {
    //for (int i = alphabet_size - 1; i >= 0; i--) {

        free(sizes_r2[i]);
        free(out2[i]);
    }

    //for (int i = alphabet_size-1; i >= 0 ; i--) {
    for(int i = 0; i < alphabet_size; i++) {
        for (int j = 0; j < alphabet_size; j++) {
            free(r2[i][j]);
        }
        free(r2[i]);
    }
    
    free(r2);
    free(sizes_r2);
    free(out2);
    
    
    return hh;
    
}

























