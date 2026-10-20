//  Created by Carl Henning Lubba on 27/09/2018.
//  Copyright © 2018 Carl Henning Lubba. All rights reserved.
//
// Based on the work of Jonas Lundgren in his Matlab Central contribution 'SPLINEFIT'.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "splinefit.h"
#include "stats.h"

#define nCoeffs 3
#define nPoints 4

#define pieces 2
#define nBreaks 3
#define deg 3
#define nSpline 4
#define piecesExt 8 //3 * deg - 1


void matrix_multiply(const int sizeA1, const int sizeA2, const double *A, const int sizeB1, const int sizeB2, const double *B, double *C){
//void matrix_multiply(int sizeA1, int sizeA2, double **A, int sizeB1, int sizeB2, double **B, double C[sizeA1][sizeB2]){
    
    if(sizeA2 != sizeB1){
        return;
    }
    
    for(int i = 0; i < sizeA1; i++){
        for(int j = 0; j < sizeB2; j++){
            
            C[i*sizeB2 + j] = 0;
            for(int k = 0; k < sizeB1; k++){
                C[i*sizeB2 + j] += A[i * sizeA2 + k]*B[k * sizeB2 + j];
            }
            
        }
    }
    
}

void matrix_times_vector(const int sizeA1, const int sizeA2, const double *A, const int sizeb, const double *b, double *c){ //c[sizeb]
    
    if(sizeA2 != sizeb){
        return;
    }
    
    // row
    for(int i = 0; i < sizeA1; i++){
        
        // column
        c[i] = 0;
        for(int k = 0; k < sizeb; k++){
            c[i] += A[i * sizeA2 + k]*b[k];
        }
        
    }
    
}

void gauss_elimination(int size, double *A, double *b, double *x){
// void gauss_elimination(int size, double A[size][size], double b[size], double x[size]){
    
    double factor;
    
    // create temp matrix and vector
    double* AElim[nSpline + 1];
    for (int i = 0; i < size; i++)
        AElim[i] = (double *)malloc(size * sizeof(double));
    double * bElim = malloc(size * sizeof(double));
    
    // -- create triangular matrix
    
    // initialise to A and b
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            AElim[i][j] = A[i*size + j];
        }
        bElim[i] = b[i];
    }
    
    // go through columns in outer loop
    for(int i = 0; i < size; i++){
    
        // go through rows to eliminate
        for(int j = i+1; j < size; j++){
            
            factor = AElim[j][i]/AElim[i][i];
            
            // subtract in vector
            bElim[j] = bElim[j] - factor*bElim[i];
            
            // go through entries of this row
            for(int k = i; k < size; k++){
                AElim[j][k] = AElim[j][k] - factor*AElim[i][k];
            }
            
        }
        
    }
    
    
    // -- go backwards through triangular matrix and solve for x
    
    // row
    double bMinusATemp;
    for(int i = size-1; i >= 0; i--){
        
        bMinusATemp = bElim[i];
        for(int j = i+1; j < size; j++){
            bMinusATemp -= x[j]*AElim[i][j];
        }
        
        x[i] = bMinusATemp/AElim[i][i];
    }
    
    for (int i = 0; i < size; i++)
        free(AElim[i]);
    free(bElim);
}

void lsqsolve_sub(const int sizeA1, const int sizeA2, const double *A, const int sizeb, const double *b, double *x)
//void lsqsolve_sub(int sizeA1, int sizeA2, double A[sizeA1][sizeA2], int sizeb, double b[sizeb], double x[sizeA1])
{
    (void) sizeb; // unused
    // create temp matrix and vector
    
    double * AT = calloc(sizeA2 * sizeA1, sizeof(double)); // zeroed: silences a false maybe-uninitialized warning
    double * ATA = malloc(sizeA2 * sizeA2 * sizeof(double));
    double * ATb = malloc(sizeA2 * sizeof(double));
    
    
    for(int i = 0; i < sizeA1; i++){
        for(int j = 0; j < sizeA2; j++){
            //AT[i,j] = A[j,i]
            AT[j * sizeA1 + i] = A[i * sizeA2 + j];
        }
    }
    
    
    matrix_multiply(sizeA2, sizeA1, AT, sizeA1, sizeA2, A, ATA);
    
    
    
    matrix_times_vector(sizeA2, sizeA1, AT, sizeA1, b, ATb);
    
    
    gauss_elimination(sizeA2, ATA, ATb, x);
    
    free(AT);
    free(ATA);
    free(ATb);
    
}

int iLimit(int x, int lim){
    return x<lim ? x : lim;
}

int splinefit(const double *y, const int size, double *yOut)
{
    // degree of spline
    //const int nSpline = 4;
    //const int deg = 3;
    
    // x-positions of spline-nodes
    //const int nBreaks = 3;
    int breaks[nBreaks];
    breaks[0] = 0;
    breaks[1] = (int)floor((double)size/2.0)-1;
    breaks[2] = size-1;
    
    // -- splinebase
    
    // spacing
    int h0[2];
    h0[0] = breaks[1] - breaks[0];
    h0[1] = breaks[2] - breaks[1];
    
    //const int pieces = 2;
    
    // repeat spacing
    int hCopy[4];
    hCopy[0] = h0[0], hCopy[1] = h0[1], hCopy[2] = h0[0], hCopy[3] = h0[1];
    
    // to the left
    int hl[deg];
    hl[0] = hCopy[deg-0];
    hl[1] = hCopy[deg-1];
    hl[2] = hCopy[deg-2];
    
    int hlCS[deg]; // cumulative sum
    icumsum(hl, deg, hlCS);
    
    int bl[deg];
    for(int i = 0; i < deg; i++){
        bl[i] = breaks[0] - hlCS[i];
    }
    
    // to the left
    int hr[deg];
    hr[0] = hCopy[0];
    hr[1] = hCopy[1];
    hr[2] = hCopy[2];
    
    int hrCS[deg]; // cumulative sum
    icumsum(hr, deg, hrCS);
    
    int br[deg];
    for(int i = 0; i < deg; i++){
        br[i] = breaks[2] + hrCS[i];
    }
    
    // add breaks
    int breaksExt[3*deg];
    for(int i = 0; i < deg; i++){
        breaksExt[i] = bl[deg-1-i];
        breaksExt[i + 3] = breaks[i];
        breaksExt[i + 6] = br[i];
    }
    int hExt[3*deg-1];
    for(int i = 0; i < deg*3-1; i++){
        hExt[i] = breaksExt[i+1] - breaksExt[i];
    }
    //const int piecesExt = 3*deg-1;
    
    // initialise polynomial coefficients
    double coefs[nSpline*piecesExt][nSpline+1];
    for(int i = 0; i < nSpline*piecesExt; i++){
        for(int j = 0; j < nSpline; j++){
        coefs[i][j] = 0;
        }
    }
    for(int i = 0; i < nSpline*piecesExt; i=i+nSpline){
        coefs[i][0] = 1;
    }
    
    // expand h using the index matrix ii
    int ii[deg+1][piecesExt];
    for(int i = 0; i < piecesExt; i++){
        ii[0][i] = iLimit(0+i, piecesExt-1);
        ii[1][i] = iLimit(1+i, piecesExt-1);
        ii[2][i] = iLimit(2+i, piecesExt-1);
        ii[3][i] = iLimit(3+i, piecesExt-1);
    }
    
    // expanded h
    double H[(deg+1)*piecesExt];
    int iiFlat;
    for(int i = 0; i < nSpline*piecesExt; i++){
        iiFlat = ii[i%nSpline][i/nSpline];
        H[i] = hExt[iiFlat];
    }
    
    //recursive generation of B-splines
    double Q[nSpline][piecesExt];
    for(int k = 1; k < nSpline; k++){
        
        //antiderivatives of splines
        for(int j = 0; j<k; j++){
            for(int l = 0; l<(nSpline*piecesExt); l++){
                coefs[l][j] *= H[l]/(k-j);
            }
        }
        
        for(int l = 0; l<(nSpline*piecesExt); l++){
            Q[l%nSpline][l/nSpline] = 0;
            for(int m = 0; m < nSpline; m++){
                Q[l%nSpline][l/nSpline] += coefs[l][m];
            }
        }
        
        //cumsum
        for(int l = 0; l<piecesExt; l++){
            for(int m = 1; m < nSpline; m++){
                Q[m][l] += Q[m-1][l];
            }
        }
        
        for(int l = 0; l<nSpline*piecesExt; l++){
            if(l%nSpline == 0)
                coefs[l][k] = 0;
            else{
                coefs[l][k] = Q[l%nSpline-1][l/nSpline]; // questionable
            }
        }
        
        // normalise antiderivatives by max value
        double fmax[piecesExt*nSpline];
        for(int i = 0; i < piecesExt; i++){
            for(int j = 0; j < nSpline; j++){
                
                fmax[i*nSpline+j] = Q[nSpline-1][i];
                
            }
        }
        
        for(int j = 0; j < k+1; j++){
            for(int l = 0; l < nSpline*piecesExt; l++){
                coefs[l][j] /= fmax[l];
            }
        }

        // diff to adjacent antiderivatives
        for(int i = 0; i < (nSpline*piecesExt)-deg; i++){
            for(int j = 0; j < k+1; j ++){
                coefs[i][j] -= coefs[deg+i][j];
            }
        }
        for(int i = 0; i < nSpline*piecesExt; i += nSpline){
            coefs[i][k] = 0;
        }
        
    }
    
    // scale coefficients
    double scale[nSpline*piecesExt];
    for(int i = 0; i < nSpline*piecesExt; i++)
    {
        scale[i] = 1;
    }
    for(int k = 0; k < nSpline-1; k++){
        for(int i = 0; i < (nSpline*piecesExt); i++){
            scale[i] /= H[i];
        }
        for(int i = 0; i < (nSpline*piecesExt); i++){
            coefs[i][(nSpline-1)-(k+1)] *= scale[i];
        }
    }
    
    // reduce pieces and sort coefficients by interval number
    int jj[nSpline][pieces];
    for(int i = 0; i < nSpline; i++){
        for(int j = 0; j < pieces; j++){
            if(i == 0)
                jj[i][j] = nSpline*(1+j);
            else
                jj[i][j] = deg;
        }
    }
    
    
    for(int i = 1; i < nSpline; i++){
        for(int j = 0; j < pieces; j++){
            jj[i][j] += jj[i-1][j];
        }
    }
    
    double coefsOut[nSpline*pieces][nSpline];
    int jj_flat;
    for(int i = 0; i < nSpline*pieces; i++){
        jj_flat = jj[i%nSpline][i/nSpline]-1;
        for(int j = 0; j < nSpline; j++){
            coefsOut[i][j] = coefs[jj_flat][j];
        }
        
    }
    
    
    // -- create first B-splines to feed into optimization
    
    // x-values for B-splines
    int * xsB = malloc((size*nSpline)* sizeof(int));
    int * indexB = malloc((size*nSpline) * sizeof(int));
    
    int breakInd = 1;
    for(int i = 0; i < size; i++){
        if(i >= breaks[breakInd] && breakInd<nBreaks-1)
            breakInd += 1;
        for(int j = 0; j < nSpline; j++){
            xsB[i*nSpline+j] = i - breaks[breakInd-1];
            indexB[i*nSpline+j] = j + (breakInd-1)*nSpline;
        }
    }
    
    double * vB = malloc((size*nSpline) * sizeof(double));
    for(int i = 0; i < size*nSpline; i++){
        vB[i] = coefsOut[indexB[i]][0];
    }
    
    for(int i = 1; i < nSpline; i ++){
        for(int j = 0; j < size*nSpline; j++){
            vB[j] = vB[j]*xsB[j] + coefsOut[indexB[j]][i];
        }
    }
    
    
    double * A = malloc(size*(nSpline+1) * sizeof(double));
    
    for(int i = 0; i < (nSpline+1)*size; i++){
        A[i] = 0;
    }
    breakInd = 0;
    for(int i = 0; i < nSpline*size; i++){
        if(i/nSpline >= breaks[1])
            breakInd = 1;
        A[(i%nSpline)+breakInd + (i/nSpline)*(nSpline+1)] = vB[i];
    }
    
    
    
    double * x = malloc((nSpline+1)*sizeof(double));
    // lsqsolve_sub(int sizeA1, int sizeA2, double *A, int sizeb, double *b, double *x)
    lsqsolve_sub(size, nSpline+1, A, size, y, x);
    
    // coeffs of B-splines to combine by optimised weighting in x
    double C[pieces+nSpline-1][nSpline*pieces];
    // initialise to 0
    for(int i = 0; i < nSpline+1; i++){
        for(int j = 0; j < nSpline*pieces; j++){
            C[i][j] = 0;
        }
    }
    
    int CRow, CCol, coefRow, coefCol;
    for(int i = 0; i < nSpline*nSpline*pieces; i++){
        
        CRow = i%nSpline + (i/nSpline)%2;
        CCol = i/nSpline;
        
        coefRow = i%(nSpline*2);
        coefCol =i/(nSpline*2);
        
        C[CRow][CCol] = coefsOut[coefRow][coefCol];
        
    }
    
    // final coefficients
    double coefsSpline[pieces][nSpline];
    for(int i = 0; i < pieces; i++){
        for(int j = 0; j < nSpline; j++){
            coefsSpline[i][j] = 0;
        }
    }
    
    //multiply with x
    for(int j = 0; j < nSpline*pieces; j++){
        coefCol = j/pieces;
        coefRow = j%pieces;
        
        for(int i = 0; i < nSpline+1; i++){
            
            coefsSpline[coefRow][coefCol] += C[i][j]*x[i];
            
        }
    }
    
    
    // compute piecewise polynomial
    
    int secondHalf = 0;
    for(int i = 0; i < size; i++){
        secondHalf = i < breaks[1] ? 0 : 1;
        yOut[i] = coefsSpline[secondHalf][0];
    }
    
    for(int i = 1; i < nSpline; i ++){
        for(int j = 0; j < size; j++){
            secondHalf = j < breaks[1] ? 0 : 1;
            yOut[j] = yOut[j]*(j - breaks[1]*secondHalf) + coefsSpline[secondHalf][i];
        }
        
    }
    
    free(xsB);
    free(indexB);
    free(vB);
    free(A);
    free(x);
    
    return 0;
    
}


