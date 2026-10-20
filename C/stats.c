#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "helper_functions.h"

double min_(const double a[], const int size)
{
    double m = a[0];
    for (int i = 1; i < size; i++) {
        if (a[i] < m) {
            m = a[i];
        }
    }
    return m;
}

double max_(const double a[], const int size)
{
    double m = a[0];
    for (int i = 1; i < size; i++) {
        if (a[i] > m) {
            m = a[i];
        }
    }
    return m;
}

double mean(const double a[], const int size)
{
    double m = 0.0;
    for (int i = 0; i < size; i++) {
        m += a[i];
    }
    m /= size;
    return m;
}

double sum(const double a[], const int size)
{
    double m = 0.0;
    for (int i = 0; i < size; i++) {
        m += a[i];
    }
    return m;
}

void icumsum(const int a[], const int size, int b[])
{
    b[0] = a[0];
    for (int i = 1; i < size; i++) {
        b[i] = a[i] + b[i-1];
    }
    
}

double median(const double a[], const int size)
{
    double m;
    double * b = malloc(size * sizeof *b);
    memcpy(b, a, size * sizeof *b);
    sort(b, size);
    if (size % 2 == 1) {
        m = b[size / 2];
    } else {
        int m1 = size / 2;
        int m2 = m1 - 1;
        m = (b[m1] + b[m2]) / (double)2.0;
    }
    free(b);
    return m;
}

// Sample standard deviation. The deviations are scaled by a power of two before
// squaring, which is exact, so the squares neither overflow nor underflow (below
// about 1.5e-154 they would lose precision as subnormals) and the result is bit for
// bit the same wherever the unscaled squares were normal doubles. The scale comes
// from max |a|, found in the mean's pass (a deviation is at most twice it; a series
// whose values are 1e154 times their spread is constant to double precision anyway).
double stddev(const double a[], const int size)
{
    double sum = 0.0, maxAbs = 0.0;
    for (int i = 0; i < size; i++) {
        sum += a[i];
        const double d = fabs(a[i]);
        if (d > maxAbs) maxAbs = d;
    }
    const double m = sum / size; // as mean()
    int e = 0;
    frexp(maxAbs, &e); // maxAbs = f * 2^e with 0.5 <= f < 1 (e = 0 for an all-zero series)
    const double scale = ldexp(1.0, -e);
    double sd = 0.0;
    for (int i = 0; i < size; i++) {
        const double d = (a[i] - m) * scale;
        sd += d * d;
    }
    return ldexp(sqrt(sd / (size - 1)), e);
}

double cov(const double x[], const double y[], const int size){
    
    double covariance = 0;
    
    double meanX = mean(x, size);
    double meanY = mean(y, size);
    
    for(int i = 0; i < size; i++){
        covariance += (x[i] - meanX) * (y[i] - meanY);
        
    }
    
    return covariance/(size-1);
    
}

double corr(const double x[], const double y[], const int size){
    
    double nom = 0;
    double denomX = 0;
    double denomY = 0;
    
    double meanX = mean(x, size);
    double meanY = mean(y, size);
    
    for(int i = 0; i < size; i++){
        nom += (x[i] - meanX) * (y[i] - meanY);
        denomX += (x[i] - meanX) * (x[i] - meanX);
        denomY += (y[i] - meanY) * (y[i] - meanY);
        
    }
    
    return nom/sqrt(denomX * denomY);
    
}

double autocorr_lag(const double x[], const int size, const int lag){
    
    return corr(x, &(x[lag]), size-lag);
    
}

void zscore_norm2(const double a[], const int size, double b[])
{
    double m = mean(a, size);
    double sd = stddev(a, size);
    for (int i = 0; i < size; i++) {
        b[i] = (a[i] - m) / sd;
    }
    return;
}

void diff(const double a[], const int size, double b[])
{
    for (int i = 1; i < size; i++) {
        b[i - 1] = a[i] - a[i - 1];
    }
}

int linreg(const int n, const double x[], const double y[], double* m, double* b)
{
    double   sumx = 0.0;                      /* sum of x     */
    double   sumx2 = 0.0;                     /* sum of x**2  */
    double   sumxy = 0.0;                     /* sum of x * y */
    double   sumy = 0.0;                      /* sum of y     */
    
    for (int i=0;i<n;i++){
        sumx  += x[i];
        sumx2 += x[i] * x[i];
        sumxy += x[i] * y[i];
        sumy  += y[i];
    }
    
    double denom = (n * sumx2 - sumx * sumx);
    if (denom == 0) {
        // singular matrix. can't solve the problem.
        *m = 0;
        *b = 0;
        return 1;
    }
    
    *m = (n * sumxy  -  sumx * sumy) / denom;
    *b = (sumy * sumx2  -  sumx * sumxy) / denom;
    
    return 0;
}
