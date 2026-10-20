#ifndef STATS_H
#define STATS_H
#include <math.h>
#include <stdlib.h>
#include <string.h>

extern double max_(const double a[], const int size);
extern double min_(const double a[], const int size);
extern double mean(const double a[], const int size);
extern double sum(const double a[], const int size);
extern void icumsum(const int a[], const int size, int b[]);
extern double median(const double a[], const int size);
extern double stddev(const double a[], const int size);
extern double corr(const double x[], const double y[], const int size);
extern double cov(const double x[], const double y[], const int size);
extern double autocorr_lag(const double x[], const int size, const int lag);
extern void zscore_norm2(const double a[], const int size, double b[]);
extern void diff(const double a[], const int size, double b[]);
extern int linreg(const int n, const double x[], const double y[], double* m, double* b);

#endif
