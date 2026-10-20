//
//  histcounts.h
//  C_polished
//
//  Created by Carl Henning Lubba on 19/09/2018.
//  Copyright © 2018 Carl Henning Lubba. All rights reserved.
//

#ifndef histcounts_h
#define histcounts_h

#include <stdlib.h>
#include <float.h>
#include <stdio.h>

extern double iqr_matlab(const double y[], const int size);
extern int hist_num_bins_auto(const double y[], const int size);
extern int hist_edges(const double y[], const int size, const int numBins, double edges[]);
extern int hist_bin_index(const double x, const double edges[], const int numBins);
extern void hist_counts_edges(const double y[], const int size, const double edges[], const int numBins, int counts[]);

#endif /* histcounts_h */
