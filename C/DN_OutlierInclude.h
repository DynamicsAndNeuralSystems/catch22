#ifndef DN_OUTLIERINCLUDE_H
#define DN_OUTLIERINCLUDE_H
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include "stats.h"

extern double DN_OutlierInclude_np_001_mdrmd(const double y[], const int size, const int sign);
extern double DN_OutlierInclude_p_001_mdrmd(const double y[], const int size);
extern double DN_OutlierInclude_n_001_mdrmd(const double y[], const int size);

#endif
