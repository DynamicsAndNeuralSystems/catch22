// catch22_mex   MATLAB gateway to catch22_compute() (compile with mexAll.m).
//
// [values, isMasked] = catch22_mex(x, doCatch24, maskNaNs)
//   x: a real double vector (the time series).
//   doCatch24: nonzero to add the mean and standard deviation (24 features).
//   maskNaNs: nonzero to set features that cannot be estimated to a fixed value
//       (see catch22.h); NaN/Inf in x is never masked (all outputs are NaN).
//   values: column vector of the 22 (or 24) feature values.
//   isMasked: logical column vector, true for each masked output.
//
// [names, shortNames, legacyNames, hctsaCalls] = catch22_mex('features', doCatch24)
//   Cell arrays (column) describing the features, in the order of values.

#include <limits.h>
#include <string.h>
#include "mex.h"
#include "catch22.h"

static int scalar_flag(const mxArray * a, const char * what)
{
    if (!(mxIsNumeric(a) || mxIsLogical(a)) || mxGetNumberOfElements(a) != 1) {
        mexErrMsgIdAndTxt("catch22:badInput", "%s must be a scalar (true or false).", what);
    }
    return mxGetScalar(a) != 0;
}

static mxArray * name_cell(const char * (*name)(const int), const int nOut)
{
    mxArray * c = mxCreateCellMatrix(nOut, 1);
    for (int i = 0; i < nOut; i++) {
        mxSetCell(c, i, mxCreateString(name(i)));
    }
    return c;
}

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    if (nrhs >= 1 && mxIsChar(prhs[0])) {
        char cmd[16];
        if (mxGetString(prhs[0], cmd, sizeof cmd) != 0 || strcmp(cmd, "features") != 0) {
            mexErrMsgIdAndTxt("catch22:badInput", "The only command is 'features'.");
        }
        if (nlhs > 4) {
            mexErrMsgIdAndTxt("catch22:badOutput", "'features' returns at most 4 outputs.");
        }
        const int catch24 = (nrhs >= 2) ? scalar_flag(prhs[1], "doCatch24") : 0;
        const int nOut = catch22_n_features(catch24);
        const char * (*names[4])(const int) = {catch22_feature_name, catch22_short_name,
                                               catch22_legacy_name, catch22_hctsa_call};
        for (int k = 0; k < 4 && k < (nlhs > 0 ? nlhs : 1); k++) {
            plhs[k] = name_cell(names[k], nOut);
        }
        return;
    }

    if (nrhs < 1 || nrhs > 3) {
        mexErrMsgIdAndTxt("catch22:badInput", "Usage: [values, isMasked] = catch22_mex(x, doCatch24, maskNaNs)");
    }
    if (nlhs > 2) {
        mexErrMsgIdAndTxt("catch22:badOutput", "catch22_mex returns at most 2 outputs.");
    }
    const mxArray * x = prhs[0];
    if (!mxIsDouble(x) || mxIsComplex(x) || mxIsSparse(x)
        || (mxGetM(x) != 1 && mxGetN(x) != 1 && !mxIsEmpty(x))) {
        mexErrMsgIdAndTxt("catch22:badInput", "x must be a real, full, double vector.");
    }
    if (mxGetNumberOfElements(x) > INT_MAX) {
        mexErrMsgIdAndTxt("catch22:badInput", "x is too long.");
    }
    const int n = (int)mxGetNumberOfElements(x);
    const int catch24 = (nrhs >= 2) ? scalar_flag(prhs[1], "doCatch24") : 0;
    const int maskNaNs = (nrhs >= 3) ? scalar_flag(prhs[2], "maskNaNs") : 1;
    const int nOut = catch22_n_features(catch24);

    plhs[0] = mxCreateDoubleMatrix(nOut, 1, mxREAL);
    int masked[CATCH24_N_FEATURES];
    catch22_compute(mxGetPr(x), n, catch24, maskNaNs, mxGetPr(plhs[0]), masked);

    if (nlhs > 1) {
        plhs[1] = mxCreateLogicalMatrix(nOut, 1);
        mxLogical * m = mxGetLogicals(plhs[1]);
        for (int i = 0; i < nOut; i++) {
            m[i] = (mxLogical)masked[i];
        }
    }
}
