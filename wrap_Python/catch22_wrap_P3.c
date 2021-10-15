#include <Python.h>

// include functions
#include "../C/CO_AutoCorr.h"
#include "../C/DN_HistogramMode_10.h"
#include "../C/DN_HistogramMode_5.h"
#include "../C/DN_Mean.h"
#include "../C/DN_Spread_Std.h"
#include "../C/DN_OutlierInclude.h"
#include "../C/FC_LocalSimple.h"
#include "../C/IN_AutoMutualInfoStats.h"
#include "../C/MD_hrv.h"
#include "../C/PD_PeriodicityWang.h"
#include "../C/SB_BinaryStats.h"
#include "../C/SB_CoarseGrain.h"
#include "../C/SB_MotifThree.h"
#include "../C/SB_TransitionMatrix.h"
#include "../C/SC_FluctAnal.h"
#include "../C/SP_Summaries.h"
#include "../C/butterworth.h"
#include "../C/fft.h"
#include "../C/helper_functions.h"
#include "../C/histcounts.h"
#include "../C/splinefit.h"
#include "../C/stats.h"


// ---------------------------------------------------------------------
// -------------------- Main wrapper function --------------------------
// ---------------------------------------------------------------------

// used for all functions to wrap via a function handle

static PyObject * python_wrapper_double(PyObject * args, double (*f) (const double*, const int), int normalize)
{

    PyObject * py_tuple;
    int n;
    double * c_array;

    // parse arguments
    if(!PyArg_ParseTuple(args, "O", &py_tuple)){
       return NULL;
    }

    if (PyList_Check(py_tuple)){

        n = (int) PyList_Size(py_tuple);

        // allocate space for input array (C data format)
        c_array= malloc(n*sizeof(double));

        // write input array to c array
        int i;
        for (i=0; i<n; i++) {
            c_array[i] = (double) PyFloat_AsDouble(PyList_GetItem(py_tuple, i));
        }

    }
    else
    {
        if (PyTuple_Check(py_tuple)){

            n = (int) PyTuple_Size(py_tuple);

            // allocate space for input array (C data format)
            c_array= malloc(n*sizeof(double));

            // write input array to c array
            int i;
            for (i=0; i<n; i++) {
                c_array[i] = (double) PyFloat_AsDouble(PyTuple_GetItem(py_tuple, i));
            }

        }
        else{
            return NULL;
        }
    }

    // result variables
    double result;
    PyObject * ret;

    // copy input vector (just in case it is altered by the function)
    double * copy = malloc(n * sizeof * copy);    
    memcpy(copy, c_array, n * sizeof * copy);

    if (normalize){
        
        double * y_zscored = malloc(n * sizeof * y_zscored);
        zscore_norm2(copy, n, y_zscored);

        result = f(y_zscored, n);

        free(y_zscored);
    } 
    else {
        result = f(copy, n);
    }   

    free(c_array);
    free(copy);

    // build the resulting string into a Python object.
    ret = Py_BuildValue("d", result);

    return ret;
};

static PyObject * python_wrapper_int(PyObject * args, int (*f) (const double*, const int), int normalize)
{

    PyObject * py_tuple;
    int n;
    double * c_array;

    // parse arguments
    if(!PyArg_ParseTuple(args, "O", &py_tuple)){
       return NULL;
    }

    if (PyList_Check(py_tuple)){

        n = (int) PyList_Size(py_tuple);

        // allocate space for input array (C data format)
        c_array= malloc(n*sizeof(double));

        // write input array to c array
        int i;
        for (i=0; i<n; i++) {
            c_array[i] = (double) PyFloat_AsDouble(PyList_GetItem(py_tuple, i));
        }

    }
    else
    {
        if (PyTuple_Check(py_tuple)){

            n = (int) PyTuple_Size(py_tuple);

            // allocate space for input array (C data format)
            c_array= malloc(n*sizeof(double));

            // write input array to c array
            int i;
            for (i=0; i<n; i++) {
                c_array[i] = (double) PyFloat_AsDouble(PyTuple_GetItem(py_tuple, i));
            }

        }
        else{
            return NULL;
        }
    }

    // result variables
    int result;
    PyObject * ret;

    // copy input vector (just in case it is altered by the function)
    double * copy = malloc(n * sizeof * copy);    
    memcpy(copy, c_array, n * sizeof * copy);

    if (normalize){
        
        double * y_zscored = malloc(n * sizeof * y_zscored);
        zscore_norm2(copy, n, y_zscored);

        result = f(y_zscored, n);

        free(y_zscored);
    } 
    else {
        result = f(copy, n);
    }   

    free(copy);
    free(c_array);

    // build the resulting string into a Python object.
    ret = Py_BuildValue("n", result);

    return ret;
};


// ---------------------------------------------------------------------
// ----------------------- Functions to wrap ---------------------------
// ---------------------------------------------------------------------

static PyObject * DN_HistogramMode_5_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_double(args, &DN_HistogramMode_5, 1);
}

static PyObject * DN_HistogramMode_10_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_double(args, &DN_HistogramMode_10, 1);
}

static PyObject * CO_f1ecac_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_int(args, &CO_f1ecac, 1);
}

static PyObject * CO_FirstMin_ac_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_int(args, &CO_FirstMin_ac, 1);
}

static PyObject * CO_HistogramAMI_even_2_5_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_double(args, &CO_HistogramAMI_even_2_5, 1);
}

static PyObject * CO_trev_1_num_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_double(args, &CO_trev_1_num, 1);
}

static PyObject * MD_hrv_classic_pnn40_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_double(args, &MD_hrv_classic_pnn40, 1);
}

static PyObject * SB_BinaryStats_mean_longstretch1_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_double(args, &SB_BinaryStats_mean_longstretch1, 1);
}

static PyObject * SB_TransitionMatrix_3ac_sumdiagcov_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_double(args, &SB_TransitionMatrix_3ac_sumdiagcov, 1);
}

static PyObject * PD_PeriodicityWang_th0_01_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_int(args, &PD_PeriodicityWang_th0_01, 1);
}

static PyObject * CO_Embed2_Dist_tau_d_expfit_meandiff_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_double(args, &CO_Embed2_Dist_tau_d_expfit_meandiff, 1);
}

static PyObject * IN_AutoMutualInfoStats_40_gaussian_fmmi_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_double(args, &IN_AutoMutualInfoStats_40_gaussian_fmmi, 1);
}

static PyObject * FC_LocalSimple_mean1_tauresrat_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_double(args, &FC_LocalSimple_mean1_tauresrat, 1);
}

static PyObject * DN_OutlierInclude_p_001_mdrmd_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_double(args, &DN_OutlierInclude_p_001_mdrmd, 1);
}

static PyObject * DN_OutlierInclude_n_001_mdrmd_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_double(args, &DN_OutlierInclude_n_001_mdrmd, 1);
}

static PyObject * SP_Summaries_welch_rect_area_5_1_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_double(args, &SP_Summaries_welch_rect_area_5_1, 1);
}

static PyObject * SB_BinaryStats_diff_longstretch0_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_double(args, &SB_BinaryStats_diff_longstretch0, 1);
}

static PyObject * SB_MotifThree_quantile_hh_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_double(args, &SB_MotifThree_quantile_hh, 1);
}

static PyObject * SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_double(args, &SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1, 1);
}

static PyObject * SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_double(args, &SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1, 1);
}

static PyObject * SP_Summaries_welch_rect_centroid_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_double(args, &SP_Summaries_welch_rect_centroid, 1);
}

static PyObject * FC_LocalSimple_mean3_stderr_wrapper(PyObject * self, PyObject * args)
{
    return python_wrapper_double(args, &FC_LocalSimple_mean3_stderr, 1);
}

// ---------------------------------------------------------------------
// ------------------------ Module creation ----------------------------
// ---------------------------------------------------------------------

// register within a module
static PyMethodDef catch22Methods[] = {
    { "DN_HistogramMode_5", DN_HistogramMode_5_wrapper, METH_VARARGS, "bla" },
    { "DN_HistogramMode_10", DN_HistogramMode_10_wrapper, METH_VARARGS, "bla" },
    { "CO_f1ecac", CO_f1ecac_wrapper, METH_VARARGS, "bla" },
    { "CO_FirstMin_ac", CO_FirstMin_ac_wrapper, METH_VARARGS, "bla" },
    { "CO_HistogramAMI_even_2_5", CO_HistogramAMI_even_2_5_wrapper, METH_VARARGS, "bla" },
    { "CO_trev_1_num", CO_trev_1_num_wrapper, METH_VARARGS, "bla" },
    { "MD_hrv_classic_pnn40", MD_hrv_classic_pnn40_wrapper, METH_VARARGS, "bla" },
    { "SB_BinaryStats_mean_longstretch1", SB_BinaryStats_mean_longstretch1_wrapper, METH_VARARGS, "bla" },
    { "SB_TransitionMatrix_3ac_sumdiagcov", SB_TransitionMatrix_3ac_sumdiagcov_wrapper, METH_VARARGS, "bla" },
    { "PD_PeriodicityWang_th0_01", PD_PeriodicityWang_th0_01_wrapper, METH_VARARGS, "bla" },
    { "CO_Embed2_Dist_tau_d_expfit_meandiff", CO_Embed2_Dist_tau_d_expfit_meandiff_wrapper, METH_VARARGS, "bla" },
    { "IN_AutoMutualInfoStats_40_gaussian_fmmi", IN_AutoMutualInfoStats_40_gaussian_fmmi_wrapper, METH_VARARGS, "bla" },
    { "FC_LocalSimple_mean1_tauresrat", FC_LocalSimple_mean1_tauresrat_wrapper, METH_VARARGS, "bla" },
    { "DN_OutlierInclude_p_001_mdrmd", DN_OutlierInclude_p_001_mdrmd_wrapper, METH_VARARGS, "bla" },
    { "DN_OutlierInclude_n_001_mdrmd", DN_OutlierInclude_n_001_mdrmd_wrapper, METH_VARARGS, "bla" },
    { "SP_Summaries_welch_rect_area_5_1", SP_Summaries_welch_rect_area_5_1_wrapper, METH_VARARGS, "bla" },
    { "SB_BinaryStats_diff_longstretch0", SB_BinaryStats_diff_longstretch0_wrapper, METH_VARARGS, "bla" },
    { "SB_MotifThree_quantile_hh", SB_MotifThree_quantile_hh_wrapper, METH_VARARGS, "bla" },
    { "SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1", SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1_wrapper, METH_VARARGS, "bla" },
    { "SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1", SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1_wrapper, METH_VARARGS, "bla" },
    { "SP_Summaries_welch_rect_centroid", SP_Summaries_welch_rect_centroid_wrapper, METH_VARARGS, "bla" },
    { "FC_LocalSimple_mean3_stderr", FC_LocalSimple_mean3_stderr_wrapper, METH_VARARGS, "bla" },
    { NULL, NULL, 0, NULL }
};


static struct PyModuleDef catch22_Cmodule = {
    PyModuleDef_HEAD_INIT,
    "catch22_C",   /* name of module */
    "22 Canonical Time-series CHaracteristcs", /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    catch22Methods,
};

PyMODINIT_FUNC PyInit_catch22_C(void)
{
    return PyModule_Create(&catch22_Cmodule);
}
