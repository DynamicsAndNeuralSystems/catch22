#include <Rcpp.h>

// include functions
extern "C" {
#include "CO_AutoCorr.h"
#include "DN_HistogramMode_10.h"
#include "DN_HistogramMode_5.h"
#include "DN_Mean.h"
#include "DN_Spread_Std.h"
#include "DN_OutlierInclude.h"
#include "FC_LocalSimple.h"
#include "IN_AutoMutualInfoStats.h"
#include "MD_hrv.h"
#include "PD_PeriodicityWang.h"
#include "SB_BinaryStats.h"
#include "SB_CoarseGrain.h"
#include "SB_MotifThree.h"
#include "SB_TransitionMatrix.h"
#include "SC_FluctAnal.h"
#include "SP_Summaries.h"
#include "butterworth.h"
#include "fft.h"
#include "helper_functions.h"
#include "histcounts.h"
#include "splinefit.h"
#include "stats.h"
}

using namespace Rcpp;

// universal wrapper for a function that takes a double array and its length
// and outputs a scalar double 
NumericVector R_wrapper_double(NumericVector x, double (*f) (const double*, const int), int normalize) {
  
  int n = x.size();
  double * arrayC = new double[n]; 
  double out;
  
  int i;
  for (i=0; i<n; i++){
    arrayC[i] = x[i];
  }
  
  if (normalize){
        
        double * y_zscored = new double[n];
        
        zscore_norm2(arrayC, n, y_zscored);

        out = f(y_zscored, n);

        // free(y_zscored);
    } 
    else {
        out = f(arrayC, n);
    }

    // free(arrayC);

  NumericVector outVec = NumericVector::create(out);

  return outVec;
  
};

// universal wrapper for a function that takes a double array and its length
// and outputs a scalar double 
NumericVector R_wrapper_int(NumericVector x, int (*f) (const double*, const int), int normalize) {
  
  int n = x.size();
  double * arrayC = new double[n]; 
  int out;
  
  int i;
  for (i=0; i<n; i++){
    arrayC[i] = x[i];
  }
  
  if (normalize){
        
        double * y_zscored = new double[n];
        
        zscore_norm2(arrayC, n, y_zscored);

        out = f(y_zscored, n);

        // free(y_zscored);
    } 
    else {
        out = f(arrayC, n);
    }

    // free(arrayC);

  NumericVector outVec = NumericVector::create(out);

  return outVec;
  
};

//-------------------------------------------------------------------------
//----------------------- Feature functions -------------------------------
//-------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector DN_HistogramMode_5(NumericVector x)
{
  return R_wrapper_double(x, &DN_HistogramMode_5, 1);
}

// [[Rcpp::export]]
NumericVector DN_HistogramMode_10(NumericVector x)
{
  return R_wrapper_double(x, &DN_HistogramMode_10, 1);
}

// [[Rcpp::export]]
NumericVector CO_f1ecac(NumericVector x)
{
  return R_wrapper_int(x, &CO_f1ecac, 1);
}

// [[Rcpp::export]]
NumericVector CO_FirstMin_ac(NumericVector x)
{
  return R_wrapper_int(x, &CO_FirstMin_ac, 1);
}

// [[Rcpp::export]]
NumericVector CO_HistogramAMI_even_2_5(NumericVector x)
{
  return R_wrapper_double(x, &CO_HistogramAMI_even_2_5, 1);
}

// [[Rcpp::export]]
NumericVector CO_trev_1_num(NumericVector x)
{
  return R_wrapper_double(x, &CO_trev_1_num, 1);
}

// [[Rcpp::export]]
NumericVector MD_hrv_classic_pnn40(NumericVector x)
{
  return R_wrapper_double(x, &MD_hrv_classic_pnn40, 1);
}

// [[Rcpp::export]]
NumericVector SB_BinaryStats_mean_longstretch1(NumericVector x)
{
  return R_wrapper_double(x, &SB_BinaryStats_mean_longstretch1, 1);
}

// [[Rcpp::export]]
NumericVector SB_TransitionMatrix_3ac_sumdiagcov(NumericVector x)
{
  return R_wrapper_double(x, &SB_TransitionMatrix_3ac_sumdiagcov, 1);
}

// [[Rcpp::export]]
NumericVector PD_PeriodicityWang_th0_01(NumericVector x)
{
  return R_wrapper_int(x, &PD_PeriodicityWang_th0_01, 1);
}

// [[Rcpp::export]]
NumericVector CO_Embed2_Dist_tau_d_expfit_meandiff(NumericVector x)
{
  return R_wrapper_double(x, &CO_Embed2_Dist_tau_d_expfit_meandiff, 1);
}

// [[Rcpp::export]]
NumericVector IN_AutoMutualInfoStats_40_gaussian_fmmi(NumericVector x)
{
  return R_wrapper_double(x, &IN_AutoMutualInfoStats_40_gaussian_fmmi, 1);
}

// [[Rcpp::export]]
NumericVector FC_LocalSimple_mean1_tauresrat(NumericVector x)
{
  return R_wrapper_double(x, &FC_LocalSimple_mean1_tauresrat, 1);
}

// [[Rcpp::export]]
NumericVector DN_OutlierInclude_p_001_mdrmd(NumericVector x)
{
  return R_wrapper_double(x, &DN_OutlierInclude_p_001_mdrmd, 1);
}

// [[Rcpp::export]]
NumericVector DN_OutlierInclude_n_001_mdrmd(NumericVector x)
{
  return R_wrapper_double(x, &DN_OutlierInclude_n_001_mdrmd, 1);
}

// [[Rcpp::export]]
NumericVector SP_Summaries_welch_rect_area_5_1(NumericVector x)
{
  return R_wrapper_double(x, &SP_Summaries_welch_rect_area_5_1, 1);
}

// [[Rcpp::export]]
NumericVector SB_BinaryStats_diff_longstretch0(NumericVector x)
{
  return R_wrapper_double(x, &SB_BinaryStats_diff_longstretch0, 1);
}

// [[Rcpp::export]]
NumericVector SB_MotifThree_quantile_hh(NumericVector x)
{
  return R_wrapper_double(x, &SB_MotifThree_quantile_hh, 1);
}

// [[Rcpp::export]]
NumericVector SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1(NumericVector x)
{
  return R_wrapper_double(x, &SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1, 1);
}

// [[Rcpp::export]]
NumericVector SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1(NumericVector x)
{
  return R_wrapper_double(x, &SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1, 1);
}

// [[Rcpp::export]]
NumericVector SP_Summaries_welch_rect_centroid(NumericVector x)
{
  return R_wrapper_double(x, &SP_Summaries_welch_rect_centroid, 1);
}

// [[Rcpp::export]]
NumericVector FC_LocalSimple_mean3_stderr(NumericVector x)
{
  return R_wrapper_double(x, &FC_LocalSimple_mean3_stderr, 1);
}

// [[Rcpp::export]]
NumericVector DN_Mean(NumericVector x)
{
  return R_wrapper_double(x, &DN_Mean, 1);
}

// [[Rcpp::export]]
NumericVector DN_Spread_Std(NumericVector x)
{
  return R_wrapper_double(x, &DN_Spread_Std, 1);
}
