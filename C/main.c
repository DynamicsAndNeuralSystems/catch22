/* Include files */
#include "main.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
//#include <dirent.h>

#include "DN_HistogramMode_5.h"
#include "DN_HistogramMode_10.h"
#include "DN_Mean.h"
#include "DN_Spread_Std.h"
#include "CO_AutoCorr.h"
#include "DN_OutlierInclude.h"
#include "FC_LocalSimple.h"
#include "IN_AutoMutualInfoStats.h"
#include "MD_hrv.h"
#include "SB_BinaryStats.h"
#include "SB_MotifThree.h"
#include "SC_FluctAnal.h"
#include "SP_Summaries.h"
#include "SB_TransitionMatrix.h"
#include "PD_PeriodicityWang.h"

#include "stats.h"

// check if data qualifies to be caught22
int quality_check(const double y[], const int size)
{
    int minSize = 10;

    if(size < minSize)
    {
        return 1;
    }
    for(int i = 0; i < size; i++)
    {
        double val = y[i];
        if(val == INFINITY || -val == INFINITY)
        {
            return 2;
        }
        if(isnan(val))
        {
            return 3;
        }
    }
    return 0;
}

void run_features(double y[], int size, FILE * outfile, bool catch24)
{
    int quality = quality_check(y, size);
    if(quality != 0)
    {
        fprintf(stdout, "Time series quality test not passed (code %i).\n", quality);
        return;
    }

    printf("Time series lenght %d.\n", size);

    double * y_zscored = malloc(size * sizeof * y_zscored);

    // variables to keep time
    clock_t begin;
    double timeTaken;
    double tot_time = 0.0;

    // output
    double result;

    // z-score first for all.
    zscore_norm2(y, size, y_zscored);

    // GOOD
    begin = clock();
    result = DN_OutlierInclude_n_001_mdrmd(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "DN_OutlierInclude_n_001_mdrmd", timeTaken);

    // GOOD
    begin = clock();
    result = DN_OutlierInclude_p_001_mdrmd(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "DN_OutlierInclude_p_001_mdrmd", timeTaken);
  
    // GOOD
    begin = clock();
    result = DN_HistogramMode_5(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "DN_HistogramMode_5", timeTaken);

    // GOOD
    begin = clock();
    result = DN_HistogramMode_10(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "DN_HistogramMode_10", timeTaken);

    //GOOD
    begin = clock();
    result = CO_Embed2_Dist_tau_d_expfit_meandiff(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "CO_Embed2_Dist_tau_d_expfit_meandiff", timeTaken);

    //GOOD (memory leak?)
    begin = clock();
    result = CO_f1ecac(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "CO_f1ecac", timeTaken);

    //GOOD
    begin = clock();
    result = CO_FirstMin_ac(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "CO_FirstMin_ac", timeTaken);

    // GOOD (memory leak?)
    begin = clock();
    result = CO_HistogramAMI_even_2_5(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "CO_HistogramAMI_even_2_5", timeTaken);

    // GOOD
    begin = clock();
    result = CO_trev_1_num(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "CO_trev_1_num", timeTaken);

    //GOOD
    begin = clock();
    result = FC_LocalSimple_mean1_tauresrat(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "FC_LocalSimple_mean1_tauresrat", timeTaken);

    //GOOD
    begin = clock();
    result = FC_LocalSimple_mean3_stderr(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "FC_LocalSimple_mean3_stderr", timeTaken);

    //GOOD (memory leak?)
    begin = clock();
    result = IN_AutoMutualInfoStats_40_gaussian_fmmi(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "IN_AutoMutualInfoStats_40_gaussian_fmmi", timeTaken);

    //GOOD
    begin = clock();
    result = MD_hrv_classic_pnn40(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "MD_hrv_classic_pnn40", timeTaken);

    //GOOD
    begin = clock();
    result = SB_BinaryStats_diff_longstretch0(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "SB_BinaryStats_diff_longstretch0", timeTaken);

    //GOOD
    begin = clock();
    result = SB_BinaryStats_mean_longstretch1(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "SB_BinaryStats_mean_longstretch1", timeTaken);

    //GOOD (memory leak?)
    begin = clock();
    result = SB_MotifThree_quantile_hh(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "SB_MotifThree_quantile_hh", timeTaken);

    //GOOD (memory leak?)
    begin = clock();
    result = SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1", timeTaken);

    //GOOD
    begin = clock();
    result = SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1", timeTaken);

    //GOOD
    begin = clock();
    result = SP_Summaries_welch_rect_area_5_1(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "SP_Summaries_welch_rect_area_5_1", timeTaken);

    //GOOD
    begin = clock();
    result = SP_Summaries_welch_rect_centroid(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "SP_Summaries_welch_rect_centroid", timeTaken);

    //OK, BUT filt in Butterworth sometimes diverges, now removed alltogether, let's see results.
    begin = clock();
    result = SB_TransitionMatrix_3ac_sumdiagcov(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "SB_TransitionMatrix_3ac_sumdiagcov", timeTaken);

    // GOOD
    begin = clock();
    result = PD_PeriodicityWang_th0_01(y_zscored, size);
    timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
    tot_time += timeTaken;
    fprintf(outfile, "%.14f, %s, %f\n", result, "PD_PeriodicityWang_th0_01", timeTaken);

    if (catch24) {
        
        // GOOD
        begin = clock();
        result = DN_Mean(y, size);
        timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
        tot_time += timeTaken;
        fprintf(outfile, "%.14f, %s, %f\n", result, "DN_Mean", timeTaken);

        // GOOD
        begin = clock();
        result = DN_Spread_Std(y, size);
        timeTaken = (double)(clock()-begin)*1000/CLOCKS_PER_SEC;
        tot_time += timeTaken;
        fprintf(outfile, "%.14f, %s, %f\n", result, "DN_Spread_Std", timeTaken);
    } else {

    }
    
    printf("Total time taken: %f\n", tot_time);
    fprintf(outfile, "\n");

    free(y_zscored);
}

void print_help(char *argv[], char msg[])
{
    if (strlen(msg) > 0) {
        fprintf(stdout, "ERROR: %s\n", msg);
    }
    fprintf(stdout, "Usage is %s <infile> <outfile>\n", argv[0]);
    fprintf(stdout, "\n\tSpecifying outfile is optional, by default it is stdout\n");
    // fprintf(stdout, "\tOutput order is:\n%s\n", HEADER);
    exit(1);
}

// memory leak check; use with valgrind.
#if 0
int main(int argc, char * argv[])
{
    double * y = malloc(1000 * sizeof(double));

    srand(42);
    for (int i = 0; i < 1000; ++i) {
        y[i] = rand() % RAND_MAX;
    }
    run_features(y, 1000, stdout);
    free(y);
}
#endif

#if 1
int main(int argc, char * argv[])
{
    FILE * infile, * outfile;
    int array_size;
    double * y;
    int size;
    double value;
    // DIR *d;
    struct dirent *dir;


    switch (argc) {
        case 1:
            print_help(argv, "");
            break;
        case 2:
            if ((infile = fopen(argv[1], "r")) == NULL) {
                print_help(argv, "Can't open input file\n");
            }
            outfile = stdout;
            break;
        case 3:
            if ((infile = fopen(argv[1], "r")) == NULL) {
                print_help(argv, "Can't open input file\n");
            }
            if ((outfile = fopen(argv[2], "w")) == NULL) {
                print_help(argv, "Can't open output file\n");
            }
            break;
    }

    /*
    // debug: fix these.
    infile = fopen("/Users/carl/PycharmProjects/catch22/C/timeSeries/tsid0244.txt", "r");
    outfile = stdout;
     */

    // fprintf(outfile, "%s", HEADER);
    array_size = 50;
    size = 0;
    y = malloc(array_size * sizeof *y);

    while (fscanf(infile, "%lf", &value) != EOF) {
        if (size == array_size) {
            y = realloc(y, 2 * array_size * sizeof *y);
            array_size *= 2;
        }
        y[size++] = value;
    }
    fclose(infile);
    y = realloc(y, size * sizeof *y);
    //printf("size=%i\n", size);

    // catch24 specification

    int catch24;
    printf("Do you want to run catch24? Enter 0 for catch22 or 1 for catch24.");
    scanf("%d", &catch24);

    if (catch24 == 1) {
        run_features(y, size, outfile, true);
    } else {
        run_features(y, size, outfile, false);
    }

    fclose(outfile);
    free(y);

    return 0;
}
#endif

#if 0
int main(int argc, char * argv[])
{
  (void)argc;
  (void)argv;

    /*
    // generate some data
    const int size = 31; // 211;

    double y[size];
    int i;
    double sinIn=0;
    for(i=0; i<size; i++)
    {
        sinIn = (double)(i-10)/10*6;
        y[i] = sin(sinIn);
        printf("%i: sinIn=%1.3f, sin=%1.3f\n", i, sinIn, y[i]);
    }*/

    // open a certain file
    FILE * infile;
    infile = fopen("C:\\Users\\Carl\\Documents\\catch22-master\\testData\\test.txt", "r");
    int array_size = 15000;
    double * y = malloc(array_size * sizeof(double));
    int size = 0;
    double value = 0;
    while (fscanf(infile, "%lf", &value) != EOF) {
        y[size++] = value;
    }

    // first, z-score.
    zscore_norm(y, size);

    double result;
  
    result = DN_OutlierInclude_n_001_mdrmd(y, size);
    printf("DN_OutlierInclude_n_001_mdrmd: %1.5f\n", result);
    result = DN_OutlierInclude_p_001_mdrmd(y, size);
    printf("DN_OutlierInclude_p_001_mdrmd: %1.5f\n", result);
    result = DN_HistogramMode_5(y, size);
    printf("DN_HistogramMode_5: %1.3f\n", result);
    result = DN_HistogramMode_10(y, size);
    printf("DN_HistogramMode_10: %1.3f\n", result);
    result = SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1(y, size);
    printf("SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1: %1.5f\n", result);
    result = SB_TransitionMatrix_3ac_sumdiagcov(y, size);
    printf("SB_TransitionMatrix_3ac_sumdiagcov: %1.5f\n", result);
    result = FC_LocalSimple_mean1_tauresrat(y, size);
    printf("FC_LocalSimple_mean1_tauresrat: %1.5f\n", result);
    result = SB_MotifThree_quantile_hh(y, size);
    printf("SB_MotifThree_quantile_hh: %1.5f\n", result);
    result = CO_HistogramAMI_even_2_5(y, size);
    printf("CO_HistogramAMI_even_2_5: %1.3f\n", result);
    result = CO_Embed2_Dist_tau_d_expfit_meandiff(y, size);
    printf("CO_Embed2_Dist_tau_d_expfit_meandiff: %1.3f\n", result);
    result = SB_BinaryStats_diff_longstretch0(y, size);
    printf("SB_BinaryStats_diff_longstretch0: %1.5f\n", result);
    result = MD_hrv_classic_pnn40(y, size);
    printf("MD_hrv_classic_pnn40: %1.5f\n", result);
    result = SB_BinaryStats_mean_longstretch1(y, size);
    printf("SB_BinaryStats_mean_longstretch1: %1.5f\n", result);
    result = FC_LocalSimple_mean3_stderr(y, size);
    printf("FC_LocalSimple_mean3_stderr: %1.5f\n", result);
    result = SP_Summaries_welch_rect_area_5_1(y, size);
    printf("SP_Summaries_welch_rect_area_5_1: %1.5f\n", result);
    result = SP_Summaries_welch_rect_centroid(y, size);
    printf("SP_Summaries_welch_rect_centroid: %1.5f\n", result);
    result = CO_f1ecac(y, size);
    printf("CO_f1ecac: %1.f\n", result);
    result = CO_FirstMin_ac(y, size);
    printf("CO_FirstMin_ac: %1.f\n", result);
    result = IN_AutoMutualInfoStats_40_gaussian_fmmi(y, size);
    printf("IN_AutoMutualInfoStats_40_gaussian_fmmi: %1.5f\n", result);
    result = PD_PeriodicityWang_th0_01(y, size);
    printf("PD_PeriodicityWang_th0_01: %1.f\n", result);
    result = SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1(y, size);
    printf("SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1: %1.5f\n", result);
    result = CO_trev_1_num(y, size);
    printf("CO_trev_1_num: %1.5f\n", result);

    int catch24;
    printf("Do you want to run catch24? Enter 0 for catch22 or 1 for catch24.");
    scanf("%d", &catch24);

    if (catch24 == 1) {
        result = DN_Mean(y, size);
        printf("DN_Mean: %1.5f\n", result);
        result = DN_Spread_Std(y, size);
        printf("DN_Spread_Std: %1.5f\n", result);
    } else {
    }

  return 0;
}
#endif
