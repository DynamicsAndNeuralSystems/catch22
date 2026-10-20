// Command-line tool: compute catch22 (or catch24) for a time series stored in a
// text file (one value per line). Usage:
//
//   run_features [--catch22 | --catch24] [--no-mask] <infile> [<outfile>]
//
// Each output line is "value, feature name, computation time (ms)". Without
// --catch22/--catch24, the tool asks interactively. With --no-mask, features
// that cannot be estimated are reported as NaN instead of being masked.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "catch22.h"

void run_features(const double y[], const int size, FILE * outfile, const int catch24, const int mask_nans)
{
    const int nOut = catch22_n_features(catch24);
    double * out = malloc(nOut * sizeof * out);
    double * time_ms = malloc(nOut * sizeof * time_ms);

    if (catch22_compute_timed(y, size, catch24, mask_nans, out, NULL, time_ms) == CATCH22_BAD_INPUT) {
        fprintf(stderr, "Warning: the time series is empty or contains NaN/Inf; all features are NaN.\n");
    }
    for (int i = 0; i < nOut; i++) {
        fprintf(outfile, "%.14f, %s, %f\n", out[i], catch22_feature_name(i), time_ms[i]);
    }
    fprintf(outfile, "\n");

    free(out);
    free(time_ms);
}

void print_help(char *argv[], char msg[])
{
    if (strlen(msg) > 0) {
        fprintf(stdout, "ERROR: %s\n", msg);
    }
    fprintf(stdout, "Usage is %s [--catch22 | --catch24] [--no-mask] <infile> [<outfile>]\n", argv[0]);
    fprintf(stdout, "\n\tSpecifying outfile is optional, by default it is stdout\n");
    fprintf(stdout, "\n       %s --list\n\n\tList the 24 features (tab-separated: name, short name, v0.x name, hctsa call)\n", argv[0]);
    exit(1);
}

int main(int argc, char * argv[])
{
    int catch24 = -1; // -1: ask
    int mask_nans = 1;
    const char * inName = NULL;
    const char * outName = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--catch24") == 0) {
            catch24 = 1;
        }
        else if (strcmp(argv[i], "--catch22") == 0) {
            catch24 = 0;
        }
        else if (strcmp(argv[i], "--no-mask") == 0) {
            mask_nans = 0;
        }
        else if (strcmp(argv[i], "--list") == 0) {
            fprintf(stdout, "name\tshort_name\tlegacy_name\thctsa_call\n");
            for (int k = 0; k < catch22_n_features(1); k++) {
                fprintf(stdout, "%s\t%s\t%s\t%s\n", catch22_feature_name(k), catch22_short_name(k),
                        catch22_legacy_name(k), catch22_hctsa_call(k));
            }
            return 0;
        }
        else if (inName == NULL) {
            inName = argv[i];
        }
        else if (outName == NULL) {
            outName = argv[i];
        }
        else {
            print_help(argv, "Too many arguments\n");
            return 1;
        }
    }
    if (inName == NULL) {
        print_help(argv, "");
        return 1;
    }

    FILE * infile = fopen(inName, "r");
    if (infile == NULL) {
        print_help(argv, "Can't open input file\n");
        return 1;
    }
    FILE * outfile = stdout;
    if (outName != NULL && (outfile = fopen(outName, "w")) == NULL) {
        fclose(infile);
        print_help(argv, "Can't open output file\n");
        return 1;
    }

    // read the time series (stops at the first entry that is not a number)
    int arraySize = 50;
    int size = 0;
    double * y = malloc(arraySize * sizeof * y);
    double value;
    while (fscanf(infile, "%lf", &value) == 1) {
        if (size == arraySize) {
            arraySize *= 2;
            y = realloc(y, arraySize * sizeof * y);
        }
        y[size++] = value;
    }
    fclose(infile);

    if (catch24 < 0) {
        printf("Do you want to run catch24? Enter 0 for catch22 or 1 for catch24.");
        if (scanf("%d", &catch24) != 1) {
            catch24 = 0; // default to catch22 if no valid answer is given
        }
    }

    run_features(y, size, outfile, catch24 == 1, mask_nans);

    if (outfile != stdout) {
        fclose(outfile);
    }
    free(y);
    return 0;
}
