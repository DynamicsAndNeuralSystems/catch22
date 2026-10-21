// Prints every feature of every series named on stdin (one path per line) as a hex
// float (bit for bit), with and without masking, one per line: file, masked, name,
// value. CI builds this against the base branch and the head of a pull request and
// diffs the two outputs, so a change in any feature value anywhere is seen on every
// platform. (Paths come on stdin because a thousand of them exceed Windows' command line.)
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "catch22.h"

int main(void)
{
    char path[4096];
    while (fgets(path, sizeof path, stdin)) {
        path[strcspn(path, "\r\n")] = '\0';
        if (!path[0]) continue;
        FILE * f = fopen(path, "r");
        if (!f) {
            fprintf(stderr, "cannot open %s\n", path);
            return 1;
        }
        int n = 0, cap = 1024;
        double * y = malloc(cap * sizeof * y);
        while (fscanf(f, "%lf", &y[n]) == 1) {
            if (++n == cap) y = realloc(y, (cap *= 2) * sizeof * y);
        }
        fclose(f);
        double out[CATCH24_N_FEATURES];
        for (int mask = 0; mask < 2; mask++) {
            catch22_compute(y, n, 1, mask, out, NULL);
            for (int i = 0; i < CATCH24_N_FEATURES; i++) {
                if (isnan(out[i])) { // sign and payload of a NaN vary by compiler
                    printf("%s\t%d\t%s\tnan\n", path, mask, catch22_feature_name(i));
                }
                else {
                    printf("%s\t%d\t%s\t%a\n", path, mask, catch22_feature_name(i), out[i]);
                }
            }
        }
        free(y);
    }
    return 0;
}
