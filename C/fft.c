#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#ifndef CMPLX
#define CMPLX(x, y) ((cplx)((double)(x) + _Imaginary_I * (double)(y)))
#endif

typedef double complex cplx;

void twiddles(cplx a[], int size)
{
    for (int i = 0; i < size; i++) {
        a[i] = cexp(-I * M_PI * i / size);
    }
}

static void _fft(cplx a[], cplx out[], int size, int step, cplx tw[])
{   
    if (step < size) {
        _fft(out, a, size, step * 2, tw);
        _fft(out + step, a + step, size, step * 2, tw);

        for (int i = 0; i < size; i += 2 * step) {
            cplx t = tw[i] * out[i + step];
            a[i / 2] = out[i] + t;
            a[(i + size) / 2] = out[i] - t;
        }
    }
}

void fft(cplx a[], int size, cplx tw[])
{
    cplx * out = malloc(size * sizeof(cplx));
    memcpy(out, a, size * sizeof(cplx));
    _fft(a, out, size, 1, tw);
    free(out);
}
