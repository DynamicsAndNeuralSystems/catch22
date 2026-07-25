/* ------------------------------------------------------------------------
 * Bit-exact, faster drop-in replacement for the recursive radix-2 FFT.
 *
 * Same signature, same twiddle table, same butterfly operands, same order
 * of operations inside every butterfly  ->  bit-identical output.
 *
 * MUST be compiled with contraction disabled, e.g.
 *     gcc -O3 -march=native -ffp-contract=off  (never -ffast-math)
 * otherwise a*c - b*d is fused into an FMA and the rounding changes.
 * (Verified: with gcc's default -ffp-contract=fast the outputs differ.)
 *
 * Anything the fast path cannot reproduce exactly is handed to a verbatim
 * copy of the original recursion:
 *   - size not a power of two (the original reads out of bounds there anyway)
 *   - any non-finite value in the input, where the inlined complex product
 *     could differ from the Annex-G product (__muldc3 / _Cmulcc) that the
 *     original calls.  Define FFT_ASSUME_FINITE to drop that O(n) screen.
 * ---------------------------------------------------------------------- */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#if __cplusplus
#   include <complex>
typedef std::complex< double > cplx;
#else
#   include <complex.h>
#if defined(__GNUC__) || defined(__GNUG__)
typedef double complex cplx;
#elif defined(_MSC_VER)
typedef _Dcomplex cplx;
#endif
#endif

#ifndef CMPLX
#define CMPLX(x, y) ((cplx)((double)(x) + _Imaginary_I * (double)(y)))
#endif

#include "helper_functions.h"

#pragma STDC FP_CONTRACT OFF

#if defined(__GNUC__) || defined(__GNUG__)
#  define RST __restrict__
#else
#  define RST __restrict
#endif

/* ===================== twiddle table (unchanged) ======================== */

void twiddles(cplx a[], int size)
{

    double PI = 3.14159265359;

    for (int i = 0; i < size; i++) {
        // cplx tmp = { 0, -PI * i / size };
        #if defined(__GNUC__) || defined(__GNUG__)
	    cplx tmp = 0.0  - PI * i / size * I;
        #elif defined(_MSC_VER)
	    cplx tmp = {0.0, -PI * i / size };
        #endif
        a[i] = cexp(tmp);
        //a[i] = cexp(-I * M_PI * i / size);
    }
}

/* ===================== exact fallback: original code ===================== */

static void _fft_rec(cplx a[], cplx out[], int size, int step, cplx tw[])
{
    if (step < size) {
        _fft_rec(out, a, size, step * 2, tw);
        _fft_rec(out + step, a + step, size, step * 2, tw);

        for (int i = 0; i < size; i += 2 * step) {
            cplx t = _Cmulcc(tw[i], out[i + step]);
            a[i / 2]          = _Caddcc(out[i], t);
            a[(i + size) / 2] = _Cminuscc(out[i], t);
        }
    }
}

static void fft_fallback(cplx a[], int size, cplx tw[])
{
    cplx *out = (cplx *)malloc((size_t)size * sizeof(cplx));
    if (!out) return;
    memcpy(out, a, (size_t)size * sizeof(cplx));
    _fft_rec(a, out, size, 1, tw);
    free(out);
}

/* ============================== fast path =============================== */

/* One decimation level, loop-interchanged so that the twiddle is loop
 * invariant and every stream is unit stride.
 *
 *   s  = 2^level,  nk = n/(2s) twiddle groups,  half = n/2
 *   dst[k*s + o]        = src[2ks + o] + tw[2ks]*src[2ks + s + o]
 *   dst[k*s + half + o] = src[2ks + o] - tw[2ks]*src[2ks + s + o]
 *
 * Exactly the butterflies the recursion performs (o was the recursion's
 * offset argument, k its loop index i/(2s)); only the visiting order of
 * independent butterflies differs, which cannot change any value.
 */
static void fft_pass(double *RST dst, const double *RST src,
                     const double *RST tw, size_t n, size_t s)
{
    const size_t nk   = n / (2 * s);
    const size_t half = n / 2;
    const size_t sd   = 2 * s;            /* block length in doubles */

    for (size_t k = 0; k < nk; ++k) {
        const double wr = tw[2 * sd * k];
        const double wi = tw[2 * sd * k + 1];

        const double *RST x = src + 2 * sd * k;   /* even half   */
        const double *RST y = x + sd;             /* odd half    */
        double *RST p = dst + sd * k;             /* "+" outputs */
        double *RST m = p + 2 * half;             /* "-" outputs */

        for (size_t o = 0; o < sd; o += 2) {
            const double yr = y[o], yi = y[o + 1];
            const double tr = wr * yr - wi * yi;
            const double ti = wr * yi + wi * yr;
            const double xr = x[o], xi = x[o + 1];
            p[o] = xr + tr;  p[o + 1] = xi + ti;
            m[o] = xr - tr;  m[o + 1] = xi - ti;
        }
    }
}

/* Top level (s = n/2, single twiddle group): source and destination indices
 * coincide, so it can run in place and save the memcpy of the original. */
static void fft_pass_top_inplace(double *RST a, const double *RST tw, size_t n)
{
    const double wr = tw[0], wi = tw[1];
    const size_t half = n;                /* n/2 complex = n doubles */

    for (size_t o = 0; o < half; o += 2) {
        const double yr = a[half + o], yi = a[half + o + 1];
        const double tr = wr * yr - wi * yi;
        const double ti = wr * yi + wi * yr;
        const double xr = a[o], xi = a[o + 1];
        a[o]            = xr + tr;  a[o + 1]            = xi + ti;
        a[half + o]     = xr - tr;  a[half + o + 1]     = xi - ti;
    }
}

/* The inlined complex product is bit-identical to the Annex-G product
 * (__muldc3 / _Cmulcc) unless both components come out NaN, which needs an
 * infinity in the operands.  Screen for that once, up front, and hand those
 * cases to the original code. */
static int all_finite(const double *p, size_t count)
{
    uint64_t bad = 0;
    for (size_t i = 0; i < count; ++i) {
        uint64_t u;
        memcpy(&u, &p[i], sizeof u);
        bad |= (uint64_t)((u & 0x7ff0000000000000ull) == 0x7ff0000000000000ull);
    }
    return !bad;
}

void fft(cplx a[], int size, cplx tw[])
{
    if (size < 2) return;                       /* original is a no-op here */

    const size_t n = (size_t)size;
    if (n & (n - 1)) { fft_fallback(a, size, tw); return; }   /* not 2^L */

    double       *A = (double *)a;
    const double *T = (const double *)tw;

#ifndef FFT_ASSUME_FINITE
    if (!all_finite(A, 2 * n)) { fft_fallback(a, size, tw); return; }
#endif

    unsigned L = 0;
    while ((n >> L) > 1) ++L;                   /* n == 1u << L */

    /* Levels run L-1 .. 0 and ping-pong between the two buffers; level 0 must
     * land in A.  If L is odd, run level L-1 in place so the parity works out
     * and the scratch buffer never needs to be initialised. */
    if (L == 1) { fft_pass_top_inplace(A, T, n); return; }   /* n == 2 */

    double *B = (double *)malloc(2 * n * sizeof(double));
    if (!B) { fft_fallback(a, size, tw); return; }            /* A untouched */

    size_t s     = n >> 1;
    unsigned rem = L;
    if (L & 1u) { fft_pass_top_inplace(A, T, n); s >>= 1; rem = L - 1; }

    double *src = A, *dst = B;
    for (; rem; --rem, s >>= 1) {
        fft_pass(dst, src, T, n, s);
        double *t = src; src = dst; dst = t;
    }
    /* rem was even, so the final write went to A */

    free(B);
}