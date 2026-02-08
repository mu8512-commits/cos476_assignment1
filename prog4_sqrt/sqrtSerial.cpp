#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <cstddef>


void sqrtSerial(int N,
                float initialGuess,
                float values[],
                float output[])
{

    static const float kThreshold = 0.00001f;

    for (int i=0; i<N; i++) {

        float x = values[i];
        float guess = initialGuess;

        float error = fabs(guess * guess * x - 1.f);

        while (error > kThreshold) {
            guess = (3.f * guess - x * guess * guess * guess) * 0.5f;
            error = fabs(guess * guess * x - 1.f);
        }

        output[i] = x * guess;
    }
}

// Helper: absolute value for __m256 (clear sign bit)
static inline __m256 mm256_abs_ps(__m256 x) {
    const __m256 signmask = _mm256_set1_ps(-0.0f); // sign bit set
    return _mm256_andnot_ps(signmask, x);
}

// Computes sqrt-like output for N floats in `values` into `output`.
// (This AVX2 version uses hardware rsqrt + 2 Newton refinements;
void sqrt_avx2(int N, float startGuess, float* values, float* output) {
    const int W = 8;
    int i = 0;

    const __m256 one   = _mm256_set1_ps(1.0f);
    const __m256 half  = _mm256_set1_ps(0.5f);
    const __m256 eps   = _mm256_set1_ps(1e-4f);   // <-- MUST match your serial epsilon
    const __m256 zero  = _mm256_set1_ps(0.0f);

    // Safety cap to avoid infinite loops if something goes wrong (match serial if it has one)
    const int MAX_ITERS = 100;

    for (; i + W <= N; i += W) {
        __m256 S = _mm256_loadu_ps(values + i);

        // x initialized from startGuess (same for all lanes)
        __m256 x = _mm256_set1_ps(startGuess);

        // Active mask: all lanes start active
        __m256 active_ps = _mm256_castsi256_ps(_mm256_set1_epi32(-1));

        for (int it = 0; it < MAX_ITERS; it++) {
            // Compute f(x) = 1/x^2 - S
            __m256 x2   = _mm256_mul_ps(x, x);
            __m256 invx2 = _mm256_div_ps(one, x2);
            __m256 f    = _mm256_sub_ps(invx2, S);

            // Check convergence: |f| > eps  => still active
            __m256 absf = mm256_abs_ps(f);
            __m256 need = _mm256_cmp_ps(absf, eps, _CMP_GT_OQ); // mask where not converged

            // If no lanes need work, break
            int need_mask = _mm256_movemask_ps(need);
            if (need_mask == 0) break;

            // Newton step:
            // f'(x) = d/dx (x^{-2} - S) = -2 x^{-3}
            // x_{new} = x - f/f' = x + 0.5 * f * x^3
            __m256 x3    = _mm256_mul_ps(x2, x);
            __m256 delta = _mm256_mul_ps(half, _mm256_mul_ps(f, x3));
            __m256 xnew  = _mm256_add_ps(x, delta);

            // Update only lanes that still need iterations
            x = _mm256_blendv_ps(x, xnew, need);
        }

        // sqrt(S) = 1/x
        __m256 sqrtS = _mm256_div_ps(one, x);
        _mm256_storeu_ps(output + i, sqrtS);
    }
    
    if ( i < N ){
       sqrtSerial(N-i, startGuess, values + i, output + i); 
    }
    
}

