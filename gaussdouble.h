# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <math.h>
# include <time.h>
#include <immintrin.h>

double gaussdouble(double mean, double variance);
__m128 ziggurat_SSE_float(void);
__m256 ziggurat_AVX_float(void);
void boxmuller_SSE_float(__m128 *data1, __m128 *data2);
void boxmuller_AVX_float(__m256 *data1, __m256 *data2);
void randominit(unsigned seed_init);
double uniformrandom(void);
double ziggurat(double mean, double variance);
void table_nor(unsigned long seed);

void sincos_ps(__m128 x, __m128 *s, __m128 *c);
void sincos256_ps(__m256 x, __m256 *s, __m256 *c);
__m128 log_ps(__m128 x) ;
__m256 log256_ps(__m256 x) ;
__m128 exp_ps(__m128 x);
__m256 exp256_ps(__m256 x);

