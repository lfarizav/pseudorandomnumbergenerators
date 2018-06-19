# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <math.h>
# include <time.h>
#include <immintrin.h>

extern float Xi2;
extern float Xi2_4;
extern float Xi2_8;

void store_data (double sample);
void store_data_SSE (__m128 sample);
void store_data_AVX (__m256 sample);
void store_data_boxmuller_SSE (__m128 a,__m128 b);
void store_data_boxmuller_AVX (__m256 a,__m256 b);
void display (FILE *gp);

void Xi_square (float a);
void Xi_square_SSE (float a[4]);
void Xi_square_AVX (float a[8]);
