# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include <string.h>
# include <unistd.h>
# include <float.h>
# include <immintrin.h>
# include "display.h"

void store_data (double sample)
{
  int i;
  FILE *fpt;
  fpt = fopen("plot.out","a");
  fprintf(fpt, "%e\n",sample);
  fflush(fpt);
  fclose(fpt);
  Xi_square (sample);
}

void store_data_SSE (__m128 sample)
{
  int i;
  float sample4[4] __attribute__((aligned(16)));
  _mm_storeu_ps(sample4,sample);
  FILE *fpt;
  fpt = fopen("plot4.out","a");
  for (i=0;i<4;i++)
  {
	  fprintf(fpt, "%e\n",sample4[i]);
	  //printf("%e\n",sample4[i]);
  }
  fflush(fpt);
  fclose(fpt);
  Xi_square_SSE (sample4);
  //printf("Xi2_4 %e\n",Xi2_4);
}
void store_data_AVX (__m256 sample)
{
  int i;
  float sample8[8] __attribute__((aligned(32)));
  _mm256_storeu_ps(sample8,sample);
  FILE *fpt;
  fpt = fopen("plot8.out","a");
  for (i=0;i<8;i++)
  {
	  fprintf(fpt, "%e\n",sample8[i]);
	  //printf("%e\n",sample4[i]);
  }
  fflush(fpt);
  fclose(fpt);
  Xi_square_AVX (sample8);
}
void store_data_boxmuller_SSE (__m128 a,__m128 b)
{
  int i;
  float a4[4] __attribute__((aligned(16)));
  float b4[4] __attribute__((aligned(16)));
  _mm_storeu_ps(a4,a);
  _mm_storeu_ps(b4,b);
  FILE *fpt;
  fpt = fopen("plot4.out","a");
  for (i=0;i<4;i++)
  {
	  fprintf(fpt, "%e\n",a4[i]);
	  //printf("%e\n",sample4[i]);
  }
  for (i=0;i<4;i++)
  {
	  fprintf(fpt, "%e\n",b4[i]);
	  //printf("%e\n",sample4[i]);
  }
  fflush(fpt);
  fclose(fpt);
  Xi_square_SSE (a4);
  Xi_square_SSE (b4);
}
void store_data_boxmuller_AVX (__m256 a,__m256 b)
{
  int i;
  float a8[8] __attribute__((aligned(32)));
  float b8[8] __attribute__((aligned(32)));
  _mm256_storeu_ps(a8,a);
  _mm256_storeu_ps(b8,b);
  FILE *fpt;
  fpt = fopen("plot8.out","a");
  for (i=0;i<8;i++)
  {
	  fprintf(fpt, "%e\n",a8[i]);
	  //printf("%e\n",sample4[i]);
  }
  for (i=0;i<8;i++)
  {
	  fprintf(fpt, "%e\n",b8[i]);
	  //printf("%e\n",sample4[i]);
  }
  fflush(fpt);
  fclose(fpt);
  Xi_square_AVX (a8);
  Xi_square_AVX (b8);
}

void display (FILE *gp)
{
  fprintf(gp,"n=100\nmax=3.\nmin=-3.\nwidth=(max-min)/n\nhist(x,width)=width*floor(x/width)+width/2.0\nset boxwidth width*0.9\nset style fill solid 0.5\nplot 'plot.out' u (hist($1,width)):(1.0) smooth freq w boxes lc rgb 'green' notitle");
}

void Xi_square (float a)
{
  	Xi2+=pow(a,2);
}
void Xi_square_SSE (float a[4])
{
  int i;
  for (i=0;i<4;i++)
  {
  	Xi2_4+=pow(a[i],2);
  }
}
void Xi_square_AVX (float a[8])
{
  int i;
  for (i=0;i<8;i++)
  {
  	Xi2_8+=pow(a[i],2);
  }
}
