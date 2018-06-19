# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
#include <immintrin.h>
# include <mkl.h>
# include "gaussdouble.h"
# include "display.h"

float Xi2=0;
float Xi2_4=0;
float Xi2_8=0;

int main ()
{
  int i,j;
  FILE *fpt;
  FILE *fpt4;
  FILE *fpt8;
  fpt = fopen("plot.out","w");
  fclose(fpt);
  fpt4 = fopen("plot4.out","w");
  fclose(fpt4);
  fpt8 = fopen("plot8.out","w");
  fclose(fpt8);
  printf("start...\n");
  printf("The initial Xi2, Xi2_4, Xi2_8 values are: %e,%e,%e\n",Xi2,Xi2_4,Xi2_8);
for (j=0;j<1;j++)
{
  randominit(1);
  table_nor(123456789);
  FILE *fpt;
  fpt = fopen("plot.out","w");
  FILE *gp;
  gp = popen ("gnuplot -persist","w");

  double sum_ziggurat=0.0;
  double sum_ziggurat_SSE=0.0;
  double sum_ziggurat_AVX=0.0;

  __m128 a4,b4;
  __m256 a8,b8;
  float a;
  /*for (i=0;i<1000000;i++)
  {
    //printf("Inside the for with i = %d\n",i);
  clock_t start_boxmuller=clock();
  gaussdouble(0.0,1.0);
  clock_t stop_boxmuller=clock();
  store_data(ziggurat(0.0,1.0));

  clock_t start_ziggurat=clock();
  ziggurat(0.0,1.0);
  clock_t stop_ziggurat=clock();
  //printf("gaussdouble %e, time %e\n",gaussdouble(0.0,1.0),(double) (stop-start)/CLOCKS_PER_SEC);
  sum_boxmuller+=(double) (stop_boxmuller-start_boxmuller)/CLOCKS_PER_SEC;
  sum_ziggurat+=(double) (stop_ziggurat-start_ziggurat)/CLOCKS_PER_SEC;
  }
  printf("Average Box-Muller time %e\n",sum_boxmuller/i);
  printf("Average Ziggurat time %e\n",sum_ziggurat/i);
  display(gp);*/
  for (i=0;i<(1000000);i++)
  {
	clock_t start=clock();
	ziggurat(0, 1);
      	clock_t stop=clock();

	clock_t start_SSE=clock();
	a4=ziggurat_SSE_float();
      	clock_t stop_SSE=clock();

	clock_t start_AVX=clock();
	a8=ziggurat_AVX_float();
      	clock_t stop_AVX=clock();

  	//store_data_SSE (ziggurat_SSE_float());
  	sum_ziggurat+=(double) (stop-start)/CLOCKS_PER_SEC;
  	sum_ziggurat_SSE+=(double) (stop_SSE-start_SSE)/CLOCKS_PER_SEC;
  	sum_ziggurat_AVX+=(double) (stop_AVX-start_AVX)/CLOCKS_PER_SEC;
        
	store_data (ziggurat(0,1));
        store_data_SSE (a4);
	store_data_AVX (a8);
	/*clock_t start_boxmuller_SSE=clock();
	boxmuller_SSE_float(&a,&b);
	clock_t stop_boxmuller_SSE=clock();

  	store_data_SSE (a,b);
  	sum_boxmuller_SSE+=(double) (stop_boxmuller_SSE-start_boxmuller_SSE)/CLOCKS_PER_SEC;*/
  }
  printf("Xi2 is %e, Xi2_4 is %e, Xi2_8 is %e\n",Xi2,Xi2_4,Xi2_8);
  printf("Average ziggurat     time %e\n",sum_ziggurat/i);
  printf("Average ziggurat SSE time %e\n",sum_ziggurat_SSE/i);
  printf("Average ziggurat AVX time %e\n",sum_ziggurat_AVX/i);
}
  printf("...end\n");
}

