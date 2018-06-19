# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include <immintrin.h>
# include "mkl.h"
# include "gaussdouble.h"
# include "display.h"

float Xi2=0;
float Xi2_4=0;
float Xi2_8=0;

int main ()
{
  int i,j,k;
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

  double sum_boxmuller=0.0;
  double sum_boxmuller_intel=0.0;
  double sum_boxmuller_SSE=0.0;
  double sum_boxmuller_AVX=0.0;

  __m128 a4,b4;
  __m256 a8,b8;
  VSLStreamStatePtr stream;
  vslNewStream( &stream, VSL_BRNG_MT19937, 777 );
  double test[1];
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
	gaussdouble(0, 1);
      	clock_t stop=clock();
	//printf("inter random box-muller\n");
	clock_t start_intel=clock();
	//vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, stream, 2, test, 0.0, 1.0);
        vdRngGaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, 1, test, 0.0, 1.0 );
      	clock_t stop_intel=clock();
	for (k=0;k<1;k++)
	{
  		//printf("%e\t",test[k]);
	}
	//printf("\n");
	clock_t start_SSE=clock();
	boxmuller_SSE_float(&a4,&b4);
      	clock_t stop_SSE=clock();

	clock_t start_AVX=clock();
	boxmuller_AVX_float(&a8,&b8);
      	clock_t stop_AVX=clock();

  	//store_data_SSE (ziggurat_SSE_float());
  	sum_boxmuller+=(double) (stop-start)/CLOCKS_PER_SEC;
  	sum_boxmuller_intel+=(double) (stop_intel-start_intel)/CLOCKS_PER_SEC;
  	sum_boxmuller_SSE+=(double) (stop_SSE-start_SSE)/CLOCKS_PER_SEC;
  	sum_boxmuller_AVX+=(double) (stop_AVX-start_AVX)/CLOCKS_PER_SEC;
        
	store_data (gaussdouble(0, 1));
        store_data_boxmuller_SSE (a4,b4);
        store_data_boxmuller_AVX (a8,b8);
	
	/*clock_t start_boxmuller_SSE=clock();
	boxmuller_SSE_float(&a,&b);
	clock_t stop_boxmuller_SSE=clock();

  	store_data_SSE (a,b);
  	sum_boxmuller_SSE+=(double) (stop_boxmuller_SSE-start_boxmuller_SSE)/CLOCKS_PER_SEC;*/
  }
  printf("Xi2 is %e, Xi2_4 is %e, Xi2_8 is %e\n",Xi2,Xi2_4,Xi2_8);
  printf("Average Boxmuller       time %e\n",sum_boxmuller/i);
  printf("Average Boxmuller intel time %e\n",sum_boxmuller_intel/i);
  printf("Average Boxmuller  SSE  time %e\n",sum_boxmuller_SSE/i);
  printf("Average Boxmuller  AVX  time %e\n",sum_boxmuller_AVX/i);
}
  printf("...end\n");
}

