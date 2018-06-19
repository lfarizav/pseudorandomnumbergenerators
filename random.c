# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
/*# include <mkl.h>*/
# include "gaussdouble.h"

int main ()
{
  int i,j;
  printf("start...\n");
for (j=0;j<1;j++)
{
  randominit(1);
  table_nor(123456789);
  FILE *fpt;
  fpt = fopen("plot.out","w");
  FILE *gp;
  gp = popen ("gnuplot -persist","w");

  double sum_boxmuller=0.0;
  double sum_ziggurat=0.0;
  double sum_ziggurat_SSE=0.0;
  double sum_boxmuller_SSE=0.0;

  __m128 a,b;
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
  for (i=0;i<(1000000>>2);i++)
  {
	clock_t start_ziggurat_SSE=clock();
	ziggurat_SSE_float();
      	clock_t stop_ziggurat_SSE=clock();

  	//store_data_SSE (ziggurat_SSE_float());
  	sum_ziggurat_SSE+=(double) (stop_ziggurat_SSE-start_ziggurat_SSE)/CLOCKS_PER_SEC;

	/*clock_t start_boxmuller_SSE=clock();
	boxmuller_SSE_float(&a,&b);
	clock_t stop_boxmuller_SSE=clock();

  	store_data_SSE (a,b);
  	sum_boxmuller_SSE+=(double) (stop_boxmuller_SSE-start_boxmuller_SSE)/CLOCKS_PER_SEC;*/
  }
  printf("Average Ziggurat SSE time %e\n",sum_ziggurat_SSE/i);
  //printf("Average Boxmuller SSE time %e\n",sum_boxmuller_SSE/i);
}
  printf("...end\n");
}

