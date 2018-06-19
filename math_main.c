# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "gaussdouble.h"

int main ()
{
  int i,j;
  printf("start testing optimized math functions ...\n");
for (j=0;j<20;j++)
{

  double sum_log_SSE=0.0;
  double sum_log=0.0;
  double sum_exp_SSE=0.0;
  double sum_exp=0.0;
  double sum_sincos_SSE=0.0;
  double sum_sincos=0.0;

  __m128 a,b,c,d;
  
  for (i=0;i<1000000;i++)
  {
	float a_float[]={i,i+1,i+2,i+3};
	float b_float[]={i,i+1,i+2,i+3};
	float c_float[]={i,i+1,i+2,i+3};
	float d_float[]={i,i+1,i+2,i+3};
	a=_mm_load_ps(a_float);
	b=_mm_load_ps(b_float);
	c=_mm_load_ps(c_float);
	d=_mm_load_ps(d_float);
	if (i==0)
	{
		/*_mm_store_ps(&a_float[0],log_ps(a));	
		_mm_store_ps(&b_float[0],exp_ps(b));
		printf("Log_ps (0,1,2,3)=(%e,%e,%e,%e)\n",a_float[0],a_float[1],a_float[2],a_float[3]);
		printf("Log (0,1,2,3)=(%e,%e,%e,%e)\n",log(0),log(1),log(2),log(3));
		printf("Exp_ps (0,1,2,3)=(%e,%e,%e,%e)\n",b_float[0],b_float[1],b_float[2],b_float[3]);
		printf("Exp (0,1,2,3)=(%e,%e,%e,%e)\n",exp(0),exp(1),exp(2),exp(3));*/

		sincos_ps(b,&c,&d);
		_mm_store_ps(&c_float[0],c);
		_mm_store_ps(&d_float[0],d);
		printf("sin_ps (0,1,2,3)=(%e,%e,%e,%e)\n",c_float[0],c_float[1],c_float[2],c_float[3]);
		printf("cos_ps (0,1,2,3)=(%e,%e,%e,%e)\n",d_float[0],d_float[1],d_float[2],d_float[3]);
		printf("sin (0,1,2,3)=(%e,%e,%e,%e)\n",sin(0),sin(1),sin(2),sin(3));
		printf("cos (0,1,2,3)=(%e,%e,%e,%e)\n",cos(0),cos(1),cos(2),cos(3));
	}
	/*clock_t start_log_SSE=clock();
	log_ps(a);
      	clock_t stop_log_SSE=clock();
  	sum_log_SSE+=(double) (stop_log_SSE-start_log_SSE)/CLOCKS_PER_SEC;

	clock_t start_log=clock();
	log(i);
      	clock_t stop_log=clock();
  	sum_log+=(double) (stop_log-start_log)/CLOCKS_PER_SEC;


	clock_t start_exp_SSE=clock();
	exp_ps(b);
      	clock_t stop_exp_SSE=clock();
  	sum_exp_SSE+=(double) (stop_exp_SSE-start_exp_SSE)/CLOCKS_PER_SEC;

	clock_t start_exp=clock();
	exp(i);
      	clock_t stop_exp=clock();
  	sum_exp+=(double) (stop_exp-start_exp)/CLOCKS_PER_SEC;*/

	clock_t start_sincos_SSE=clock();
	sincos_ps(b,&c,&d);
      	clock_t stop_sincos_SSE=clock();
  	sum_sincos_SSE+=(double) (stop_sincos_SSE-start_sincos_SSE)/CLOCKS_PER_SEC;

	clock_t start_sincos=clock();
	sin(i);
	cos(i);
      	clock_t stop_sincos=clock();
  	sum_sincos+=(double) (stop_sincos-start_sincos)/CLOCKS_PER_SEC;

  }
  /*printf("Average log SSE time %e (%d times)\n",sum_log_SSE/i,i);
  printf("Average log time %e\n\n",sum_log/i);
  printf("Average exp SSE time %e (%d times)\n",sum_exp_SSE/i,i);
  printf("Average exp time %e\n",sum_exp/i);*/
  printf("Average exp SSE time %e (%d times)\n",sum_sincos_SSE/i,i);
  printf("Average exp time %e\n",sum_sincos/i);
}
  printf("...end\n");
}

