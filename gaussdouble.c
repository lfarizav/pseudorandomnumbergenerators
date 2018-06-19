# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <math.h>
# include <time.h>
#include <immintrin.h>

# include "gaussdouble.h"
void sincos256_ps(__m256 x, __m256 *s, __m256 *c) {
  __m256 xmm1, xmm2, xmm3 = _mm256_setzero_ps(), sign_bit_sin, y;
  __m256i imm0, imm2, imm4;
  
  sign_bit_sin = x;
  /* take the absolute value */
  x = _mm256_and_ps(x, _mm256_castsi256_ps(_mm256_set1_epi32(~0x80000000)));//_ps_inv_sign_mask
  /* extract the sign bit (upper one) */
  sign_bit_sin = _mm256_and_ps(sign_bit_sin, _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000)));//_ps_sign_mask
  
  /* scale by 4/Pi */
  y = _mm256_mul_ps(x, _mm256_set1_ps(1.27323954473516f));
    

  /* store the integer part of y in imm2 */
  imm2 = _mm256_cvttps_epi32(y);

  /* j=(j+1) & (~1) (see the cephes sources) */
  imm2 = _mm256_add_epi32(imm2, _mm256_set1_epi32(1));//_pi32_1
  imm2 = _mm256_and_si256(imm2, _mm256_set1_epi32(~1));//_pi32_inv1
  y = _mm256_cvtepi32_ps(imm2);

  imm4 = imm2;

  /* get the swap sign flag for the sine */
  imm0 = _mm256_and_si256(imm2, _mm256_set1_epi32(4));//_pi32_4
  imm0 = _mm256_slli_epi32(imm0, 29);

  /* get the polynom selection mask for the sine*/
  imm2 = _mm256_and_si256(imm2, _mm256_set1_epi32(2));//_pi32_2
  imm2 = _mm256_cmpeq_epi32(imm2, _mm256_setzero_si256());

  __m256 swap_sign_bit_sin = _mm256_castsi256_ps(imm0);
  __m256 poly_mask = _mm256_castsi256_ps(imm2);

  /* The magic pass: "Extended precision modular arithmetic" 
     x = ((x - y * DP1) - y * DP2) - y * DP3; */

  xmm1 = _mm256_set1_ps(-0.78515625f);//_ps_minus_cephes_DP1
  xmm2 = _mm256_set1_ps(-2.4187564849853515625e-4f);//_ps_minus_cephes_DP2
  xmm3 = _mm256_set1_ps(-3.77489497744594108e-8f);//_ps_minus_cephes_DP3
  xmm1 = _mm256_mul_ps(y, xmm1);
  xmm2 = _mm256_mul_ps(y, xmm2);
  xmm3 = _mm256_mul_ps(y, xmm3);
  x = _mm256_add_ps(x, xmm1);
  x = _mm256_add_ps(x, xmm2);
  x = _mm256_add_ps(x, xmm3);


  imm4 = _mm256_sub_epi32(imm4, _mm256_set1_epi32(2));//_pi32_2
  imm4 = _mm256_andnot_si256(imm4, _mm256_set1_epi32(4));//_pi32_4
  imm4 = _mm256_slli_epi32(imm4, 29);
  __m256 sign_bit_cos = _mm256_castsi256_ps(imm4);

  sign_bit_sin = _mm256_xor_ps(sign_bit_sin, swap_sign_bit_sin);

  
  /* Evaluate the first polynom  (0 <= x <= Pi/4) */
  __m256 z = _mm256_mul_ps(x,x);
  y = _mm256_set1_ps( 2.443315711809948E-005f);//_ps_coscof_p0

  y = _mm256_mul_ps(y, z);
  y = _mm256_add_ps(y, _mm256_set1_ps(-1.388731625493765E-003f));//_ps_coscof_p1
  y = _mm256_mul_ps(y, z);
  y = _mm256_add_ps(y, _mm256_set1_ps( 4.166664568298827E-002f));//_ps_coscof_p2
  y = _mm256_mul_ps(y, z);
  y = _mm256_mul_ps(y, z);
  __m256 tmp = _mm256_mul_ps(z, _mm256_set1_ps(0.5f));//_ps_0p5
  y = _mm256_sub_ps(y, tmp);
  y = _mm256_add_ps(y, _mm256_set1_ps(1.f));//_ps_1
  
  /* Evaluate the second polynom  (Pi/4 <= x <= 0) */

  __m256 y2 = _mm256_set1_ps(-1.9515295891E-4f);//*(__m128*)_ps_sincof_p0;
  y2 = _mm256_mul_ps(y2, z);
  y2 = _mm256_add_ps(y2, _mm256_set1_ps( 8.3321608736E-3f));//*(__m128*)_ps_sincof_p1);
  y2 = _mm256_mul_ps(y2, z);
  y2 = _mm256_add_ps(y2, _mm256_set1_ps(-1.6666654611E-1f));//_ps_sincof_p2
  y2 = _mm256_mul_ps(y2, z);
  y2 = _mm256_mul_ps(y2, x);
  y2 = _mm256_add_ps(y2, x);

  /* select the correct result from the two polynoms */  
  xmm3 = poly_mask;
  __m256 ysin2 = _mm256_and_ps(xmm3, y2);
  __m256 ysin1 = _mm256_andnot_ps(xmm3, y);
  y2 = _mm256_sub_ps(y2,ysin2);
  y = _mm256_sub_ps(y, ysin1);

  xmm1 = _mm256_add_ps(ysin1,ysin2);
  xmm2 = _mm256_add_ps(y,y2);
 
  /* update the sign */
  *s = _mm256_xor_ps(xmm1, sign_bit_sin);
  *c = _mm256_xor_ps(xmm2, sign_bit_cos);
}


void sincos_ps(__m128 x, __m128 *s, __m128 *c) {
  __m128 xmm1, xmm2, xmm3 = _mm_setzero_ps(), sign_bit_sin, y;
  __m128i emm0, emm2, emm4;

  sign_bit_sin = x;
  /* take the absolute value */
  x = _mm_and_ps(x, _mm_castsi128_ps(_mm_set1_epi32(~0x80000000)));//_ps_inv_sign_mask
  /* extract the sign bit (upper one) */
  sign_bit_sin = _mm_and_ps(sign_bit_sin, _mm_castsi128_ps(_mm_set1_epi32(0x80000000)));//_ps_sign_mask
  
  /* scale by 4/Pi */
  y = _mm_mul_ps(x, _mm_set1_ps(1.27323954473516f));
    

  /* store the integer part of y in emm2 */
  emm2 = _mm_cvttps_epi32(y);

  /* j=(j+1) & (~1) (see the cephes sources) */
  emm2 = _mm_add_epi32(emm2, _mm_set1_epi32(1));//_pi32_1
  emm2 = _mm_and_si128(emm2, _mm_set1_epi32(~1));//_pi32_inv1
  y = _mm_cvtepi32_ps(emm2);

  emm4 = emm2;

  /* get the swap sign flag for the sine */
  emm0 = _mm_and_si128(emm2, _mm_set1_epi32(4));//_pi32_4
  emm0 = _mm_slli_epi32(emm0, 29);
  __m128 swap_sign_bit_sin = _mm_castsi128_ps(emm0);

  /* get the polynom selection mask for the sine*/
  emm2 = _mm_and_si128(emm2, _mm_set1_epi32(2));//_pi32_2
  emm2 = _mm_cmpeq_epi32(emm2, _mm_setzero_si128());
  __m128 poly_mask = _mm_castsi128_ps(emm2);

  /* The magic pass: "Extended precision modular arithmetic" 
     x = ((x - y * DP1) - y * DP2) - y * DP3; */

  xmm1 = _mm_set1_ps(-0.78515625f);//_ps_minus_cephes_DP1
  xmm2 = _mm_set1_ps(-2.4187564849853515625e-4f);//_ps_minus_cephes_DP2
  xmm3 = _mm_set1_ps(-3.77489497744594108e-8f);//_ps_minus_cephes_DP3
  xmm1 = _mm_mul_ps(y, xmm1);
  xmm2 = _mm_mul_ps(y, xmm2);
  xmm3 = _mm_mul_ps(y, xmm3);
  x = _mm_add_ps(x, xmm1);
  x = _mm_add_ps(x, xmm2);
  x = _mm_add_ps(x, xmm3);


  emm4 = _mm_sub_epi32(emm4, _mm_set1_epi32(2));//_pi32_2
  emm4 = _mm_andnot_si128(emm4, _mm_set1_epi32(4));//_pi32_4
  emm4 = _mm_slli_epi32(emm4, 29);
  __m128 sign_bit_cos = _mm_castsi128_ps(emm4);

  sign_bit_sin = _mm_xor_ps(sign_bit_sin, swap_sign_bit_sin);

  
  /* Evaluate the first polynom  (0 <= x <= Pi/4) */
  __m128 z = _mm_mul_ps(x,x);
  y = _mm_set1_ps( 2.443315711809948E-005f);//_ps_coscof_p0

  y = _mm_mul_ps(y, z);
  y = _mm_add_ps(y, _mm_set1_ps(-1.388731625493765E-003f));//_ps_coscof_p1
  y = _mm_mul_ps(y, z);
  y = _mm_add_ps(y, _mm_set1_ps( 4.166664568298827E-002f));//_ps_coscof_p2
  y = _mm_mul_ps(y, z);
  y = _mm_mul_ps(y, z);
  __m128 tmp = _mm_mul_ps(z, _mm_set1_ps(0.5f));//_ps_0p5
  y = _mm_sub_ps(y, tmp);
  y = _mm_add_ps(y, _mm_set1_ps(1.f));//_ps_1
  
  /* Evaluate the second polynom  (Pi/4 <= x <= 0) */

  __m128 y2 = _mm_set1_ps(-1.9515295891E-4f);//*(__m128*)_ps_sincof_p0;
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_add_ps(y2, _mm_set1_ps( 8.3321608736E-3f));//*(__m128*)_ps_sincof_p1);
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_add_ps(y2, _mm_set1_ps(-1.6666654611E-1f));//_ps_sincof_p2
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_mul_ps(y2, x);
  y2 = _mm_add_ps(y2, x);

  /* select the correct result from the two polynoms */  
  xmm3 = poly_mask;
  __m128 ysin2 = _mm_and_ps(xmm3, y2);
  __m128 ysin1 = _mm_andnot_ps(xmm3, y);
  y2 = _mm_sub_ps(y2,ysin2);
  y = _mm_sub_ps(y, ysin1);

  xmm1 = _mm_add_ps(ysin1,ysin2);
  xmm2 = _mm_add_ps(y,y2);
 
  /* update the sign */
  *s = _mm_xor_ps(xmm1, sign_bit_sin);
  *c = _mm_xor_ps(xmm2, sign_bit_cos);
}


__m128 log_ps(__m128 x) {

  __m128i emm0 __attribute__((aligned(16)));
  __m128 one __attribute__((aligned(16)))=_mm_set1_ps(1.f);
  __m128 invalid_mask __attribute__((aligned(16))) = _mm_cmple_ps(x, _mm_setzero_ps());

  x = _mm_max_ps(x, _mm_castsi128_ps(_mm_set1_epi32(0x00800000)));  // cut off denormalized stuff
  emm0 = _mm_srli_epi32(_mm_castps_si128(x), 23);

  // keep only the fractional part 
  x = _mm_and_ps(x,_mm_castsi128_ps(_mm_set1_epi32(~0x7f800000)));
  //printf("ln inside x is %e,%e,%e,%e\n",x[0],x[1],x[2],x[3]);
  x = _mm_or_ps(x, _mm_set1_ps(0.5f));
  //printf("ln inside x is %e,%e,%e,%e\n",x[0],x[1],x[2],x[3]);


  // now e=mm0:mm1 contain the really base-2 exponent 
  emm0 = _mm_sub_epi32(emm0, _mm_set1_epi32(0x7f));
  __m128 e = _mm_cvtepi32_ps(emm0); 


  e = _mm_add_ps(e, one);


  __m128 mask = _mm_cmplt_ps(x, _mm_set1_ps(0.707106781186547524f));

  __m128 tmp = _mm_and_ps(x, mask);
  x = _mm_sub_ps(x, one);
  e = _mm_sub_ps(e, _mm_and_ps(one, mask));

  x = _mm_add_ps(x, tmp);


  __m128 z = _mm_mul_ps(x,x);

  __m128 y = _mm_set1_ps(7.0376836292E-2f);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, _mm_set1_ps(- 1.1514610310E-1f));
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, _mm_set1_ps(1.1676998740E-1f));
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, _mm_set1_ps(- 1.2420140846E-1f));
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, _mm_set1_ps(1.4249322787E-1f));
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, _mm_set1_ps(- 1.6668057665E-1f));
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, _mm_set1_ps(2.0000714765E-1f));
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, _mm_set1_ps(- 2.4999993993E-1f));
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, _mm_set1_ps(3.3333331174E-1f));
  y = _mm_mul_ps(y, x);

  y = _mm_mul_ps(y, z);
  

  tmp = _mm_mul_ps(e, _mm_set1_ps(-2.12194440e-4f));
  y = _mm_add_ps(y, tmp);


  tmp = _mm_mul_ps(z, _mm_set1_ps(0.5f));
  y = _mm_sub_ps(y, tmp);

  tmp = _mm_mul_ps(e, _mm_set1_ps(0.693359375f));
  x = _mm_add_ps(x, y);
  x = _mm_add_ps(x, tmp);
  x = _mm_or_ps(x, invalid_mask); // negative arg will be NAN

  return x;
}
__m256 log256_ps(__m256 x) {

  __m256i imm0 __attribute__((aligned(32)));
  __m256 one __attribute__((aligned(32)))=_mm256_set1_ps(1.f);
  __m256 invalid_mask __attribute__((aligned(32))) = _mm256_cmp_ps(x, _mm256_setzero_ps(),_CMP_LE_OS);

  x = _mm256_max_ps(x, _mm256_castsi256_ps(_mm256_set1_epi32(0x00800000)));  // cut off denormalized stuff
  imm0 = _mm256_srli_epi32(_mm256_castps_si256(x), 23);

  // keep only the fractional part 
  x = _mm256_and_ps(x,_mm256_castsi256_ps(_mm256_set1_epi32(~0x7f800000)));
  //printf("ln inside x is %e,%e,%e,%e\n",x[0],x[1],x[2],x[3]);
  x = _mm256_or_ps(x, _mm256_set1_ps(0.5f));
  //printf("ln inside x is %e,%e,%e,%e\n",x[0],x[1],x[2],x[3]);


  // now e=mm0:mm1 contain the really base-2 exponent 
  imm0 = _mm256_sub_epi32(imm0, _mm256_set1_epi32(0x7f));
  __m256 e = _mm256_cvtepi32_ps(imm0); 


  e = _mm256_add_ps(e, one);


  __m256 mask = _mm256_cmp_ps(x, _mm256_set1_ps(0.707106781186547524f),_CMP_LT_OS);

  __m256 tmp = _mm256_and_ps(x, mask);
  x = _mm256_sub_ps(x, one);
  e = _mm256_sub_ps(e, _mm256_and_ps(one, mask));

  x = _mm256_add_ps(x, tmp);


  __m256 z = _mm256_mul_ps(x,x);

  __m256 y = _mm256_set1_ps(7.0376836292E-2f);
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, _mm256_set1_ps(- 1.1514610310E-1f));
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, _mm256_set1_ps(1.1676998740E-1f));
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, _mm256_set1_ps(- 1.2420140846E-1f));
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, _mm256_set1_ps(1.4249322787E-1f));
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, _mm256_set1_ps(- 1.6668057665E-1f));
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, _mm256_set1_ps(2.0000714765E-1f));
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, _mm256_set1_ps(- 2.4999993993E-1f));
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, _mm256_set1_ps(3.3333331174E-1f));
  y = _mm256_mul_ps(y, x);

  y = _mm256_mul_ps(y, z);
  

  tmp = _mm256_mul_ps(e, _mm256_set1_ps(-2.12194440e-4f));
  y = _mm256_add_ps(y, tmp);


  tmp = _mm256_mul_ps(z, _mm256_set1_ps(0.5f));
  y = _mm256_sub_ps(y, tmp);

  tmp = _mm256_mul_ps(e, _mm256_set1_ps(0.693359375f));
  x = _mm256_add_ps(x, y);
  x = _mm256_add_ps(x, tmp);
  x = _mm256_or_ps(x, invalid_mask); // negative arg will be NAN

  return x;
}
__m128 exp_ps(__m128 x) {
  __m128 tmp = _mm_setzero_ps(), fx;

  __m128i emm0;

  __m128 one = _mm_set1_ps(1.f);

  x = _mm_min_ps(x, _mm_set1_ps(88.3762626647949f));
  x = _mm_max_ps(x, _mm_set1_ps(-88.3762626647949f));

  /* express exp(x) as exp(g + n*log(2)) */
  fx = _mm_mul_ps(x, _mm_set1_ps(1.44269504088896341f));
  fx = _mm_add_ps(fx, _mm_set1_ps(0.5f));

  /* how to perform a floorf with SSE: just below */

  emm0 = _mm_cvttps_epi32(fx);
  tmp  = _mm_cvtepi32_ps(emm0);

  /* if greater, substract 1 */
  __m128 mask = _mm_cmpgt_ps(tmp, fx);    
  mask = _mm_and_ps(mask, one);
  fx = _mm_sub_ps(tmp, mask);

  tmp = _mm_mul_ps(fx, _mm_set1_ps(0.693359375f));
  __m128 z = _mm_mul_ps(fx, _mm_set1_ps(-2.12194440e-4f));
  x = _mm_sub_ps(x, tmp);
  x = _mm_sub_ps(x, z);

  z = _mm_mul_ps(x,x);
  
  __m128 y = _mm_set1_ps(1.9875691500E-4f);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, _mm_set1_ps(1.3981999507E-3f));
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, _mm_set1_ps(8.3334519073E-3f));
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, _mm_set1_ps(4.1665795894E-2f));
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, _mm_set1_ps(1.6666665459E-1f));
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, _mm_set1_ps(5.0000001201E-1f));
  y = _mm_mul_ps(y, z);
  y = _mm_add_ps(y, x);
  y = _mm_add_ps(y, one);

  /* build 2^n */
  emm0 = _mm_cvttps_epi32(fx);
  emm0 = _mm_add_epi32(emm0, _mm_set1_epi32(0x7f));
  emm0 = _mm_slli_epi32(emm0, 23);
  __m128 pow2n = _mm_castsi128_ps(emm0);
  y = _mm_mul_ps(y, pow2n);
  return y;
}

__m256 exp256_ps(__m256 x) {
  __m256 tmp = _mm256_setzero_ps(), fx;

  __m256i imm0;

  __m256 one = _mm256_set1_ps(1.f);

  x = _mm256_min_ps(x, _mm256_set1_ps(88.3762626647949f));
  x = _mm256_max_ps(x, _mm256_set1_ps(-88.3762626647949f));

  /* express exp(x) as exp(g + n*log(2)) */
  fx = _mm256_mul_ps(x, _mm256_set1_ps(1.44269504088896341f));
  fx = _mm256_add_ps(fx, _mm256_set1_ps(0.5f));

  /* how to perform a floorf with SSE: just below */

  /*emm0 = _mm_cvttps_epi32(fx);
  tmp  = _mm_cvtepi32_ps(emm0);*/
  tmp = _mm256_floor_ps(fx);

  /* if greater, substract 1 */
  __m256 mask = _mm256_cmp_ps(tmp, fx, _CMP_GT_OS);    
  mask = _mm256_and_ps(mask, one);
  fx = _mm256_sub_ps(tmp, mask);

  tmp = _mm256_mul_ps(fx, _mm256_set1_ps(0.693359375f));
  __m256 z = _mm256_mul_ps(fx, _mm256_set1_ps(-2.12194440e-4f));
  x = _mm256_sub_ps(x, tmp);
  x = _mm256_sub_ps(x, z);

  z = _mm256_mul_ps(x,x);
  
  __m256 y = _mm256_set1_ps(1.9875691500E-4f);
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, _mm256_set1_ps(1.3981999507E-3f));
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, _mm256_set1_ps(8.3334519073E-3f));
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, _mm256_set1_ps(4.1665795894E-2f));
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, _mm256_set1_ps(1.6666665459E-1f));
  y = _mm256_mul_ps(y, x);
  y = _mm256_add_ps(y, _mm256_set1_ps(5.0000001201E-1f));
  y = _mm256_mul_ps(y, z);
  y = _mm256_add_ps(y, x);
  y = _mm256_add_ps(y, one);

  /* build 2^n */
  imm0 = _mm256_cvttps_epi32(fx);
  imm0 = _mm256_add_epi32(imm0, _mm256_set1_epi32(0x7f));
  imm0 = _mm256_slli_epi32(imm0, 23);
  __m256 pow2n = _mm256_castsi256_ps(imm0);
  y = _mm256_mul_ps(y, pow2n);
  return y;
}

static double wn[128],fn[128];
static uint32_t iz,jz,jsr=123456789,kn[128];
static int32_t hz;
static unsigned int seed, iy, ir[98];

#define SHR3 (jz=jsr, jsr^=(jsr<<13),jsr^=(jsr>>17),jsr^=(jsr<<5),jz+jsr)
#define UNI (0.5+(signed) SHR3 * 0.2328306e-9)
#define NOR (hz=SHR3,iz=(hz&127),(abs(hz)<kn[iz])? hz*wn[iz] : nfix())


#define a 1664525lu
#define mod 4294967296.0                /* is 2**32 */

#if 1
void randominit(unsigned seed_init)
{
  int i;
  // this need to be integrated with the existing rng, like taus: navid
  printf("Initializing random number generator, seed %d\n",seed_init);

  if (seed_init == 0) {
    srand((unsigned)time(NULL));

    seed = (unsigned long) rand();
  } else {
    seed = seed_init;
  }

  if (seed % 2 == 0) seed += 1; /* seed and mod are relative prime */

  for (i=1; i<=97; i++) {
    seed = a*seed;                 /* mod 2**32  */
    ir[i]= seed;                   /* initialize the shuffle table    */
  }

 iy=1;
}
#endif

double uniformrandom(void)
{
#define a 1664525lu
#define mod 4294967296.0                /* is 2**32 */

  int j;

  j=1 + 97.0*iy/mod;
  iy=ir[j];
  seed = a*seed;                          /* mod 2**32 */
  ir[j] = seed;
  return( (double) iy/mod );
}

double gaussdouble(double mean, double variance)
{
  static long iset=0;
  static double gset;
  double fac,r,v1,v2;
  static double max=-1000000;
  static double min=1000000;

  if (iset == 0) {
    do {
      v1 = 2.0*uniformrandom()-1.0;
      v2 = 2.0*uniformrandom()-1.0;
      r = v1*v1+v2*v2;
    }  while (r >= 1.0);
    fac = sqrt(-2.0*log(r)/r);
    gset= v1*fac;
    iset=1;
    return(sqrt(variance)*v2*fac + mean);
  } else {
    iset=0;
    if (max<sqrt(variance)*gset + mean)
	max=sqrt(variance)*gset + mean;
    if (min>sqrt(variance)*gset + mean)
	min=sqrt(variance)*gset + mean;

    return(sqrt(variance)*gset + mean);
  }
}

double nfix(void)
{
  const double r = 3.442620; 
  static double x, y;
  for (;;)
  {
      x=hz *  wn[iz];
      if (iz==0)
      {   
        do
        {
          x = - 0.2904764 * log (UNI);
          y = - log (UNI);
	} 
        while (y+y < x*x);
        return (hz>0)? r+x : -r-x;
      }
      if (fn[iz]+UNI*(fn[iz-1]-fn[iz])<exp(-0.5*x*x)){
        return x;
      }
      hz = SHR3;
      iz = hz&127;
      if (abs(hz) < kn[iz]){
        return ((hz)*wn[iz]);
      }
  }
}
/*!\Procedure to create tables for normal distribution kn,wn and fn. */
void table_nor(unsigned long seed)
{
  jsr=seed;
  double dn = 3.442619855899;
  int i;
  const double m1 = 2147483648.0;
  double q;
  double tn = 3.442619855899;
  const double vn = 9.91256303526217E-03;

  q = vn/exp(-0.5*dn*dn);

  kn[0] = ((dn/q)*m1);
  kn[1] = 0;

  wn[0] =  ( q / m1 );
  wn[127] = ( dn / m1 );

  fn[0] = 1.0;
  fn[127] = ( exp ( - 0.5 * dn * dn ) );
  for ( i = 126; 1 <= i; i-- )
  {
    dn = sqrt (-2.0 * log ( vn/dn + exp(-0.5*dn*dn)));
    kn[i+1] = ((dn / tn)*m1);
    tn = dn;
    fn[i] = (exp (-0.5*dn*dn));
    wn[i] = (dn / m1);
  }

  return;
}
double ziggurat(double mean, double variance)
{
  return NOR;
}

//This initialization depends on the seed for nor_table function in oaisim_functions.c file.
static uint32_t jsr4[4] __attribute__((aligned(16))) = {123456789,112548569,985584512,452236879};
static uint32_t jsr8[8] __attribute__((aligned(32))) = {123456789,112548569,985584512,452236879,536457891,654354657,765645365,254676590};
static uint32_t iz4[4] __attribute__((aligned(16)));
static uint32_t iz8[8] __attribute__((aligned(32)));
static int32_t hz4[4] __attribute__((aligned(16)));
static int32_t hz8[8] __attribute__((aligned(32)));
static __m128i jsr_128 __attribute__((aligned(16)));
static __m256i jsr_256 __attribute__((aligned(32)));
static __m128i jz_128 __attribute__((aligned(16)));
static __m256i jz_256 __attribute__((aligned(32)));
static __m128i hz_128 __attribute__((aligned(16)));
static __m256i hz_256 __attribute__((aligned(32)));
static __m128i hz1_128 __attribute__((aligned(16)));
static __m128i hz2_128 __attribute__((aligned(16)));
static __m128i abs_hz_128 __attribute__((aligned(16)));
static __m256i abs_hz_256 __attribute__((aligned(32)));
static __m128i abs_hz1_128 __attribute__((aligned(16)));
static __m128i abs_hz2_128 __attribute__((aligned(16)));
static __m128i iz_128 __attribute__((aligned(16)));
static __m256i iz_256 __attribute__((aligned(32)));
static __m128i iz1_128 __attribute__((aligned(16)));
static __m128i iz2_128 __attribute__((aligned(16)));
static __m128i cmplt_option0_128 __attribute__((aligned(16)));
static __m256i cmplt_option0_256 __attribute__((aligned(32)));
static int count99=0;
static int count0=0;
static int option=0;
static int flag=0;
static int nfix_first_run=0;
static __m128 x __attribute__((aligned(16)));
static __m256 x256 __attribute__((aligned(32)));

#define SHR3_SSE (jsr_128=_mm_loadu_si128((__m128i *)jsr4),jz_128=jsr_128, jsr_128=_mm_xor_si128(_mm_slli_epi32(jsr_128,13),jsr_128),jsr_128=_mm_xor_si128(_mm_srli_epi32(jsr_128,17),jsr_128),jsr_128=_mm_xor_si128(_mm_slli_epi32(jsr_128,5),jsr_128),_mm_storeu_si128((__m128i *)jsr4,jsr_128),_mm_add_epi32(jz_128,jsr_128))

#define UNI_SSE (_mm_add_ps(_mm_mul_ps(_mm_set1_ps(0.2328306e-9),_mm_cvtepi32_ps(SHR3_SSE)),_mm_set1_ps(0.5)))

#define NOR_SSE (hz_128=SHR3_SSE,_mm_storeu_si128((__m128i *)hz4,hz_128),iz_128=_mm_and_si128(hz_128,_mm_set1_epi32(127)),_mm_storeu_si128((__m128i *)iz4,iz_128),abs_hz_128=_mm_and_si128(hz_128, _mm_set1_epi32(~0x80000000)),cmplt_option0_128 = _mm_cmplt_epi32(abs_hz_128,_mm_setr_epi32(kn[iz4[0]],kn[iz4[1]],kn[iz4[2]],kn[iz4[3]])),count99=(count99>99)?0:count99+4,nfix_first_run=(count99>99)?0:nfix_first_run,(_mm_testc_si128(cmplt_option0_128,_mm_setr_epi32(0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF)))?_mm_mul_ps(_mm_cvtepi32_ps(hz_128),_mm_setr_ps(wn[iz4[0]],wn[iz4[1]],wn[iz4[2]],wn[iz4[3]])):nfix_SSE(iz_128))

/*#define SHR3 (jz=jsr, jsr^=(jsr<<13),jsr^=(jsr>>17),jsr^=(jsr<<5),jz+jsr)
#define UNI (0.5+(signed) SHR3 * 0.2328306e-9)
#define NOR (hz=SHR3,iz=(hz&127),(abs(hz)<kn[iz])? hz*wn[iz] : nfix())*/


#define SHR3_AVX (jsr_256=_mm256_loadu_si256((__m256i *)jsr8),jz_256=jsr_256, jsr_256=_mm256_xor_si256(_mm256_slli_epi32(jsr_256,13),jsr_256),jsr_256=_mm256_xor_si256(_mm256_srli_epi32(jsr_256,17),jsr_256),jsr_256=_mm256_xor_si256(_mm256_slli_epi32(jsr_256,5),jsr_256),_mm256_storeu_si256((__m256i *)jsr8,jsr_256),_mm256_add_epi32(jz_256,jsr_256))

#define UNI_AVX (_mm256_add_ps(_mm256_mul_ps(_mm256_set1_ps(0.2328306e-9),_mm256_cvtepi32_ps(SHR3_AVX)),_mm256_set1_ps(0.5)))

#define NOR_AVX (hz_256=SHR3_AVX,_mm256_storeu_si256((__m256i *)hz8,hz_256),iz_256=_mm256_and_si256(hz_256,_mm256_set1_epi32(127)),_mm256_storeu_si256((__m256i *)iz8,iz_256),abs_hz_256=_mm256_and_si256(hz_256, _mm256_set1_epi32(~0x80000000)),cmplt_option0_256 = _mm256_cmpgt_epi32(_mm256_setr_epi32(kn[iz8[0]],kn[iz8[1]],kn[iz8[2]],kn[iz8[3]],kn[iz8[4]],kn[iz8[5]],kn[iz8[6]],kn[iz8[7]]),abs_hz_256),count99=(count99>99)?0:count99+8,nfix_first_run=(count99>99)?0:nfix_first_run,(_mm256_testc_si256(cmplt_option0_256,_mm256_setr_epi32(0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF)))?_mm256_mul_ps(_mm256_cvtepi32_ps(hz_256),_mm256_setr_ps(wn[iz8[0]],wn[iz8[1]],wn[iz8[2]],wn[iz8[3]],wn[iz8[4]],wn[iz8[5]],wn[iz8[6]],wn[iz8[7]])):nfix_AVX(iz_256))


__m256 nfix_AVX(__m256i iz)
{
  __m256 y __attribute__((aligned(32)));
  __m256i cmplt_option1_256 __attribute__((aligned(32)));
  __m256i cmplt_option2_256 __attribute__((aligned(32)));
  int32_t cmplt_option0[8] __attribute__((aligned(32)));
  int32_t cmplt_option1[8] __attribute__((aligned(32)));
  int32_t cmplt_option2[8] __attribute__((aligned(32)));
  float output[24] __attribute__((aligned(32)));
  float x8_option0[8] __attribute__((aligned(32)));
  float x8[8] __attribute__((aligned(32)));
  int i;
  static float r = 3.442620; 
  uint32_t iz8_i[8] __attribute__((aligned(16))) ;

    //x=hz *  wn[iz];
    _mm256_storeu_ps(x8_option0,_mm256_mul_ps(_mm256_cvtepi32_ps(hz_256),_mm256_setr_ps(wn[iz8[0]],wn[iz8[1]],wn[iz8[2]],wn[iz8[3]],wn[iz8[4]],wn[iz8[5]],wn[iz8[6]],wn[iz8[7]])));
    count0=0;
    for (i=0;i<8;i++)
    {
	    if ((int)cmplt_option0_256[i]==0xFFFFFFFF)
	    {
		output[count0]=hz8[i]*wn[iz8[i]];
		count0++;
	    }  
    }
    if ((iz8[0]==0||iz8[1]==0||iz8[2]==0||iz8[3]==0||iz8[4]==0||iz8[5]==0||iz8[6]==0||iz8[7]==0)&&nfix_first_run==0&&count0>0&&flag==1)
    {
		nfix_first_run=1;
		option=1;
		do
		{
		    //x = - 0.2904764 * log (UNI);
		    x256 = _mm256_mul_ps(_mm256_set1_ps(-0.2904764f), log256_ps(UNI_AVX));
		    _mm256_storeu_ps(x8,x256);
		    //y = - log (UNI);
		    y = _mm256_mul_ps(_mm256_set1_ps(-1.0f), log256_ps(UNI_AVX));
		   //(y+y < x*x)?
		    cmplt_option1_256 = _mm256_cvtps_epi32(_mm256_cmp_ps(_mm256_add_ps(y,y),_mm256_mul_ps(x256,x256),_CMP_LT_OS));
		    _mm256_storeu_si256((__m256i *)cmplt_option1,cmplt_option1_256);
		    for (i=0;i<8;i++)
		    {
			    if (cmplt_option1[i]==0x80000000)
			    {
				output[7]=(hz8[i]>0)? x8[i]+r:-x8[i]-r;
			        break;
			    }  
		    }
		}
		while (cmplt_option1[0]!=0x80000000 && cmplt_option1[1]!=0x80000000 && cmplt_option1[2]!=0x80000000 && cmplt_option1[3]!=0x80000000 && cmplt_option1[4]!=0x80000000 && cmplt_option1[5]!=0x80000000 && cmplt_option1[6]!=0x80000000 && cmplt_option1[7]!=0x80000000);
		//return _mm_setr_ps(output[0],output[1],output[2],output[3]);	
		flag=2;   
    }
    else if (iz8[0]>0&&iz8[1]>0&&iz8[2]>0&&iz8[3]>0&&iz8[4]>0&&iz8[5]>0&&iz8[6]>0&&iz8[7]>0&&nfix_first_run==0&&count0>0&&flag==2)
    {
        nfix_first_run=1;
	option=1;
	cmplt_option2_256 = _mm256_cvtps_epi32(_mm256_cmp_ps(_mm256_add_ps(_mm256_setr_ps(fn[iz8[0]],fn[iz8[1]],fn[iz8[2]],fn[iz8[3]],fn[iz8[4]],fn[iz8[5]],fn[iz8[6]],fn[iz8[7]]),_mm256_mul_ps(UNI_AVX,_mm256_sub_ps(_mm256_setr_ps(fn[iz8[0]-1],fn[iz8[1]-1],fn[iz8[2]-1],fn[iz8[3]-1],fn[iz8[4]-1],fn[iz8[5]-1],fn[iz8[6]-1],fn[iz8[7]-1]),_mm256_setr_ps(fn[iz8[0]],fn[iz8[1]],fn[iz8[2]],fn[iz8[3]],fn[iz8[4]],fn[iz8[5]],fn[iz8[6]],fn[iz8[7]])))),exp256_ps(_mm256_mul_ps(_mm256_mul_ps(x256,x256),_mm256_set1_ps(-0.5f))),_CMP_LT_OS));
	_mm256_storeu_si256((__m256i *)cmplt_option2,cmplt_option2_256);
	for (i=0;i<8;i++)
	{
		if (cmplt_option2[i]==0x80000000)
		{
			output[7]=x8_option0[i];
			break; 
		} 
	}
	//return _mm_setr_ps(output[0],output[1],output[2],output[3]);
	flag=1;
    }
    do      {
	    hz_256=SHR3_AVX;
    	    _mm256_storeu_si256((__m256i *)hz8,hz_256);
    	    iz_256=_mm256_and_si256(hz_256,_mm256_set1_epi32(127));
    	    _mm256_storeu_si256((__m256i *)iz8,iz_256);
	    abs_hz_256=_mm256_and_si256(hz_256, _mm256_set1_epi32(~0x80000000));
    	    _mm256_storeu_si256((__m256i *)cmplt_option0,_mm256_cmpgt_epi32(_mm256_setr_epi32(kn[iz8[0]],kn[iz8[1]],kn[iz8[2]],kn[iz8[3]],kn[iz8[4]],kn[iz8[5]],kn[iz8[6]],kn[iz8[7]]),abs_hz_256));
	    for (i=0;i<7-option;i++)
	    {
		    if ((int)cmplt_option0_256[i]==0xFFFFFFFF)
		    {
			output[count0]=hz8[i]*wn[iz8[i]];
			count0++;
		    }  
	    }  
	    }while( count0 < 8 ); 
	    option=0;
	    return _mm256_setr_ps(output[0],output[1],output[2],output[3],output[4],output[5],output[6],output[7]);
}
__m128 nfix_SSE(__m128i iz)
{
  //printf("nfix_SSE count99 = %d\n",count99);
  __m128 y __attribute__((aligned(16)));
  __m128i cmplt_option1_128 __attribute__((aligned(16)));
  __m128i cmplt_option2_128 __attribute__((aligned(16)));
  int32_t cmplt_option0[4] __attribute__((aligned(16)));
  int32_t cmplt_option1[4] __attribute__((aligned(16)));
  int32_t cmplt_option2[4] __attribute__((aligned(16)));
  float output[12] __attribute__((aligned(16)));
  float x4_option0[4] __attribute__((aligned(16)));
  float x4[4] __attribute__((aligned(16)));
  int i;
  static float r = 3.442620; 
  /*for (i=0;i<4;i++)
  {
    printf("cmplt_option[%d] %x\n",i,(int)cmplt_option0_128[i]);
  }*/
    //x=hz *  wn[iz];
    //_mm_storeu_si128((__m128i *)cmplt_option0,cmplt_option0_128);
    _mm_storeu_ps(x4_option0,_mm_mul_ps(_mm_cvtepi32_ps(hz_128),_mm_setr_ps(wn[iz4[0]],wn[iz4[1]],wn[iz4[2]],wn[iz4[3]])));
    count0=0;
    for (i=0;i<4;i++)
    {
	    if ((int)cmplt_option0_128[i]==0xFFFFFFFF)
	    {
		output[count0]=hz4[i]*wn[iz4[i]];
		count0++;
		//printf("count0 %d\n",count0);
	    }  
    }
    
    if ((iz4[0]==0||iz4[1]==0||iz4[2]==0||iz4[3]==0)&&nfix_first_run==0&&count0>0&&flag==2)
    {
		//printf("option 1\n");
		option=1;
		nfix_first_run=1;
		do
		{
		    //x = - 0.2904764 * log (UNI);
		    x = _mm_mul_ps(_mm_set1_ps(-0.2904764f), log_ps(UNI_SSE));
		    _mm_storeu_ps(x4,x);
		    //y = - log (UNI);
		    y = _mm_mul_ps(_mm_set1_ps(-1.0f), log_ps(UNI_SSE));
		   //(y+y < x*x)?
		    cmplt_option1_128 = _mm_cvtps_epi32(_mm_cmplt_ps(_mm_add_ps(y,y),_mm_mul_ps(x,x)));
		    _mm_storeu_si128((__m128i *)cmplt_option1,cmplt_option1_128);
		    for (i=0;i<4;i++)
		    {
			    if (cmplt_option1[i]==0x80000000)
			    {
				output[3]=(hz4[i]>0)? x4[i]+r:-x4[i]-r;
			        break;
			    }  
		    }
		}
		while (cmplt_option1[0]!=0x80000000 && cmplt_option1[1]!=0x80000000 && cmplt_option1[2]!=0x80000000 && cmplt_option1[3]!=0x80000000);
		//return _mm_setr_ps(output[0],output[1],output[2],output[3]);	 
		flag=1;  
    }
    else if (iz4[0]>0&&iz4[1]>0&&iz4[2]>0&&iz4[3]>0&&nfix_first_run==0&&count0>0&&flag==1)
    {
	//printf("option 2\n");
	option=1;
        nfix_first_run=1;
	cmplt_option2_128 = _mm_cvtps_epi32(_mm_cmplt_ps(_mm_add_ps(_mm_setr_ps(fn[iz4[0]],fn[iz4[1]],fn[iz4[2]],fn[iz4[3]]),_mm_mul_ps(UNI_SSE,_mm_sub_ps(_mm_setr_ps(fn[iz4[0]-1],fn[iz4[1]-1],fn[iz4[2]-1],fn[iz4[3]-1]),_mm_setr_ps(fn[iz4[0]],fn[iz4[1]],fn[iz4[2]],fn[iz4[3]])))),exp_ps(_mm_mul_ps(_mm_mul_ps(x,x),_mm_set1_ps(-0.5f)))));
	_mm_storeu_si128((__m128i *)cmplt_option2,cmplt_option2_128);

	for (i=0;i<4;i++)
	{
		if (cmplt_option2[i]==0x80000000)
		{
			output[3]=x4_option0[i];
			break; 
		} 
	}
	//return _mm_setr_ps(output[0],output[1],output[2],output[3]);
	flag=2;
    }
    do {
/*#define NOR_SSE (hz_128=SHR3_SSE,_mm_storeu_si128((__m128i *)hz4,hz_128),iz_128=_mm_and_si128(hz_128,_mm_set1_epi32(127)),_mm_storeu_si128((__m128i *)iz4,iz_128),abs_hz_128=_mm_and_si128(hz_128, _mm_set1_epi32(~0x80000000)),cmplt_option0_128 = _mm_cmplt_epi32(abs_hz_128,_mm_setr_epi32(kn[iz4[0]],kn[iz4[1]],kn[iz4[2]],kn[iz4[3]])),count99=(count99>99)?0:count99+4,nfix_first_run=(count99>99)?0:1,printf("test FFF = %d\n",(_mm_testc_si128(cmplt_option0_128,_mm_setr_epi32(0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF)))),(_mm_testc_si128(cmplt_option0_128,_mm_setr_epi32(0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF)))?_mm_mul_ps(_mm_cvtepi32_ps(hz_128),_mm_setr_ps(wn[iz4[0]],wn[iz4[1]],wn[iz4[2]],wn[iz4[3]])):nfix_SSE(iz_128))*/

	  //printf("filling remaining\n");
	  hz_128=SHR3_SSE;
    	  _mm_storeu_si128((__m128i *)hz4,hz_128);
    	  iz_128=_mm_and_si128(hz_128,_mm_set1_epi32(127));
    	  _mm_storeu_si128((__m128i *)iz4,iz_128);
	  abs_hz_128=_mm_and_si128(hz_128, _mm_set1_epi32(~0x80000000));
    	  cmplt_option0_128 = _mm_cmplt_epi32(abs_hz_128,_mm_setr_epi32(kn[iz4[0]],kn[iz4[1]],kn[iz4[2]],kn[iz4[3]]));

	  /*for (i=0;i<4;i++)
	  { 
	      printf("cmplt_option filling[%d] %x\n",i,(int)cmplt_option0_128[i]);
	  }*/
	  //printf("\n");
	  for (i=0;i<4-option;i++)
	  {
		    if ((int)cmplt_option0_128[i]==0xFFFFFFFF)
		    {
				output[count0]=hz4[i]*wn[iz4[i]];
				count0++;
				//printf("count00 %d\n",count0);
		    }  
	  }
	  }while( count0 < 4 );  
	  option=0;
	  return _mm_setr_ps(output[0],output[1],output[2],output[3]);
}
__m128 ziggurat_SSE_float(void)
{
  return   NOR_SSE;
}
__m256 ziggurat_AVX_float(void)
{
  return   NOR_AVX;
}

void boxmuller_SSE_float(__m128 *data1, __m128 *data2) {
	__m128 twopi = _mm_set1_ps(2.0f * 3.14159265358979323846f);
	__m128 minustwo = _mm_set1_ps(-2.0f);
	__m128 u1_ps,u2_ps;
	__m128 radius,theta,sintheta,costheta;

	u1_ps = UNI_SSE;
	u2_ps = UNI_SSE;
	radius = _mm_sqrt_ps(_mm_mul_ps(minustwo, log_ps(u1_ps)));
	theta = _mm_mul_ps(twopi, u2_ps);
        sincos_ps(theta, &sintheta, &costheta);
	*data1=_mm_mul_ps(radius, costheta);
	*data2=_mm_mul_ps(radius, sintheta);
}
void boxmuller_AVX_float(__m256 *data1, __m256 *data2) {
	__m256 twopi = _mm256_set1_ps(2.0f * 3.14159265358979323846f);
	__m256 minustwo = _mm256_set1_ps(-2.0f);
	__m256 u1_ps,u2_ps;
	__m256 radius,theta,sintheta,costheta;
	u1_ps = UNI_AVX;
	u2_ps = UNI_AVX;
	radius = _mm256_sqrt_ps(_mm256_mul_ps(minustwo, log256_ps(u1_ps)));
	theta = _mm256_mul_ps(twopi, u2_ps);
        sincos256_ps(theta, &sintheta, &costheta);
	*data1=_mm256_mul_ps(radius, costheta);
	*data2=_mm256_mul_ps(radius, sintheta);
}
