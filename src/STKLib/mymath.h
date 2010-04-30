#ifndef _MYMATH_H
#define _MYMATH_H

// NOTE: this module does not solve speed in dspc.cpp for complex sqrt and abs

#include <math.h>
#include <float.h>
#include <assert.h>

// 64 bits integer type
#ifndef long_long
  typedef long long long_long;
#endif

#ifdef WIN32
  #if defined(__BORLANDC__) && ! defined(__MATHF__)
    #define __MATHF__
    inline float logf(float arg) { return log(arg); }
    inline float sqrtf(float arg) { return sqrt(arg); }
    inline float cosf(float arg) { return cos(arg); }
    inline float expf(float arg) { return exp(arg); }
    inline float fabsf(float arg) { return fabs(arg); }
  #endif
#endif

#ifndef SWAP
#define SWAP(a,b) tempr=(a); (a)=(b); (b)=tempr
#endif

#ifndef M_PI
#define M_PI  3.1415926535897932384626433832795
#endif

#ifndef M_EXP_MINUS_PI
#define M_EXP_MINUS_PI   0.043213918263772
#endif

#ifndef M_LOG_2
#define M_LOG_2  0.69314718055994530941723212145818
#endif

#ifndef M_SQRT_2
#define M_SQRT_2  1.4142135623730950488016887242097
#endif

#ifndef M_LOG_0
#define M_LOG_0   -1.0e10
#endif

#ifndef M_LOG_MIN
#define M_LOG_MIN  0.5 * M_LOG_0
#endif

#ifndef M_FLT_EPSILON
#define M_FLT_EPSILON  1.192093e-7
#endif

#ifndef M_FLT_MIN_LOG_DIFF
#define M_FLT_MIN_LOG_DIFF -15.942385065
#endif

#if (defined __GNUC__ && defined __FAST_MATH__)

  #define my_log(v) __builtin_log(v)
  #define my_log10(v) __builtin_log10(v)
  #define my_sqrt(v) __builtin_sqrt(v)
  #define my_cos(v) __builtin_cos(v)
  #define my_sin(v) __builtin_sin(v)
  #define my_exp(v) __builtin_exp(v)
  #define my_fabs(v) __builtin_fabs(v)
  #define my_pow(v1, v2) __builtin_pow(v1, v2)

  #ifdef USE_FMATH

    #include "fmath.h"

    #define my_logf(v)  fmath::log(v)
    #define my_log10f(v)  __builtin_log10f(v)
    #define my_sqrtf(v)  __builtin_sqrtf(v)
    #define my_cosf(v) __builtin_cosf(v)
    #define my_sinf(v) __builtin_sinf(v)
    #define my_expf(v) fmath::exp(v)
    #define my_fabsf(v) __builtin_fabsf(v)
    #define my_powf(v1, v2) __builtin_powf(v1, v2)

  #else

    #define my_logf(v)  __builtin_logf(v)
    #define my_log10f(v)  __builtin_log10f(v)
    #define my_sqrtf(v)  __builtin_sqrtf(v)
    #define my_cosf(v) __builtin_cosf(v)
    #define my_sinf(v) __builtin_sinf(v)
    #define my_expf(v) __builtin_expf(v)
    #define my_fabsf(v) __builtin_fabsf(v)
    #define my_powf(v1, v2) __builtin_powf(v1, v2)

  #endif

#else

  #define my_log(v) log(v)
  #define my_log10(v) log10(v)
  #define my_sqrt(v) sqrt(v)
  #define my_cos(v) cos(v)
  #define my_sin(v) sin(v)
  #define my_exp(v) exp(v)
  #define my_fabs(v) fabs(v)
  #define my_pow(v1, v2) pow(v1, v2)

  #ifdef USE_FMATH

    #include "fmath.h"

    #define my_logf(v) fmath::log(v)
    #define my_log10f(v) log10f(v)
    #define my_sqrtf(v) sqrtf(v)
    #define my_cosf(v) cosf(v)
    #define my_sinf(v) sinf(v)
    #define my_expf(v) fmath::exp(v)
    #define my_fabsf(v) fabsf(v)
    #define my_powf(v1, v2) powf(v1, v2)

  #else

    #define my_logf(v) logf(v)
    #define my_log10f(v) log10f(v)
    #define my_sqrtf(v) sqrtf(v)
    #define my_cosf(v) cosf(v)
    #define my_sinf(v) sinf(v)
    #define my_expf(v) expf(v)
    #define my_fabsf(v) fabsf(v)
    #define my_powf(v1, v2) powf(v1, v2)

  #endif

#endif

static inline float my_logaddf(float x, float y)
{
  float diff;

  if (x < y) 
  {
    diff = x - y;
    x = y;
  }
  else
  {
    diff = y - x;
  }

  if (x < M_LOG_MIN)
    return M_LOG_0;

  if (diff <  M_FLT_MIN_LOG_DIFF) 
    return  x;

  return x + my_logf(1.0 + my_expf(diff));
}

static inline float my_logsubf(float x, float y)
{
  float diff = y - x;
  assert(diff <= 0.0);

  if (x < M_LOG_MIN || diff > -M_FLT_EPSILON) 
    return M_LOG_0;

  if(diff <  M_FLT_MIN_LOG_DIFF)
    return  x;

  return x + my_logf(1.0 - my_expf(diff));
}


#endif
