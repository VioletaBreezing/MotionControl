#ifndef __LIE_SOPHUS_H__
#define __LIE_SOPHUS_H__

#define _USE_DOUBLE_PRECISION_ (0)

#if _USE_DOUBLE_PRECISION_
typedef double real_T;
#define func_sin sin
#define func_cos cos
#else
typedef float real_T;
#define func_sin sinf
#define func_cos cosf
#endif

typedef struct
{
	real_T x, y, z;
} vector3_T;

typedef struct
{
    vector3_T rho, phi;
} se3_T;


#endif