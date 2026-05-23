#ifndef __MATRIX_LIB_H__
#define __MATRIX_LIB_H__ 

#ifndef _USE_DOUBLE_PRECISION_
#define _USE_DOUBLE_PRECISION_ (0)
#endif

#if _USE_DOUBLE_PRECISION_
typedef double real_T;
#else
typedef float real_T;
#endif

typedef real_T vector3_T[3];
typedef real_T matrix3_T[3][3];
typedef real_T matrix4_T[4][4];

#endif
