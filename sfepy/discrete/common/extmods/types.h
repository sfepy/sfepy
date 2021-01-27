#ifndef _SFEPY_TYPES_H_
#define _SFEPY_TYPES_H_

#include <Python.h>
#include <numpy/npy_common.h>

#define FI32 "%"NPY_INT32_FMT
#define FF64 "%"NPY_FLOAT64_FMT

typedef unsigned char uchar;
typedef npy_int16 int16;
typedef npy_uint16 uint16;
typedef npy_int32 int32;
typedef npy_uint32 uint32;
typedef npy_float32 float32;
typedef npy_float64 float64;

/*!
  @name Utility macros
  Inspired by umfpack.
  @par Revision history:
  - 04.06.2001
  - 03.04.2003
*/
/*@{*/ 
#define Abs(x) ((x) >= 0 ? (x) : -(x))
#define Max(a,b) (((a) > (b)) ? (a) : (b))
#define Min(a,b) (((a) < (b)) ? (a) : (b))
#define Sgn(a) (((a)>0) ? 1 : (((a)<0) ? -1 : 0))
#define StringMatch(s1,s2) (strcmp ((s1), (s2)) == 0)
#define Implies(p,q) (!(p) || (q))
#define FLG_TestStatus( obj, flag ) (((obj)->status) & (flag))
#define FLG_SetStatus( obj, flag ) (((obj)->status) |= flag)
#define FLG_ClearStatus( obj, flag ) (((obj)->status) &= ~flag)
/*@}*/ 


#endif /* Header */
