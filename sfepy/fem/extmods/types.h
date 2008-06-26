#ifndef _TYPES_H_
#define _TYPES_H_

#include <string.h>

#define CONVF64 (double)
#define CONVPF64 (double *)
#define PRINTF64 "%f"
#define SCANF64 "%lf"
#define SCANI16 "%hd"

#ifndef Bool
//typedef enum { False, True } Bool;
#endif

/*!
  @par Revision history:
  21.05.2002, c
*/
typedef enum {
  TYPE_Char,
  TYPE_UChar,
  TYPE_Int16,
  TYPE_UInt16,
  TYPE_Int32,
  TYPE_UInt32,
  TYPE_Float32,
  TYPE_Float64,
} BasicType;

typedef unsigned char uchar;
typedef short int16;
typedef unsigned short uint16;
typedef int int32;
typedef unsigned int uint32;
typedef float float32;
typedef double float64;

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
