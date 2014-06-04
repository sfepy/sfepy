#ifndef _NURBS_H_
#define _NURBS_H_ 1

#ifdef __cplusplus
#  define BEGIN_C_DECLS         extern "C" {
#  define END_C_DECLS           }
#else
#  define BEGIN_C_DECLS
#  define END_C_DECLS
#endif

#include "types.h"
#include "version.h"

#include "fmfield.h"

int32 eval_bernstein_basis(FMField *funs, FMField *ders,
                           float64 x, uint32 degree);

#endif /* !NURBS_H */
