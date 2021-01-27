#ifndef _LOBATTO_H_
#define _LOBATTO_H_ 1

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

extern const int32 max_order;

typedef float64 (*fun)(float64 x);

extern fun lobatto[];
extern fun d_lobatto[];
extern fun kernel[];
extern fun d_kernel[];
extern fun legendre[];
extern fun d_legendre[];

int32 eval_lobatto1d(FMField *out, FMField *coors, int32 order);

int32 eval_lobatto_tensor_product(FMField *out, FMField *coors,
                                  int32 *nodes,
                                  float64 cmin, float64 cmax,
                                  int32 diff);

#endif /* !LOBATTO_H */
