#ifndef _LOBATTO1D_H_
#define _LOBATTO1D_H_ 1

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

typedef float64 (*fun)(float64 x);

// Start of generated code.
// REPLACE_TEXT
// End of generated code.

#endif /* !LOBATTO1D_H */
