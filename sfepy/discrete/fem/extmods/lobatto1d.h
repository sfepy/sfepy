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

// Lobatto functions.
float64 lobatto_000(float64 x);
float64 lobatto_001(float64 x);
float64 lobatto_002(float64 x);
float64 lobatto_003(float64 x);
float64 lobatto_004(float64 x);
float64 lobatto_005(float64 x);
float64 lobatto_006(float64 x);
float64 lobatto_007(float64 x);
float64 lobatto_008(float64 x);
float64 lobatto_009(float64 x);
float64 lobatto_010(float64 x);

// Derivatives of Lobatto functions.
float64 d_lobatto_000(float64 x);
float64 d_lobatto_001(float64 x);
float64 d_lobatto_002(float64 x);
float64 d_lobatto_003(float64 x);
float64 d_lobatto_004(float64 x);
float64 d_lobatto_005(float64 x);
float64 d_lobatto_006(float64 x);
float64 d_lobatto_007(float64 x);
float64 d_lobatto_008(float64 x);
float64 d_lobatto_009(float64 x);
float64 d_lobatto_010(float64 x);

// Kernel functions.
float64 kernel_002(float64 x);
float64 kernel_003(float64 x);
float64 kernel_004(float64 x);
float64 kernel_005(float64 x);
float64 kernel_006(float64 x);
float64 kernel_007(float64 x);
float64 kernel_008(float64 x);
float64 kernel_009(float64 x);
float64 kernel_010(float64 x);

// Derivatives of kernel functions.
float64 d_kernel_002(float64 x);
float64 d_kernel_003(float64 x);
float64 d_kernel_004(float64 x);
float64 d_kernel_005(float64 x);
float64 d_kernel_006(float64 x);
float64 d_kernel_007(float64 x);
float64 d_kernel_008(float64 x);
float64 d_kernel_009(float64 x);
float64 d_kernel_010(float64 x);

// Legendre functions.
float64 legendre_000(float64 x);
float64 legendre_001(float64 x);
float64 legendre_002(float64 x);
float64 legendre_003(float64 x);
float64 legendre_004(float64 x);
float64 legendre_005(float64 x);
float64 legendre_006(float64 x);
float64 legendre_007(float64 x);
float64 legendre_008(float64 x);
float64 legendre_009(float64 x);
float64 legendre_010(float64 x);

// Derivatives of Legendre functions.
float64 d_legendre_000(float64 x);
float64 d_legendre_001(float64 x);
float64 d_legendre_002(float64 x);
float64 d_legendre_003(float64 x);
float64 d_legendre_004(float64 x);
float64 d_legendre_005(float64 x);
float64 d_legendre_006(float64 x);
float64 d_legendre_007(float64 x);
float64 d_legendre_008(float64 x);
float64 d_legendre_009(float64 x);
float64 d_legendre_010(float64 x);

// End of generated code.

#endif /* !LOBATTO1D_H */
