/*!
  @par Revision history:
  - 08.06.2006, c
*/
#ifndef _GEOMTRANS_H_
#define _GEOMTRANS_H_

#include "common.h"
BEGIN_C_DECLS

#define CONST_MachEps   1e-16

int32 gtr_cross_product( float64 obj[3], float64 obj1[3], float64 obj2[3] );
int32 gtr_normalize_v3( float64 obj[3], float64 obj1[3], int32 dim,
                        int32 verbose );
int32 gtr_dot_v3( float64 *p_val, float64 obj1[3], float64 obj2[3], int32 dim );


END_C_DECLS

#endif /* Header */
