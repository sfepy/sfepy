#include "geomtrans.h"

#undef __FUNC__
#define __FUNC__ "gtr_cross_product"
/*!
  @par Revision history:
  - 22.03.2002, c
*/
int32 gtr_cross_product( float64 obj[3], float64 obj1[3], float64 obj2[3] )
{
  obj[0] = obj1[1]*obj2[2]-obj1[2]*obj2[1];
  obj[1] = obj1[2]*obj2[0]-obj1[0]*obj2[2];
  obj[2] = obj1[0]*obj2[1]-obj1[1]*obj2[0];

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "gtr_normalize_v3"
/*!
  @par Revision history:
  - 22.03.2002, c
*/
int32 gtr_normalize_v3( float64 obj[3], float64 obj1[3], int32 dim,
                        int32 verbose )
{
  float64 norm;

  if (dim == 3) {
    norm = sqrt( obj1[0] * obj1[0] + obj1[1] * obj1[1] + obj1[2] * obj1[2] );
    if (norm > CONST_MachEps) {
      obj[0] = obj1[0] / norm;
      obj[1] = obj1[1] / norm;
      obj[2] = obj1[2] / norm;
    } else {
      if (verbose) errput( "zero norm!\n" );
      obj[0] = obj[1] = obj[2] = 0.0;
    }

  } else { // dim == 2
    norm = sqrt( obj1[0] * obj1[0] + obj1[1] * obj1[1] );
    if (norm > CONST_MachEps) {
      obj[0] = obj1[0] / norm;
      obj[1] = obj1[1] / norm;
    } else {
      if (verbose) errput( "zero norm!\n" );
      obj[0] = obj[1] = 0.0;
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "gtr_dot_v3"
/*!
  @par Revision history:
  - 08.06.2006, c
*/
int32 gtr_dot_v3( float64 *p_val, float64 obj1[3], float64 obj2[3], int32 dim )
{
  if (dim == 3) {
    *p_val = obj1[0] * obj2[0] + obj1[1] * obj2[1] + obj1[2] * obj2[2];

  } else { // dim == 2
    *p_val = obj1[0] * obj2[0] + obj1[1] * obj2[1];

  }

  return( RET_OK );
}
