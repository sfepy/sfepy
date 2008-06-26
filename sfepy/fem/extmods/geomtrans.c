#include "geomtrans.h"

#undef __FUNC__
#define __FUNC__ "gtr_crossProduct"
/*!
  @par Revision history:
  - 22.03.2002, c
*/
int32 gtr_crossProduct( float64 obj[3], float64 obj1[3], float64 obj2[3] )
{
  obj[0] = obj1[1]*obj2[2]-obj1[2]*obj2[1];
  obj[1] = obj1[2]*obj2[0]-obj1[0]*obj2[2];
  obj[2] = obj1[0]*obj2[1]-obj1[1]*obj2[0];

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "gtr_normalizeV3"
/*!
  @par Revision history:
  - 22.03.2002, c
*/
int32 gtr_normalizeV3( float64 obj[3], float64 obj1[3] )
{
  float64 norm;

  norm = sqrt( obj1[0] * obj1[0] + obj1[1] * obj1[1] + obj1[2] * obj1[2] );
  if (norm > CONST_MachEps) {
    obj[0] = obj1[0] / norm;
    obj[1] = obj1[1] / norm;
    obj[2] = obj1[2] / norm;
  } else {
    errput( "zero norm!\n" );
    obj[0] = obj[1] = obj[2] = 0.0;
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "gtr_dotV3"
/*!
  @par Revision history:
  - 08.06.2006, c
*/
int32 gtr_dotV3( float64 *p_val, float64 obj1[3], float64 obj2[3] )
{
  *p_val = obj1[0] * obj2[0] + obj1[1] * obj2[1] + obj1[2] * obj2[2];

  return( RET_OK );
}
