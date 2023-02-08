/*!
  @par Revision history:
  - 11.10.2005, c
*/
#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include "common.h"
BEGIN_C_DECLS

typedef struct Mapping {
  int32 nEl;
  int32 nQP;
  int32 dim;
  int32 nEP;
  FMField *bf;
  FMField *bfGM; // Volume or SurfaceExtra only.
  FMField *det; // detJMR or detJSR.

  FMField *normal; // Surface only.

  FMField *volume;
  float64 totalVolume;
} Mapping;

END_C_DECLS

#endif /* Header */
