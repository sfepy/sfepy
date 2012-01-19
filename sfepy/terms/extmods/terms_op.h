#ifndef _TERMS_OP_H_
#define _TERMS_OP_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "geometry.h"

int32 mulATB_integrate(FMField *out, FMField *A, FMField *B,
		       VolumeGeometry *vg);

END_C_DECLS

#endif /* Header */
