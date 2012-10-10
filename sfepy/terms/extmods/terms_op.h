#ifndef _TERMS_OP_H_
#define _TERMS_OP_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "refmaps.h"

int32 mulATB_integrate(FMField *out, FMField *A, FMField *B,
		       Mapping *vg);

END_C_DECLS

#endif /* Header */
