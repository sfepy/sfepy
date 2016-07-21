#ifndef _TERMS_OP_H_
#define _TERMS_OP_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "refmaps.h"

int32 mulAB_integrate(FMField *out, FMField *A, FMField *B,
                      Mapping *vg, int32 mode);
int32 actBfT(FMField *out, FMField *bf, FMField *A);
int32 sym2nonsym(FMField *out, FMField *A);

END_C_DECLS

#endif /* Header */
