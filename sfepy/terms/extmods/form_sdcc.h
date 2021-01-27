/*!
  @par Revision history:
  - 03.08.2006, c
*/
#ifndef _FORMSDCC_H_
#define _FORMSDCC_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"

int32 form_sdcc_strainCauchy_VS( FMField *strain, FMField *dg );
int32 form_sdcc_actOpGT_VS3( FMField *diff, FMField *gc, FMField *vec );
int32 form_sdcc_actOpGT_M3( FMField *diff, FMField *gc, FMField *mtx );
int32 form_sdcc_actOpG_RM3( FMField *diff, FMField *mtx, FMField *gc );
int32 build_nonsym_grad(FMField *out, FMField *gc);

END_C_DECLS

#endif /* Header */
