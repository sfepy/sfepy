/*!
  @par Revision history:
  - 31.03.2008, c
*/
#ifndef _TERMSBIOT_H_
#define _TERMSBIOT_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "refmaps.h"

int32 dw_biot_grad( FMField *out, float64 coef, FMField *pressure_qp,
		    FMField *mtxD, Mapping *svg, Mapping *vvg,
		    int32 isDiff );

int32 dw_biot_div( FMField *out, float64 coef, FMField *strain,
		   FMField *mtxD, Mapping *svg, Mapping *vvg,
                   int32 isDiff );

int32 d_biot_div( FMField *out, float64 coef, FMField *state, FMField *strain,
		  FMField *mtxD, Mapping *vg );

END_C_DECLS

#endif /* Header */
