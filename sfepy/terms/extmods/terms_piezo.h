/*!
  @par Revision history:
  - 21.09.2008, c
*/
#ifndef _TERMSPIEZO_H_
#define _TERMSPIEZO_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "refmaps.h"

int32 dw_piezo_coupling( FMField *out, FMField *strain, FMField *charge_grad,
			 FMField *mtxG, Mapping *vg,
			 int32 mode );

int32 d_piezo_coupling( FMField *out, FMField *strain, FMField *charge_grad,
			FMField *mtxG, Mapping *vg );

END_C_DECLS

#endif /* Header */

