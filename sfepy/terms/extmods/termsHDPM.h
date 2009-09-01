/*!
  @par Revision history:
  - 03.08.2006, c
*/
#ifndef _TERMSHDPM_H_
#define _TERMSHDPM_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"
#include "geometry.h"

int32 d_hdpm_surfdvel( FMField *out, FMField *state, int32 offset,
		       FMField *mtxD, SurfaceGeometry *sg,
		       int32 *fis, int32 nFa, int32 nFP,
		       int32 *faceList, int32 faceList_nRow,
		       int32 *conn, int32 nEl, int32 nEP,
		       int32 *elList, int32 elList_nRow );

END_C_DECLS

#endif /* Header */
