/*!
  @par Revision history:
  - 18.09.2006, c
*/
#ifndef _GEOMMECH_H_
#define _GEOMMECH_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"

extern int32 t2i[];
extern int32 t2j[];
extern int32 t4s[];

int32 geme_invert3x3( FMField *mtxI, FMField *mtx );
int32 geme_tensor2vectorS3( FMField *vec, FMField *mtx );
int32 geme_det3x3( float64 *det, FMField *mtx );
int32 geme_trace3x3( float64 *tr, FMField *mtx );
int32 geme_norm3( float64 *out, FMField *mtx );
int32 geme_eig3x3( float64 *out, FMField *mtx );
int32 geme_mulAVSB3( FMField *out, FMField *vs, FMField *in );

int32 geme_mulT2ST2S_T4S_ikjl( FMField *t4, FMField *t21, FMField *t22 );
int32 geme_mulT2ST2S_T4S_iljk( FMField *t4, FMField *t21, FMField *t22 );

int32 geme_elementVolume( float64 *volume, float64 *jacobian, int32 nQP );

int32 geme_buildOpOmega_VS3( float64 *pomega, float64 *pdir,
                             int32 nItem, int32 dim, int32 sym );
int32 geme_projectToDir( float64 *pdef, float64 *pomega,
                         float64 *pstrain, int32 nItem, int32 size );
int32 geme_bfMtx( FMField *out, FMField *bf, FMField *in );

int32 bf_act( FMField *out, FMField *bf, FMField *in );
int32 bf_ract( FMField *out, FMField *bf, FMField *in );
int32 bf_actt( FMField *out, FMField *bf, FMField *in );
int32 bf_actt_c1( FMField *out, FMField *bf, FMField *in );
int32 bf_buildFTF( FMField *ftf, FMField *ftf1 );

END_C_DECLS

#endif /* Header */
