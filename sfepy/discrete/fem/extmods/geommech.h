/*!
  @par Revision history:
  - 18.09.2006, c
*/
#ifndef _GEOMMECH_H_
#define _GEOMMECH_H_

#include "common.h"
BEGIN_C_DECLS

#include "fmfield.h"

#define sym2dim( sym ) ((int32)(sym/3+1))

extern int32 t2i1D[];
extern int32 t2j1D[];
extern int32 t4s1D[];
extern int32 t2i2D[];
extern int32 t2j2D[];
extern int32 t4s2D[];
extern int32 t2i3D[];
extern int32 t2j3D[];
extern int32 t4s3D[];

int32 geme_invert3x3( FMField *mtxI, FMField *mtx );
int32 geme_invert4x4( FMField *mtxI, FMField *mtx );
int32 geme_tensor2vectorS3( FMField *vec, FMField *mtx );
int32 geme_det3x3( float64 *det, FMField *mtx );
int32 geme_trace3x3( float64 *tr, FMField *mtx );
int32 geme_norm3( float64 *out, FMField *mtx );
int32 geme_eig3x3( float64 *out, FMField *mtx );
int32 geme_mulAVSB3( FMField *out, FMField *vs, FMField *in );

int32 geme_mulT2ST2S_T4S_ikjl( FMField *t4, FMField *t21, FMField *t22 );
int32 geme_mulT2ST2S_T4S_iljk( FMField *t4, FMField *t21, FMField *t22 );
int32 geme_mulT2S_AA( FMField *R, FMField *A );

int32 geme_elementVolume( float64 *volume, float64 *jacobian, int32 nQP );

int32 geme_buildOpOmega_VS3( float64 *pomega, float64 *pdir,
                             int32 nItem, int32 dim, int32 sym );
int32 geme_projectToDir( float64 *pdef, float64 *pomega,
                         float64 *pstrain, int32 nItem, int32 size );

int32 bf_act( FMField *out, FMField *bf, FMField *in );
int32 bf_ract( FMField *out, FMField *bf, FMField *in );
int32 bf_actt( FMField *out, FMField *bf, FMField *in );
int32 bf_actt_c1( FMField *out, FMField *bf, FMField *in );
int32 bf_buildFTF( FMField *ftf, FMField *ftf1 );

int32 geme_invar1( float64 *invar, FMField *mtx );
int32 geme_invar2( float64 *invar, FMField *mtx );

void debug_printConn( int32 *conn, int32 nEP );

int32 ele_extractNodalValuesNBN( FMField *out, FMField *in,
				 int32 *conn );

int32 ele_extractNodalValuesDBD( FMField *out, FMField *in,
				 int32 *conn );

END_C_DECLS

#endif /* Header */
