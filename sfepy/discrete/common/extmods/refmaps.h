/*!
  @par Revision history:
  - 11.10.2005, c
*/
#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include "common.h"
BEGIN_C_DECLS

#include "geommech.h"

typedef enum MappingMode {
  MM_Volume,
  MM_Surface,
  MM_SurfaceExtra,
} MappingMode;

typedef struct Mapping {
  MappingMode mode;
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

int32 map_print( Mapping *obj, FILE *file, int32 mode );

int32 map_describe( Mapping *obj,
                    float64 *coorIn, int32 nNod, int32 dim,
                    int32 *conn, int32 nEl, int32 nEP,
                    FMField *bfGR, FMField *ebfGR, FMField *weight );
int32 _v_describe( Mapping *obj,
                   float64 *coorIn, int32 nNod, int32 dim,
                   int32 *conn, int32 nEl, int32 nEP,
                   FMField *bfGR, FMField *ebfGR, FMField *weight );
int32 _s_describe( Mapping *obj,
                   float64 *coorIn, int32 nNod, int32 dim,
                   int32 *fconn, int32 nFa, int32 nFP,
                   FMField *bfGR, FMField *weight );

int32 map_integrate( Mapping *obj, FMField *out, FMField *in,
                     int32 mode );

int32 map_getElementDiameters( Mapping *obj, FMField *out,
                               int32 *edges, int32 edges_nRow, int32 edges_nCol,
                               float64 *coorIn, int32 nNod, int32 dim,
                               int32 *conn, int32 nEl, int32 nEP,
                               int32 *elList, int32 elList_nRow,
                               int32 mode );
int32 map_evaluateBFBGM( Mapping *obj, FMField *bfBGR, FMField *ebfBGR,
                         float64 *coorIn, int32 nNod, int32 dim,
                         int32 *fis, int32 nFa, int32 nFP,
                         int32 *conn, int32 nEl, int32 nEP );

END_C_DECLS

#endif /* Header */
