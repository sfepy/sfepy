/*!
  @par Revision history:
  - 11.10.2005, c
*/
#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include "common.h"
BEGIN_C_DECLS

#include "geommech.h"

/*!
  @par Revision history:
  - 11.10.2005, c
*/
typedef enum GeometryMode {
  GM_Material,
  GM_Spatial,
} GeometryMode;

/*!
  @par Revision history:
  - 11.10.2005, c
  - 12.10.2005
  - 12.10.2006
*/
typedef struct VolumeGeometry {
  GeometryMode mode;
  int32 nEl;
  int32 nQP;
  int32 dim;
  int32 nEP;
  FMField *bf;
  FMField *bfGM;
  FMField *det; // detJMR or detJSR.
  FMField *volume;
  float64 totalVolume;
} VolumeGeometry;

int32 vg_createAlloc( VolumeGeometry **p_obj,
		      int32 nEl, int32 nQP, int32 dim, int32 nEP );
int32 vg_freeDestroy( VolumeGeometry **p_obj );
int32 vg_print( VolumeGeometry *obj, FILE *file, int32 mode );

int32 vg_describe( VolumeGeometry *obj,
		   float64 *coorIn, int32 nNod, int32 dim,
		   int32 *conn, int32 nEl, int32 nEP,
		   FMField *bfGR, FMField *ebfGR, FMField *weight );
int32 vg_integrate( VolumeGeometry *obj, FMField *out, FMField *in,
                    int32 mode );
int32 vg_getElementDiameters( VolumeGeometry *obj, FMField *out,
			      int32 *edges, int32 edges_nRow, int32 edges_nCol,
			      float64 *coorIn, int32 nNod, int32 dim,
			      int32 *conn, int32 nEl, int32 nEP,
			      int32 *elList, int32 elList_nRow,
			      int32 mode );

/*!
  @par Revision history:
  - 21.12.2005, c
  - 12.10.2006
*/
typedef struct SurfaceGeometry {
  GeometryMode mode;
  int32 nFa;
  int32 nQP;
  int32 dim;
  int32 nFP;
  FMField *normal;
  FMField *det; // detJMR.

  FMField *bf;
  FMField *bfBGM;
  FMField *detF;
  FMField *mtxFI;

  FMField *area;
  float64 totalArea;
} SurfaceGeometry;

int32 sg_createAlloc( SurfaceGeometry **p_obj,
		      int32 nFa, int32 nQP, int32 dim, int32 nSP );
int32 sg_freeDestroy( SurfaceGeometry **p_obj );
int32 sg_print( SurfaceGeometry *obj, FILE *file, int32 mode );

int32 sg_describe( SurfaceGeometry *obj,
		   float64 *coorIn, int32 nNod, int32 dim,
		   int32 *fconn, int32 nFa, int32 nSP,
		   FMField *bfGR, FMField *weight );
int32 sg_integrate( SurfaceGeometry *obj, FMField *out, FMField *in,
                    int32 mode );
int32 sg_evaluateBFBGM( SurfaceGeometry *obj, FMField *bfBGR, FMField *ebfBGR,
			float64 *coorIn, int32 nNod, int32 dim,
			int32 *fis, int32 nFa, int32 nFP,
			int32 *conn, int32 nEl, int32 nEP );

END_C_DECLS

#endif /* Header */
