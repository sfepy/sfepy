#include "geometry.h"

#undef __FUNC__
#define __FUNC__ "vg_createAlloc"
/*!
  @par Revision history:
  - 11.10.2005, c
  - 12.10.2005
  - 15.12.2005
  - 12.10.2006
*/
int32 vg_createAlloc( VolumeGeometry **p_obj,
		      int32 nEl, int32 nQP, int32 dim, int32 nEP )
{
  VolumeGeometry *obj;

  obj = allocMem( VolumeGeometry, 1 );
  fmf_createAlloc( &(obj->bfGM), nEl, nQP, dim, nEP );
  fmf_createAlloc( &(obj->det), nEl, nQP, 1, 1 );
  fmf_createAlloc( &(obj->volume), nEl, 1, 1, 1 );
  obj->mode = -1;

  obj->nEl = obj->bfGM->nCell;
  obj->nQP = obj->bfGM->nLev;
  obj->dim = obj->bfGM->nRow;
  obj->nEP = obj->bfGM->nCol;


  *p_obj = obj;

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "vg_freeDestroy"
/*!
  @par Revision history:
  - 11.10.2005, c
  - 12.10.2005
*/
int32 vg_freeDestroy( VolumeGeometry **p_obj )
{
  VolumeGeometry *obj = *p_obj;

  if (!obj) return( RET_OK );

  if (obj->bfGM) 
    fmf_freeDestroy( &(obj->bfGM) );
  if (obj->det) 
    fmf_freeDestroy( &(obj->det) );
  if (obj->volume) 
    fmf_freeDestroy( &(obj->volume) );
  freeMem( *p_obj );

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "vg_print"
/*!
  @par Revision history:
  - 11.10.2005, c
  - 12.10.2005
  - 12.10.2006
*/
int32 vg_print( VolumeGeometry *obj, FILE *file, int32 mode )
{
  int32 ii;

  fprintf( file, "VolumeGeometry: mode %d, nEl %d, nQP %d, dim: %d, nEP: %d\n",
	   obj->mode, obj->nEl, obj->nQP, obj->dim, obj->nEP );
  fprintf( file, "totalVolume: %.5f\n", obj->totalVolume );

  for (ii = 0; ii < obj->det->nCell; ii++) {
    FMF_SetCell( obj->bfGM, ii );
    FMF_SetCell( obj->det, ii );
    FMF_SetCell( obj->volume, ii );
    
    fprintf( file, "%d bfGM:\n", ii );
    fmf_print( obj->bfGM, file, mode );
    
    fprintf( file, "%d det:\n", ii );
    fmf_print( obj->det, file, mode );
    
    fprintf( file, "%d volume:\n", ii );
    fmf_print( obj->volume, file, mode );
  }


  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "vg_describe"
/*!
  @par Revision history:
  - 11.10.2005, c
  - 12.10.2005
  - 05.06.2006
  - 06.06.2006
*/
int32 vg_describe( VolumeGeometry *obj,
		   float64 *coorIn, int32 nNod, int32 dim,
		   int32 *conn, int32 nEl, int32 nEP,
		   FMField *bfGR, FMField *weight )
{
  int32 iel, inod, idim, pos, iqp, nQP, ret = RET_OK;
  FMField *mtxMR = 0, *mtxMRI = 0, *coor = 0;
  
  nQP = bfGR->nLev;
  if (!((nEl == obj->nEl) &&
	(dim == obj->dim) &&
	(nQP == obj->nQP) &&
	(nEP == obj->nEP))) {
    output( "nNod: %d, dim: %d, nEl: %d, nEP: %d\n",  nNod, dim, nEl, nEP );
    fmf_print( obj->bfGM, stdout, 1 );
    fmf_print( bfGR, stdout, 1 );
    fmf_print( weight, stdout, 1 );
    errput( "size mismatch!\n" );
    return( RET_Fail );
  }

  fmf_createAlloc( &mtxMR, 1, nQP, dim, dim );
  fmf_createAlloc( &mtxMRI, 1, nQP, dim, dim );
  fmf_createAlloc( &coor, 1, 1, nEP, dim );

  obj->totalVolume = 0.0;
/*   output( "nCell %d\n",  obj->bfGM->nCell ); */
  for (iel = 0; iel < obj->bfGM->nCell; iel++) {
    FMF_SetCell( obj->bfGM, iel );
    FMF_SetCell( obj->det, iel );
    FMF_SetCell( obj->volume, iel );

    for (inod = 0; inod < nEP; inod++) {
      pos = dim*conn[inod];
      for (idim = 0; idim < dim; idim++ ) {
	coor->val[dim*inod+idim] = coorIn[idim+pos];
      }
    }

    // Jacobi matrix from reference to material system.
    fmf_mulATBT_1n( mtxMR, coor, bfGR );
    // Its determinant, preweighted.
    geme_det3x3( obj->det->val, mtxMR );
    for (iqp = 0; iqp < nQP; iqp++) {
      if (obj->det->val[iqp] <= MachEps) {
	errput( "warp violation %e at (iel: %d, iqp: %d)!\n",
		obj->det->val[iqp], iel, iqp );
      }
    }
    fmf_mul( obj->det, weight->val );

    // Element volume.
    geme_elementVolume( obj->volume->val, obj->det->val, nQP );
    obj->totalVolume += obj->volume->val[0];

    // Inverse of Jacobi matrix reference to material system.
    geme_invert3x3( mtxMRI, mtxMR );
    // Base function gradient w.r.t. material system.
    fmf_mulATB_nn( obj->bfGM, mtxMRI, bfGR );

    conn += nEP;
    
/*     output( "cell %d\n", iel ); */
/*     fmf_print( coor, stdout, 0 ); */
/*     fmf_print( obj->det, stdout, 0 ); */

    ERR_CheckGo( ret );
  }
 end_label:
  fmf_freeDestroy( &mtxMR ); 
  fmf_freeDestroy( &mtxMRI ); 
  fmf_freeDestroy( &coor ); 

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "vg_integrate"
/*!
  @par Revision history:
  - 15.12.2005, c
*/
int32 vg_integrate( VolumeGeometry *obj, FMField *out, FMField *in )
{
  int32 iel;

  for (iel = 0; iel < obj->bfGM->nCell; iel++) {
    FMF_SetCell( obj->det, iel );
    FMF_SetCell( in, iel );
    FMF_SetCell( out, iel );
    fmf_sumLevelsMulF( out, in, obj->det->val );
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "vg_integrateChunk"
/*!
  @par Revision history:
  - 01.11.2007, c
*/
int32 vg_integrateChunk( VolumeGeometry *obj, FMField *out, FMField *in,
			 int32 *elList, int32 elList_nRow )
{
  int32 ii, iel;

  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];
    FMF_SetCell( obj->det, iel );
    FMF_SetCell( out, ii );
    FMF_SetCell( in, ii );
    fmf_sumLevelsMulF( out, in, obj->det->val );
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "vg_getElementDiameters"
/*!
  @par Revision history:
  - 09.01.2006, c
  - 10.01.2006
*/
int32 vg_getElementDiameters( VolumeGeometry *obj, FMField *out,
			      int32 *edges, int32 edges_nRow, int32 edges_nCol,
			      float64 *coorIn, int32 nNod, int32 dim,
			      int32 *conn, int32 nEl, int32 nEP,
			      int32 *elList, int32 elList_nRow,
			      int32 mode )
{
  int32 ii, ie, id, iel, nd;
  float64 val0 = 0.0, val1 = 0.0, vv, aux = 0.0, exponent;

  if ((mode < 0) && (mode > 2)) {
    errput( ErrHead "ERR_Switch\n" );
    return( RET_Fail );
  }

  output( "%d %d %d %d %d %d %d\n",
	  edges_nRow, edges_nCol, nNod, dim, nEl, nEP, elList_nRow );

  nd = obj->bfGM->nRow; // Can be <> dim.
  exponent = 1.0 / ((float64) dim);
  for (ii = 0; ii < elList_nRow; ii++) {
    iel = elList[ii];
    FMF_SetCell( out, ii );

    if ((mode == 0) || (mode == 2)) {
      val0 = 0.0;
      for (ie = 0; ie < edges_nRow; ie++) {
	vv = 0.0;
	for (id = 0; id < nd; id++) {
	  aux = coorIn[dim*conn[nEP*iel+edges[2*ie+1]]+id]
	    - coorIn[dim*conn[nEP*iel+edges[2*ie+0]]+id];
/* 	  output( "%d %d %d %d %f %f %f\n", ii, iel, ie, id, aux, */
/* 		  coorIn[dim*conn[nEP*iel+edges[2*ie+1]]+id], */
/* 		  coorIn[dim*conn[nEP*iel+edges[2*ie+0]]+id] ); */
/* 	  sys_pause(); */
	  vv += aux * aux;
	}
	val0 = Max( val0, aux );
	out->val[0] = val0;
      }
    }
    if ((mode == 1) || (mode == 2)) {
      FMF_SetCell( obj->volume, iel );
      val1 = pow( 0.16 * obj->volume->val[0], exponent );
      out->val[0] = val1;
    }
    if (mode == 2) {
      out->val[0] = Max( val0, val1 );
    }
  }
  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "sg_createAlloc"
/*!
  @par Revision history:
  - 21.12.2005, c
*/
int32 sg_createAlloc( SurfaceGeometry **p_obj,
		      int32 nFa, int32 nQP, int32 dim, int32 nFP )
{
  SurfaceGeometry *obj;

  obj = allocMem( SurfaceGeometry, 1 );
  obj->nFa = nFa;
  obj->nQP = nQP;
  obj->dim = dim;
  obj->nFP = nFP;
  fmf_createAlloc( &(obj->normal), nFa, nQP, dim, 1 );
  fmf_createAlloc( &(obj->det), nFa, nQP, 1, 1 );
  fmf_createAlloc( &(obj->area), nFa, 1, 1, 1 );
  obj->mode = -1;

  *p_obj = obj;

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "sg_freeDestroy"
/*!
  @par Revision history:
  - 21.12.2005, c
  - 04.05.2007
*/
int32 sg_freeDestroy( SurfaceGeometry **p_obj )
{
  SurfaceGeometry *obj = *p_obj;

  if (!obj) return( RET_OK );

  if (obj->normal) 
    fmf_freeDestroy( &(obj->normal) );
  if (obj->det) 
    fmf_freeDestroy( &(obj->det) );
  if (obj->area) 
    fmf_freeDestroy( &(obj->area) );
  if (obj->bfBGM) 
    fmf_freeDestroy( &(obj->bfBGM) );
  freeMem( *p_obj );

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "sg_print"
/*!
  @par Revision history:
  - c: 21.12.2005, r: 18.01.2008
*/
int32 sg_print( SurfaceGeometry *obj, FILE *file, int32 mode )
{
  int32 ii;

  fprintf( file, "SurfaceGeometry: mode %d, nFa %d, nQP %d, dim: %d, nFP: %d\n",
	   obj->mode, obj->nFa, obj->nQP, obj->dim, obj->nFP );
  fprintf( file, "totalArea: %.5f\n", obj->totalArea );

  for (ii = 0; ii < obj->det->nCell; ii++) {
    FMF_SetCell( obj->normal, ii );
    FMF_SetCell( obj->det, ii );
    FMF_SetCell( obj->area, ii );
    
    fprintf( file, "%d normal:\n", ii );
    fmf_print( obj->normal, file, mode );
    
    fprintf( file, "%d det:\n", ii );
    fmf_print( obj->det, file, mode );
    
    fprintf( file, "%d area:\n", ii );
    fmf_print( obj->area, file, mode );

    if (obj->bfBGM) {
      FMF_SetCell( obj->bfBGM, ii );
      fprintf( file, "%d bfBGM:\n", ii );
      fmf_print( obj->bfBGM, file, mode );
    }

  }

  return( RET_OK );
}


#undef __FUNC__
#define __FUNC__ "sg_describe"
/*!
  @par Revision history:
  - 21.12.2005, c
  - 05.09.2006
*/
int32 sg_describe( SurfaceGeometry *obj,
		   float64 *coorIn, int32 nNod, int32 dim,
		   int32 *fconn, int32 nFa, int32 nFP,
		   FMField *bfGR, FMField *weight )
{
  int32 ii, pos, inod, idim, iqp, nQP, ret = RET_OK;
  float64 c1, c2, c3, det;
  float64 *jmat;
  FMField *faceCoor = 0, *mtxRMS = 0;
  
  nQP = bfGR->nLev;
  if (!((nFa == obj->nFa) &&
	(dim == obj->dim) &&
	(nQP == obj->nQP) &&
	(nFP == obj->nFP))) {
    output( "nNod: %d, dim: %d, nFa: %d, nFP: %d\n",  nNod, dim, nFa, nFP );
    fmf_print( obj->normal, stdout, 1 );
    fmf_print( bfGR, stdout, 1 );
    fmf_print( weight, stdout, 1 );
    errput( "size mismatch!\n" );
    return( RET_Fail );
  }

/*    output( "%d %d %d %d\n", dim, nQP, nFP, nNod ); */
  fmf_createAlloc( &faceCoor, 1, 1, nFP, dim );
  fmf_createAlloc( &mtxRMS, 1, nQP, dim - 1, dim );

  for (ii = 0; ii < nFa; ii++) {
    FMF_SetCell( obj->normal, ii );
    FMF_SetCell( obj->det, ii );
    FMF_SetCell( obj->area, ii );

    for (inod = 0; inod < nFP; inod++) {
      pos = dim*fconn[inod];
      for (idim = 0; idim < dim; idim++ ) {
	faceCoor->val[dim*inod+idim] = coorIn[idim+pos];
      }
    }

/*      fmf_print( faceCoor, stdout, 0 ); */
    fmf_mulAB_n1( mtxRMS, bfGR, faceCoor );
    /* Surface jacobian and normal. */
    switch (dim) {
    case 2:
      /* dl = \norma{dx} = sqrt( dx^2 + dy^2 ) */
      for (iqp = 0; iqp < nQP; iqp++) {
	jmat = FMF_PtrLevel( mtxRMS, iqp );
	c1 = jmat[0];
	c2 = jmat[1];
	det = sqrt( c1*c1 + c2*c2 );
	obj->det->val[iqp] = det * weight->val[iqp];
	/* Unit outward normal. */
	obj->normal->val[2*(iqp)+0] = c2 / det;
	obj->normal->val[2*(iqp)+1] = -c1 / det;
      }
      break;
    case 3:
      /* dS = \norma{dx x dy} */
      for (iqp = 0; iqp < nQP; iqp++) {
	jmat = FMF_PtrLevel( mtxRMS, iqp );
	c1 = jmat[1] * jmat[5] - jmat[4] * jmat[2];
	c2 = jmat[0] * jmat[5] - jmat[3] * jmat[2];
	c3 = jmat[0] * jmat[4] - jmat[1] * jmat[3];
	det = sqrt( c1*c1 + c2*c2 + c3*c3 );
	/*  	  printf( "s: %f %f %f %f\n", c1, -c2, c3, det ); */
	obj->det->val[iqp] = det * weight->val[iqp];
	/* Unit outward normal. */
	obj->normal->val[3*(iqp)+0] = c1 / det;
	obj->normal->val[3*(iqp)+1] = -c2 / det;
	obj->normal->val[3*(iqp)+2] = c3 / det;
      }
      break;
    default:
      errput( ErrHead "ERR_Switch\n" );
    }

    // Element area.
    geme_elementVolume( obj->area->val, obj->det->val, nQP );
    obj->totalArea += obj->area->val[0];

  
    fconn += nFP;
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &faceCoor );
  fmf_freeDestroy( &mtxRMS );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "sg_integrate"
/*!
  For scalar input p: \int_{\Gamma} p dS
  For vector input v: \int_{\Gamma} v n dS

  @par Revision history:
  - 24.04.2007, c
*/
int32 sg_integrate( SurfaceGeometry *obj, FMField *out, FMField *in )
{
  int32 dim, nQP, iel, ret = RET_OK;
  FMField *vn = 0;

  dim = obj->normal->nRow;
  nQP = obj->normal->nLev;

  if (in->nRow == 1) {
    for (iel = 0; iel < obj->det->nCell; iel++) {
      FMF_SetCell( obj->det, iel );
      FMF_SetCell( in, iel );
      FMF_SetCell( out, iel );

      fmf_sumLevelsMulF( out, in, obj->det->val );
      ERR_CheckGo( ret );
    }
  } else if (in->nRow == dim) {
    fmf_createAlloc( &vn, 1, nQP, 1, 1 );

    for (iel = 0; iel < obj->det->nCell; iel++) {
      FMF_SetCell( obj->normal, iel );
      FMF_SetCell( obj->det, iel );
      FMF_SetCell( in, iel );
      FMF_SetCell( out, iel );

      fmf_mulATB_nn( vn, in, obj->normal );
/*       fmf_mulC( vn, -1.0 ); */

      fmf_sumLevelsMulF( out, vn, obj->det->val );
      ERR_CheckGo( ret );
    }
  } else {
    errput( ErrHead "ERR_Switch\n" );
  }

 end_label:
  fmf_freeDestroy( &vn );

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "sg_integrateChunk"
/*!
  For scalar input p: \int_{\Gamma} p dS
  For vector input v: \int_{\Gamma} v n dS

  @par Revision history:
  - 01.11.2007, c
*/
int32 sg_integrateChunk( SurfaceGeometry *obj, FMField *out, FMField *in,
			 int32 *elList, int32 elList_nRow )
{
  int32 dim, nQP, iel, ii, ret = RET_OK;
  FMField *vn = 0;

  dim = obj->normal->nRow;
  nQP = obj->normal->nLev;

  if (in->nRow == 1) {
    for (ii = 0; ii < elList_nRow; ii++) {
      iel = elList[ii];
      FMF_SetCell( obj->det, iel );
      FMF_SetCell( in, ii );
      FMF_SetCell( out, ii );

      fmf_sumLevelsMulF( out, in, obj->det->val );
      ERR_CheckGo( ret );
    }
  } else if (in->nRow == dim) {
    fmf_createAlloc( &vn, 1, nQP, 1, 1 );

    for (ii = 0; ii < elList_nRow; ii++) {
      iel = elList[ii];
      FMF_SetCell( obj->normal, iel );
      FMF_SetCell( obj->det, iel );
      FMF_SetCell( in, ii );
      FMF_SetCell( out, ii );

      fmf_mulATB_nn( vn, in, obj->normal );
/*       fmf_mulC( vn, -1.0 ); */

      fmf_sumLevelsMulF( out, vn, obj->det->val );
/*       fmf_print( vn, stdout, 0 ); */
/*       fmf_print( out, stdout, 0 ); */
/*       sys_pause(); */
      ERR_CheckGo( ret );
    }
  } else {
    errput( ErrHead "ERR_Switch\n" );
  }

 end_label:
  fmf_freeDestroy( &vn );

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "sg_evaluateBFBGM"
/*!
  @par Revision history:
  - 04.05.2007, c
  - 30.05.2007
*/
int32 sg_evaluateBFBGM( SurfaceGeometry *obj, FMField *bfBGR,
			float64 *coorIn, int32 nNod, int32 dim,
			int32 *fis, int32 nFa, int32 nFP,
			int32 *conn, int32 nEl, int32 nEP )
{
  int32 ii, iel, ifa, inod, idim, pos, nQP, ret = RET_OK;
  FMField *volCoor0 = 0, *mtxRM = 0, *mtxRMI = 0;
  
  nQP = obj->normal->nLev;

/*    output( "%d %d %d %d\n", dim, nQP, nEP, nNod ); */

  fmf_createAlloc( &volCoor0, 1, 1, nEP, dim );
  fmf_createAlloc( &mtxRM, 1, nQP, dim, dim );
  fmf_createAlloc( &mtxRMI, 1, nQP, dim, dim );

  for (ii = 0; ii < nFa; ii++) {
    iel = fis[ii*nFP+0];
    ifa = fis[ii*nFP+1];
    
    FMF_SetCell( obj->bfBGM, ii );
    FMF_SetCell( bfBGR, ifa );

    for (inod = 0; inod < nEP; inod++) {
      pos = dim*conn[nEP*iel+inod];
      for (idim = 0; idim < dim; idim++ ) {
	volCoor0->val[dim*inod+idim] = coorIn[idim+pos];
      }
    }
    fmf_mulAB_n1( mtxRM, bfBGR, volCoor0 );
    geme_invert3x3( mtxRMI, mtxRM );
    fmf_mulAB_nn( obj->bfBGM, mtxRMI, bfBGR );
/*     fmf_mulATBT_1n( mtxRM, volCoor0, bfBGR ); */
/*     geme_invert3x3( mtxRMI, mtxRM ); */
/*     fmf_mulATB_nn( obj->bfBGM, mtxRMI, bfBGR ); */

/*     output( "%d %d %d\n", ii, iel, ifa); */
/*     fmf_print( bfBGR, stdout, 0 ); */
/*     fmf_print( volCoor0, stdout, 0 ); */
/*     fmf_print( obj->bfBGM, stdout, 0 ); */
/*     sys_pause(); */
    
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &volCoor0 );
  fmf_freeDestroy( &mtxRM );
  fmf_freeDestroy( &mtxRMI );

  return( RET_OK );
}
