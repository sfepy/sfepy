#include "refmaps.h"

#undef __FUNC__
#define __FUNC__ "map_print"
/*!
  @par Revision history:
  - 11.10.2005, c
  - 12.10.2005
  - 12.10.2006
*/
int32 map_print( Mapping *obj, FILE *file, int32 mode )
{
  int32 ii;
  char *modes[] = {"volume", "surface", "surface_extra"};

  fprintf( file, "Mapping: mode %s, nEl "FI32", nQP "FI32", dim: "
           FI32", nEP: "FI32"\n",
           modes[obj->mode], obj->nEl, obj->nQP, obj->dim, obj->nEP );
  fprintf( file, "totalVolume: %.5f\n", obj->totalVolume );

  for (ii = 0; ii < obj->det->nCell; ii++) {
    FMF_SetCell( obj->det, ii );
    FMF_SetCell( obj->volume, ii );

    fprintf( file, FI32" det:\n", ii );
    fmf_print( obj->det, file, mode );

    fprintf( file, FI32" volume:\n", ii );
    fmf_print( obj->volume, file, mode );

    if ((obj->mode == MM_Volume) || (obj->mode == MM_SurfaceExtra)){
      FMF_SetCell( obj->bfGM, ii );

      fprintf( file, FI32 " bfGM:\n", ii );
      fmf_print( obj->bfGM, file, mode );
    } else {
      FMF_SetCell( obj->normal, ii );

      fprintf( file, FI32 " normal:\n", ii );
      fmf_print( obj->normal, file, mode );
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "map_describe"
/*!
  @par Revision history:
  - 11.10.2005, c
  - 12.10.2005
  - 05.06.2006
  - 06.06.2006
*/
int32 map_describe( Mapping *obj,
                    float64 *coorIn, int32 nNod, int32 dim,
                    int32 *conn, int32 nEl, int32 nEP,
                    FMField *bfGR, FMField *ebfGR, FMField *weight )
{
  int32 nQP, ret = RET_OK;

  nQP = bfGR->nLev;
  if (!((nEl == obj->nEl) &&
        (dim == obj->dim) &&
        (nQP == obj->nQP) &&
        (nEP == bfGR->nCol) &&
        (ebfGR->nCol == obj->nEP))) {
    map_print( obj, stdout, 1 );
    errput( "size mismatch!\n" );
    return( RET_Fail );
  }

  if (obj->mode == MM_Volume) {
    ret = _v_describe(obj, coorIn, nNod, dim, conn, nEl, nEP,
                      bfGR, ebfGR, weight);
  } else {
    ret = _s_describe(obj, coorIn, nNod, dim, conn, nEl, nEP,
                      bfGR, weight);
  }

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "_v_describe"
/*!
*/
int32 _v_describe( Mapping *obj,
                   float64 *coorIn, int32 nNod, int32 dim,
                   int32 *conn, int32 nEl, int32 nEP,
                   FMField *bfGR, FMField *ebfGR, FMField *weight )
{
  int32 iel, inod, idim, pos, iqp, nQP, ret = RET_OK;
  FMField *mtxMR = 0, *mtxMRI = 0, *coor = 0;

  nQP = bfGR->nLev;

  fmf_createAlloc( &mtxMR, 1, nQP, dim, dim );
  fmf_createAlloc( &mtxMRI, 1, nQP, dim, dim );
  fmf_createAlloc( &coor, 1, 1, nEP, dim );

  obj->totalVolume = 0.0;

  for (iel = 0; iel < obj->bfGM->nCell; iel++) {
    FMF_SetCell( obj->bfGM, iel );
    FMF_SetCell( obj->det, iel );
    FMF_SetCell( obj->volume, iel );
    FMF_SetCellX1( ebfGR, iel );

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
        errput( "warp violation %e at (iel: "FI32", iqp: "FI32")!\n",
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
    fmf_mulATB_nn( obj->bfGM, mtxMRI, ebfGR );

    conn += nEP;

    ERR_CheckGo( ret );
  }
 end_label:
  fmf_freeDestroy( &mtxMR );
  fmf_freeDestroy( &mtxMRI );
  fmf_freeDestroy( &coor );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "_s_describe"
/*!
  @par Revision history:
  - 21.12.2005, c
  - 05.09.2006
*/
int32 _s_describe( Mapping *obj,
                   float64 *coorIn, int32 nNod, int32 dim,
                   int32 *fconn, int32 nFa, int32 nFP,
                   FMField *bfGR, FMField *weight )
{
  int32 ii, pos, inod, idim, iqp, nQP, ret = RET_OK;
  float64 c1, c2, c3, det;
  float64 *jmat;
  FMField *faceCoor = 0, *mtxRMS = 0;

  nQP = bfGR->nLev;

/*    output( "%d %d %d %d\n", dim, nQP, nFP, nNod ); */
  fmf_createAlloc( &faceCoor, 1, 1, nFP, dim );
  fmf_createAlloc( &mtxRMS, 1, nQP, dim - 1, dim );

  for (ii = 0; ii < nFa; ii++) {
    FMF_SetCell( obj->normal, ii );
    FMF_SetCell( obj->det, ii );
    FMF_SetCell( obj->volume, ii );

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
        /*        printf( "s: %f %f %f %f\n", c1, -c2, c3, det ); */
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

    // Face area.
    geme_elementVolume( obj->volume->val, obj->det->val, nQP );
    obj->totalVolume += obj->volume->val[0];

    fconn += nFP;
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &faceCoor );
  fmf_freeDestroy( &mtxRMS );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "map_integrate"
/*!
  Volume mapping: \int_{\Omega} u dV

  Surface mapping:
  - for scalar input p: \int_{\Gamma} p dS
  - for vector input v or mode < 3: \int_{\Gamma} v n dS

  @par Revision history:
  - 15.12.2005, c
*/
int32 map_integrate( Mapping *obj, FMField *out, FMField *in,
                     int32 mode )
{
  int32 nQP, iel, ret = RET_OK;
  FMField *vn = 0;

  if ((obj->mode == MM_Volume) || (mode < 3) || (in->nRow == 1)) {
    for (iel = 0; iel < obj->bfGM->nCell; iel++) {
      FMF_SetCell( obj->det, iel );
      FMF_SetCell( in, iel );
      FMF_SetCell( out, iel );
      fmf_sumLevelsMulF( out, in, obj->det->val );
      if (mode == 1) {
        FMF_SetCell( obj->volume, iel );
        fmf_mulC( out, 1.0 / obj->volume->val[0] );
      }
      ERR_CheckGo( ret );
    }
  } else if (in->nRow == obj->dim) {
    nQP = obj->normal->nLev;

    fmf_createAlloc( &vn, 1, nQP, 1, 1 );

    for (iel = 0; iel < obj->det->nCell; iel++) {
      FMF_SetCell( obj->normal, iel );
      FMF_SetCell( obj->det, iel );
      FMF_SetCell( in, iel );
      FMF_SetCell( out, iel );

      fmf_mulATB_nn( vn, in, obj->normal );

      fmf_sumLevelsMulF( out, vn, obj->det->val );
      if (mode == 4) {
        FMF_SetCell( obj->volume, iel );
        fmf_mulC( out, 1.0 / obj->volume->val[0] );
      }
      ERR_CheckGo( ret );
    }
  } else {
    errput( ErrHead "ERR_Switch\n" );
    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &vn );

  return( ret );
}

#undef __FUNC__
#define __FUNC__ "map_getElementDiameters"
/*!
  @par Revision history:
  - 09.01.2006, c
  - 10.01.2006
*/
int32 map_getElementDiameters( Mapping *obj, FMField *out,
                               int32 *edges, int32 edges_nRow, int32 edges_nCol,
                               float64 *coorIn, int32 nNod, int32 dim,
                               int32 *conn, int32 nEl, int32 nEP,
                               int32 *elList, int32 elList_nRow,
                               int32 mode )
{
  int32 ii, ie, id, iel, nd;
  float64 val0 = 0.0, val1 = 0.0, vv, aux = 0.0, exponent;

  if (obj->mode != MM_Volume) {
    errput( ErrHead "only for volume mappings!\n" );
    return( RET_Fail );
  }

  if ((mode < 0) && (mode > 2)) {
    errput( ErrHead "ERR_Switch\n" );
    return( RET_Fail );
  }

/*   output( "%d %d %d %d %d %d %d\n", */
/*        edges_nRow, edges_nCol, nNod, dim, nEl, nEP, elList_nRow ); */

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
/*        output( "%d %d %d %d %f %f %f\n", ii, iel, ie, id, aux, */
/*                coorIn[dim*conn[nEP*iel+edges[2*ie+1]]+id], */
/*                coorIn[dim*conn[nEP*iel+edges[2*ie+0]]+id] ); */
/*        sys_pause(); */
          vv += aux * aux;
        }
/*      output("%f\n", sqrt(vv)); */
        val0 = Max( val0, vv );
        out->val[0] = val0;
      }
    }
    if ((mode == 1) || (mode == 2)) {
      FMF_SetCell( obj->volume, ii );
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
#define __FUNC__ "map_evaluateBFBGM"
/*!
  @par Revision history:
  - 04.05.2007, c
  - 30.05.2007
*/
int32 map_evaluateBFBGM( Mapping *obj, FMField *bfBGR, FMField *ebfBGR,
                         float64 *coorIn, int32 nNod, int32 dim,
                         int32 *fis, int32 nFa, int32 nFP,
                         int32 *conn, int32 nEl, int32 nEP )
{
  int32 ii, iel, ifa, inod, idim, pos, nQP, ret = RET_OK;
  FMField *volCoor0 = 0, *mtxRM = 0, *mtxRMI = 0;

  if (obj->mode != MM_SurfaceExtra) {
    errput( ErrHead "only for surface extra mappings!\n" );
    return( RET_Fail );
  }

  nQP = obj->normal->nLev;

  fmf_createAlloc( &volCoor0, 1, 1, nEP, dim );
  fmf_createAlloc( &mtxRM, 1, nQP, dim, dim );
  fmf_createAlloc( &mtxRMI, 1, nQP, dim, dim );

  for (ii = 0; ii < nFa; ii++) {
    iel = fis[ii*nFP+0];
    ifa = fis[ii*nFP+1];

    FMF_SetCell( obj->bfGM, ii );
    FMF_SetCell( bfBGR, ifa );
    FMF_SetCell( ebfBGR, ifa );

    for (inod = 0; inod < nEP; inod++) {
      pos = dim*conn[nEP*iel+inod];
      for (idim = 0; idim < dim; idim++ ) {
        volCoor0->val[dim*inod+idim] = coorIn[idim+pos];
      }
    }
    fmf_mulAB_n1( mtxRM, bfBGR, volCoor0 );
    geme_invert3x3( mtxRMI, mtxRM );
    fmf_mulAB_nn( obj->bfGM, mtxRMI, ebfBGR );

    ERR_CheckGo( ret );
  }

 end_label:
  fmf_freeDestroy( &volCoor0 );
  fmf_freeDestroy( &mtxRM );
  fmf_freeDestroy( &mtxRMI );

  return( ret );
}
