#include "fem.h"
#include "geommech.h"
#include "sort.h"

#undef __FUNC__
#define __FUNC__ "assemble_vector"
/*!
  @par Revision history:
  - 21.11.2005, c
  - 26.11.2005
  - 27.11.2005
*/
int32 assemble_vector( FMField *vec, FMField *vecInEls,
		       int32 *iels, int32 iels_len,
		       float64 sign, int32 *conn, int32 nEl, int32 nEP )
{
  int32 ii, iel, ir, irg;
  int32 *pconn;
  float64 *val;

/*   output( "%f %d %d %d\n", sign, iels_len, nEl, nEP ); */

  val = FMF_PtrFirst( vec );

  for (ii = 0; ii < iels_len; ii++) {
    iel = iels[ii];
    FMF_SetCell( vecInEls, ii );

    pconn = conn + nEP * iel;
    for (ir = 0; ir < nEP; ir++) {
      irg = pconn[ir];
      if (irg < 0) continue;
      
/*       output( "%d %d %d\n", iel, ir, irg ); */
      val[irg] += sign * vecInEls->val[ir];
    }
/*     fmf_print( vecInEls, stdout, 0 ); */
/*     sys_pause(); */
  }
/*   fmf_print( vec, stdout, 0 ); */

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "assemble_vector_complex"
int32 assemble_vector_complex( FMField *vec_r, FMField *vec_i,
			       FMField *vecInEls_r, FMField *vecInEls_i,
			       int32 *iels, int32 iels_len,
			       float64 sign_r, float64 sign_i,
			       int32 *conn, int32 nEl, int32 nEP )
{
  int32 ii, iel, ir, irg, stride, stride2;
  int32 *pconn;
  float64 *val_r, *val_i;

  stride = vec_r->stride;
  stride2 = vecInEls_r->stride;

  val_r = FMF_PtrFirst( vec_r );
  val_i = FMF_PtrFirst( vec_i );

  for (ii = 0; ii < iels_len; ii++) {
    iel = iels[ii];
    FMF_SetCell( vecInEls_r, ii );
    FMF_SetCell( vecInEls_i, ii );

    pconn = conn + nEP * iel;
    for (ir = 0; ir < nEP; ir++) {
      irg = pconn[ir];
      if (irg < 0) continue;
      
/*       output( "%d %d %d\n", iel, ir, irg ); */
      val_r[stride*irg] += sign_r * vecInEls_r->val[stride2*ir]
   	- sign_i * vecInEls_i->val[stride2*ir];
      val_i[stride*irg] += sign_r * vecInEls_i->val[stride2*ir]
	+ sign_i * vecInEls_r->val[stride2*ir];
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "assemble_matrix"
/*!
  Requires a CSR matrix.

  @par Revision history:
  - 27.11.2005, c
  - 15.12.2005
*/
int32 assemble_matrix( FMField *mtx,
		       int32 *prows, int32 prows_len,
		       int32 *cols, int32 cols_len,
		       FMField *mtxInEls,
		       int32 *iels, int32 iels_len, float64 sign,
		       int32 *connR, int32 nElR, int32 nEPR,
		       int32 *connC, int32 nElC, int32 nEPC )
{
  int32 ii, iel, ir, ic, irg, icg, is, iloc, found;
  int32 *pconnR, *pconnC;
  float64 *val;

/*   output( "%f %d %d %d %d %d %d %d\n", */
/* 	  sign, iels_len, prows_len, cols_len, nElR, nEPR, nElC, nEPC ); */

  val = FMF_PtrFirst( mtx );

  for (ii = 0; ii < iels_len; ii++) {
    iel = iels[ii];
    FMF_SetCell( mtxInEls, ii );

    pconnR = connR + nEPR * iel;
    pconnC = connC + nEPC * iel;
    
    for (ir = 0; ir < nEPR; ir++) {
      irg = pconnR[ir];
      if (irg < 0) continue;

      for (ic = 0; ic < nEPC; ic++) {
	icg = pconnC[ic];
	if (icg < 0) continue;
	iloc = nEPC * ir + ic;

/* 	output( "%d %d %d %d %d %d\n", iel, ir, ic, irg, icg, iloc ); */
	/* Try bsearch instead here... */
	found = 0;
	for (is = prows[irg]; is < prows[irg+1]; is++) {
	  if (cols[is] == icg) {
	    val[is] += sign * mtxInEls->val[iloc];
	    found = 1;
	    break;
	  }
	}
	if (!found) {
	  errput( "matrix item (%d,%d) does not exist\n", irg, icg );
	  return( RET_Fail );
	}
      }
    }
/*     fmf_print( mtxInEls, stdout, 0 ); */
/*     sys_pause(); */
  }
/*   fmf_print( mtx, stdout, 0 ); */

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "assemble_matrix_complex"
/*!
  Requires a CSR matrix.
*/
int32 assemble_matrix_complex( FMField *mtx_r, FMField *mtx_i,
			       int32 *prows, int32 prows_len,
			       int32 *cols, int32 cols_len,
			       FMField *mtxInEls_r, FMField *mtxInEls_i,
			       int32 *iels, int32 iels_len,
			       float64 sign_r, float64 sign_i,
			       int32 *connR, int32 nElR, int32 nEPR,
			       int32 *connC, int32 nElC, int32 nEPC )
{
  int32 ii, iel, ir, ic, irg, icg, is, iloc, found, stride, stride2;
  int32 *pconnR, *pconnC;
  float64 *val_r, *val_i;

  stride = mtx_r->stride;
  stride2 = mtxInEls_r->stride;

  val_r = FMF_PtrFirst( mtx_r );
  val_i = FMF_PtrFirst( mtx_i );

  for (ii = 0; ii < iels_len; ii++) {
    iel = iels[ii];
    FMF_SetCell( mtxInEls_r, ii );
    FMF_SetCell( mtxInEls_i, ii );

    pconnR = connR + nEPR * iel;
    pconnC = connC + nEPC * iel;
    
    for (ir = 0; ir < nEPR; ir++) {
      irg = pconnR[ir];
      if (irg < 0) continue;

      for (ic = 0; ic < nEPC; ic++) {
	icg = pconnC[ic];
	if (icg < 0) continue;
	iloc = stride2 * (nEPC * ir + ic);

/* 	output( "%d %d %d %d %d %d\n", iel, ir, ic, irg, icg, iloc ); */
	/* Try bsearch instead here... */
	found = 0;
	for (is = prows[irg]; is < prows[irg+1]; is++) {
	  if (cols[is] == icg) {
	    val_r[stride*is] += sign_r * mtxInEls_r->val[iloc]
	      - sign_i * mtxInEls_i->val[iloc];
	    val_i[stride*is] += sign_r * mtxInEls_i->val[iloc]
	      + sign_i * mtxInEls_r->val[iloc];
	    found = 1;
	    break;
	  }
	}
	if (!found) {
	  errput( "matrix item (%d,%d) does not exist\n", irg, icg );
	  return( RET_Fail );
	}
      }
    }
/*     fmf_print( mtxInEls, stdout, 0 ); */
/*     sys_pause(); */
  }
/*   fmf_print( mtx, stdout, 0 ); */

  return( RET_OK );
}

int32 compareI32( const void *a, const void *b )
{
  int32 i1, i2;

  i1 = *((int32 *) a);
  i2 = *((int32 *) b);

  return( i1 - i2 );
}

#undef __FUNC__
#define __FUNC__ "mesh_nodInElCount"
/*!
  @par Revision history:
  - 21.11.2003, c
  - 23.11.2003
*/
int32 mesh_nodInElCount( int32 *p_niecMax, int32 *niec,
			 int32 nNod, int32 nGr, int32 *nEl,
			 int32 *nEP, int32 **conn )
{
  int32 ig, iel, iep, in, niecMax;
  int32 *pconn;

  memset( niec, 0, (nNod + 1) * sizeof( int32 ) );
  for (ig = 0; ig < nGr; ig++) {
    for (iel = 0; iel < nEl[ig]; iel++) {
      pconn = conn[ig] + nEP[ig] * iel;
      for (iep = 0; iep < nEP[ig]; iep++) {
	niec[1+pconn[iep]]++;
/* 	output( "%d %d %d\n", iep, niec[1+pconn[iep]], pconn[iep] ); */
      }
    }
  }

  niec[0] = 0;
  niecMax = 0;
  for (in = 0; in <= nNod; in++) {
/*     output( "%d %d\n", in, niec[in] ); */
    niecMax = Max( niecMax, niec[in] );
  }
  *p_niecMax = niecMax;

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "mesh_graph"
/*!
  @par Revision history:
  - 23.05.2003, c
  - 26.05.2003
  - 27.05.2003
  - 28.05.2003
  - 21.11.2003 former mesh_meshGraph()
  - 23.11.2003
  - 01.03.2004
  - 03.03.2005
  - 07.02.2006
*/
int32 mesh_graph( int32 *p_nnz, int32 **p_prow, int32 **p_icol,
		  int32 nRow, int32 nCol, int32 nGr, int32 *nEl,
		  int32 *nEPR, int32 **connR, int32 *nEPC, int32 **connC )
{
  int32 in, ii, ip, ig, iel, iep, ir, ic, nn, np, pr,
    niecMaxR, nEPMaxC, nUnique, iir, iic, found;
  int32 *niec, *pconnR, *pconnC, *eonlist, *nir, *nods, *icol;


/*   output( "%d %d %d %d %d %d\n", nRow, nCol, nGr, nEl[0], nEPR[0], nEPC[0] ); */

  /* Get niec (= nodes in elements count) for rows. */
  niec = alloc_mem( int32, nRow + 1 );
  mesh_nodInElCount( &niecMaxR, niec, nRow, nGr, nEl, nEPR, connR );
/*   output( "%d\n", niecMaxR ); */

  /* Cummulative sum. */
  for (in = 0; in < nRow; in++) {
    niec[in+1] += niec[in];
  }

/*    output( "00\n" ); */

  /* eon = elements of nodes */
  nn = 0;
  nEPMaxC = 0;
  for (ig = 0; ig < nGr; ig++) {
    nn += nEPR[ig] * nEl[ig];
    nEPMaxC = Max( nEPMaxC, nEPC[ig] );
  }
  eonlist = alloc_mem( int32, 2 * nn );

  /* nir is just a buffer here. */
  nir = alloc_mem( int32, nRow + 1 );
  memset( nir, 0, (nRow + 1) * sizeof( int32 ) );

/*    output( "1\n" ); */

  /* Get list of elements each row node is in. */
  for (ig = 0; ig < nGr; ig++) {
    for (iel = 0; iel < nEl[ig]; iel++) {
      pconnR = connR[ig] + nEPR[ig] * iel;
      for (iep = 0; iep < nEPR[ig]; iep++) {
	np = pconnR[iep];
	if (np >= 0) {
	  eonlist[2*(niec[np]+nir[np])+0] = iel;
	  eonlist[2*(niec[np]+nir[np])+1] = ig;
/*  	output( "  %d %d %d %d\n", np, eonlist[2*(niec[np]+nir[np])+0], */
/*  		   eonlist[2*(niec[np]+nir[np])+1], nir[np] ); */
	  nir[np]++;
	}
      }
    }
  }

/*    output( "2\n" ); */
 
  /* nir = number in row. */
  memset( nir, 0, (nRow + 1) * sizeof( int32 ) );

  /* List of column nodes for each row node. */
/*   output( "%d, %d\n", nEPMaxC, niecMaxR * nEPMaxC ); */
  nods = alloc_mem( int32, niecMaxR * nEPMaxC );

  nn = 0;
  for (in = 0; in < nRow; in++) {
    ii = 0;
/*      output( "%d\n", in ); */
    for (ip = niec[in]; ip < niec[in+1]; ip++) {
      iel = eonlist[2*(ip)+0];
      ig = eonlist[2*(ip)+1];
/*        output( " %d %d %d\n", ip, ig, iel ); */
      for (iep = 0; iep < nEPC[ig]; iep++) {
	np = connC[ig][nEPC[ig]*iel+iep];
	if (np >= 0) {
	  nods[ii] = np;
/*  	output( "  %d %d\n", ii, nods[ii] ); */
	  ii++;
	}
      }
    }
/*     output( "%d\n", ii ); */

    if (ii > 0) {
/*       qsort( nods, ii, sizeof( int32 ), &compareI32 ); */
      int32_quicksort( nods, ii, 0 );
      nUnique = 1;
      for (ir = 0; ir < (ii - 1); ir++) {
	if (nods[ir] != nods[ir+1]) {
	  nUnique++;
	}
      }
    } else {
      nUnique = 0;
    }
    nn += nUnique;
/*      output( " -> %d\n", nUnique ); */

    nir[in] = nUnique;
  }
  
/*    output( "3\n" ); */

  *p_nnz = nn;
  *p_prow = niec;
  icol = *p_icol = alloc_mem( int32, nn );

  /* Fill in *p_prow. */
  niec[0] = 0;
  for (in = 0; in < nRow; in++) {
    niec[in+1] = niec[in] + nir[in];
/*      output( " %d\n", niec[in+1] ); */
  }

/*   output( "4\n" ); */
  /* Fill in *p_icol (sorted). */
  memset( nir, 0, (nRow + 1) * sizeof( int32 ) );
  for (ig = 0; ig < nGr; ig++) {
/*     output( "ig %d\n", ig ); */
    for (iel = 0; iel < nEl[ig]; iel++) {
      pconnR = connR[ig] + nEPR[ig] * iel;
      pconnC = connC[ig] + nEPC[ig] * iel;
      for (ir = 0; ir < nEPR[ig]; ir++) {
	iir = pconnR[ir];
	if (iir < 0) continue;
	pr = niec[iir];
/*  	output( " %d %d %d\n", iir, pr, niec[iir+1] - pr ); */
	for (ic = 0; ic < nEPC[ig]; ic++) {
	  iic = pconnC[ic];
	  if (iic < 0) continue;
/*  	  output( "   %d %d\n", iic, nir[iir] ); */
	  /* This is a bottle-neck! */
	  found = 0;
	  for (ii = pr; ii < (pr + nir[iir]); ii++) {
	    if (icol[ii] == iic) {
	      found = 1;
	      break;
	    }
	  }
/*  	  output( "  ? %d\n", found ); */
	  if (!found) {
	    if (nir[iir] < (niec[iir+1] - pr)) {
	      icol[pr+nir[iir]] = iic;
	      nir[iir]++;
/*  	      output( "  + %d %d\n", nir[iir], niec[iir+1] - pr ); */
	    } else {
	      output( "  %d %d\n", nir[iir], niec[iir+1] - pr );
	      errput( "ERR_VerificationFail\n" );
	    }
	  }
	}
/* 	qsort( icol + pr, nir[iir], sizeof( int32 ), &compareI32 ); */
	int32_quicksort( icol + pr, nir[iir], 0 );
      }
    }
  }

/*   output( "5\n" ); */

  free_mem( nods );
  free_mem( nir );
  free_mem( eonlist );

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "raw_graph"
/*!
  @par Revision history:
  - 27.11.2005, c
  - 19.02.2007
*/
int32 raw_graph( int32 *p_nRow, int32 **p_prow,
		int32 *p_nnz, int32 **p_icol,
		int32 nRow, int32 nCol, int32 nGr,
		int32 *nElR, int32 *nEPR, int32 **connR,
		int32 *nElC, int32 *nEPC, int32 **connC )
{
  int32 ii;

  for (ii = 0; ii < nGr; ii++) {
    if (nElR[ii] != nElC[ii]) {
      errput( "row and col connectivities nEl: %d == %d\n",
	      nElR[ii], nElC[ii] );
      return( RET_Fail );
    }
  }

  mesh_graph( p_nnz, p_prow, p_icol,
	      nRow, nCol, nGr, nElR, nEPR, connR, nEPC, connC );
  *p_nRow = nRow + 1;

  return( RET_OK );
}

const int32 mapTo2D1[8] = {0, 0,
			   1, 0,
			   1, 1,
			   0, 1};
const int32 mapTo3D1[3*8] = {0, 0, 0,
			    1, 0, 0,
			    1, 1, 0,
			    0, 1, 0,
			    0, 0, 1,
			    1, 0, 1,
			    1, 1, 1,
			    0, 1, 1};

const float64  lag1x[2] = {-0.5, 0.5};

float64 lagrange1( float64 x, int32 which )
{
  switch( which ) {
  case 0:
    return( 0.5 * (1.0 - x) );
    break;
  case 1:
    return( 0.5 * (1.0 + x) );
    break;
  default:
    errput( "lagrange1(): wrong function number!" );
    return( 0.0 );
  }
}

float64 lagrange1x( float64 x, int32 which )
{
  switch( which ) {
  case 0:
  case 1:
    return( lag1x[which] );
    break;
  default:
    errput( "lagrange1x(): wrong function number!" );
    return( 0.0 );
  }
}

float64 baseBiL( float64 x, float64 y, int32 which )
{
  int32 i, j;

  if ((which < 0) || (which > 3)) {
    errput( "baseBiL(): wrong function number!" );
    return( 0.0 );
  }
  i = mapTo2D1[2*which];
  j = mapTo2D1[2*which+1];
  return( lagrange1( x, i ) * lagrange1( y, j ) );
}

float64 baseBiLx( float64 x, float64 y, int32 which )
{
  int32 i, j;

  if ((which < 0) || (which > 3)) {
    errput( "baseBiLx(): wrong function number!" );
    return( 0.0 );
  }
  i = mapTo2D1[2*which];
  j = mapTo2D1[2*which+1];
  return( lag1x[i] * lagrange1( y, j ) );
}

float64 baseBiLy( float64 x, float64 y, int32 which )
{
  int32 i, j;

  if ((which < 0) || (which > 3)) {
    errput( "baseBiLy(): wrong function number!" );
    return( 0.0 );
  }
  i = mapTo2D1[2*which];
  j = mapTo2D1[2*which+1];
  return( lagrange1( x, i ) * lag1x[j] );
}

float64 baseTriL( float64 x, float64 y, float64 z, int32 which )
{
  int32 i, j, k;

  if ((which < 0) || (which > 7)) {
    errput( "baseTriL(): wrong function number!" );
    return( 0.0 );
  }
  i = mapTo3D1[3*which];
  j = mapTo3D1[3*which+1];
  k = mapTo3D1[3*which+2];
  return( lagrange1( x, i ) * lagrange1( y, j ) * lagrange1( z, k ) );
}

float64 baseTriLx( float64 x, float64 y, float64 z, int32 which )
{
  int32 i, j, k;

  if ((which < 0) || (which > 7)) {
    errput( "baseTriL(): wrong function number!" );
    return( 0.0 );
  }
  i = mapTo3D1[3*which];
  j = mapTo3D1[3*which+1];
  k = mapTo3D1[3*which+2];
  return( lag1x[i] * lagrange1( y, j ) * lagrange1( z, k ) );
}

float64 baseTriLy( float64 x, float64 y, float64 z, int32 which )
{
  int32 i, j, k;

  if ((which < 0) || (which > 7)) {
    errput( "baseTriLy(): wrong function number!" );
    return( 0.0 );
  }
  i = mapTo3D1[3*which];
  j = mapTo3D1[3*which+1];
  k = mapTo3D1[3*which+2];
  return( lagrange1( x, i ) * lag1x[j] * lagrange1( z, k ) );
}

float64 baseTriLz( float64 x, float64 y, float64 z, int32 which )
{
  int32 i, j, k;

  if ((which < 0) || (which > 7)) {
    errput( "baseTriLz(): wrong function number!" );
    return( 0.0 );
  }
  i = mapTo3D1[3*which];
  j = mapTo3D1[3*which+1];
  k = mapTo3D1[3*which+2];
  return( lagrange1( x, i ) * lagrange1( y, j ) * lag1x[k] );
}

int32 eval_lagrange_simplex( FMField *out, FMField *coors,
			     int32 *nodes, int32 nNod, int32 nCol,
			     int32 order, int32 diff,
			     FMField *mtx_i, FMField *bc,
			     int32 suppress_errors, float64 eps )
{
  int32 ii, ir, ic, i1, i2, in, error, n_i1, n_ii;
  int32 n_coor, n_v, dim, cdim, ret = RET_OK;
  float64 val, dval, dd, vv;
  float64 *pout;

  n_coor = coors->nRow;
  n_v = bc->nRow;
  dim = n_v - 1;
  cdim = coors->nCol;

  // Barycentric coordinates.
  for (ic = 0; ic < n_coor; ic++) {
    for (ir = 0; ir < n_v; ir++) {
      val = 0.0;
      for (ii = 0; ii < dim; ii++) {
	val += mtx_i->val[n_v*ir+ii] * coors->val[cdim*ic+ii];
      }
      val += mtx_i->val[n_v*ir+dim];

      error = 0;
      if (val < 0.0) {
	if (val > (-eps)) {
	  val = 0.0;
	} else {
	  error = 1;
	}
      }
      if (val > 1.0) {
	if (val < (1.0 + eps)) {
	  val = 1.0;
	} else {
	  error = 1;
	}
      }

      if ((error) && (!(suppress_errors))) {
	errput("quadrature point outside of element!\n");
      }

      bc->val[n_coor*ir+ic] = val;

      ERR_CheckGo( ret );
    }
  } 

  if (!diff) {
    fmf_fillC(out, 1.0);
    
    for (ic = 0; ic < n_coor; ic++) {
      pout = FMF_PtrLevel(out, ic);

      for (in = 0; in < nNod; in++) {

	for (i1 = 0; i1 < n_v; i1++) {
	  n_i1 = nodes[nCol*in+i1];
	  /* printf("%d %d \n", in, n_i1); */
	  for (i2 = 0; i2 < n_i1; i2++) {
	    pout[in] *= (order * bc->val[n_coor*i1+ic] - i2) / (i2 + 1.0);
	  }
	}
      }
    }
  } else {
    fmf_fillC(out, 0.0);

    for (ic = 0; ic < n_coor; ic++) {
      pout = FMF_PtrLevel(out, ic);

      for (in = 0; in < nNod; in++) {

	for (ii = 0; ii < n_v; ii++) {
	  vv = 1.0;
	  
	  for (i1 = 0; i1 < n_v; i1++) {
	    if (i1 == ii) continue;
	    n_i1 = nodes[nCol*in+i1];
	    
	    for (i2 = 0; i2 < n_i1; i2++) {
	      vv *= (order * bc->val[n_coor*i1+ic] - i2) / (i2 + 1.0);
	    }

	  }

	  dval = 0.0;
	  n_ii = nodes[nCol*in+ii];
	  for (i1 = 0; i1 < n_ii; i1++) {
	    dd = 1.0;

	    for (i2 = 0; i2 < n_ii; i2++) {
	      if (i1 == i2) continue;

	      dd *= (order * bc->val[n_coor*ii+ic] - i2) / (i2 + 1.0);
	    }
	    dval += dd * order / (i1 + 1.0);
	  }

	  for (ir = 0; ir < dim; ir++) {
	    pout[nNod*ir+in] += vv * dval * mtx_i->val[n_v*ii+ir];
	  }
	}
      }
    }
  }

 end_label:
 
  return( ret );
}

int32 eval_lagrange_tensor_product( FMField *out, FMField *coors,
				    int32 *nodes, int32 nNod, int32 nCol,
				    int32 order, int32 diff,
				    FMField *mtx_i, FMField *bc, FMField *base1d,
				    int32 suppress_errors, float64 eps )
{
  int32 ii, id, im, ic, nr, nc, dim, ret = RET_OK;
  int32 *pnodes = 0;
  FMField c1[1];
    
  c1->nAlloc = -1;
  dim = coors->nCol;

  fmf_fillC( out, 1.0 );

  if (!diff) {
    for (ii = 0; ii < dim; ii++) {
      // slice [:,2*ii:2*ii+2]
      pnodes = nodes + 2 * ii;
      // slice [:,ii:ii+1]
      fmf_pretend( c1, 1, 1, coors->nRow, coors->nCol, coors->val + ii );

      eval_lagrange_simplex( base1d, c1, pnodes, nNod, nCol, order, diff,
			     mtx_i, bc, suppress_errors, eps );

      for (im = 0; im < out->cellSize; im++) {
	out->val[im] *= base1d->val[im];
      }

      ERR_CheckGo( ret );
    }

  } else {

    nr = out->nRow;
    nc = out->nCol;

    for (ii = 0; ii < dim; ii++) {
      // slice [:,2*ii:2*ii+2]
      pnodes = nodes + 2 * ii;
      // slice [:,ii:ii+1]
      fmf_pretend( c1, 1, 1, coors->nRow, coors->nCol, coors->val + ii );

      for (id = 0; id < dim; id++) {
	if (ii == id) {
	  eval_lagrange_simplex( base1d, c1, pnodes, nNod, nCol, order, diff,
				 mtx_i, bc, suppress_errors, eps );
	} else {
	  eval_lagrange_simplex( base1d, c1, pnodes, nNod, nCol, order, 0,
				 mtx_i, bc, suppress_errors, eps );
	}

	// slice [:,id:id+1,:]
	for (im = 0; im < out->nLev; im++) {
	  for (ic = 0; ic < nc; ic++) {
	    out->val[nr*nc*im + nc*id + ic] *= base1d->val[nc*im + ic];
	  }
	}
      }

      ERR_CheckGo( ret );
    }
  }
 
 end_label:

  return( ret );
}


void rezidual( FMField *res, FMField *xi, FMField *coors, FMField *e_coors,
	       FMField *bf, FMField *xint )
{
  int32 iep;

  if (xi->nCol == 3) {
    for (iep = 0; iep < bf->nCol; iep++) {
      bf->val[iep] = baseTriL( xi->val[0], xi->val[1], xi->val[2], iep );
    }
  } else {
    for (iep = 0; iep < bf->nCol; iep++) {
      bf->val[iep] = baseBiL( xi->val[0], xi->val[1], iep );
    }
  }
  // X(xi).
  fmf_mulAB_n1( xint, bf, e_coors );
  // Rezidual.
  fmf_subAB_nn( res, coors, xint );
}

void matrix( FMField *mtx, FMField *xi, FMField *e_coors, FMField *bfg )
{
  int32 nEP = bfg->nCol, dim = bfg->nRow, iep;

  if (dim == 3) {
    for (iep = 0; iep < nEP; iep++) {
      bfg->val[0*nEP+iep] = baseTriLx( xi->val[0], xi->val[1], xi->val[2], iep );
      bfg->val[1*nEP+iep] = baseTriLy( xi->val[0], xi->val[1], xi->val[2], iep );
      bfg->val[2*nEP+iep] = baseTriLz( xi->val[0], xi->val[1], xi->val[2], iep );
    }
  } else {
    for (iep = 0; iep < nEP; iep++) {
      bfg->val[0*nEP+iep] = baseBiLx( xi->val[0], xi->val[1], iep );
      bfg->val[1*nEP+iep] = baseBiLy( xi->val[0], xi->val[1], iep );
    }
  }
  // - Matrix.
  fmf_mulAB_n1( mtx, bfg, e_coors );
}

#undef __FUNC__
#define __FUNC__ "inverse_element_mapping"
int32 inverse_element_mapping( FMField *out,
			       FMField *coors, FMField *e_coors,
			       FMField *ref_coors, int32 i_max, float64 eps )
{
  int32 ii, id, dim, nEP, ret = RET_OK;
  float64 err;
  FMField *bf = 0, *bfg = 0, *res = 0, *mtx = 0, *imtx = 0, *xint = 0;

  dim = ref_coors->nCol;
  nEP = ref_coors->nRow;

  fmf_createAlloc( &bf, 1, 1, 1, nEP );
  fmf_createAlloc( &bfg, 1, 1, dim, nEP );
  fmf_createAlloc( &res, 1, 1, 1, dim );
  fmf_createAlloc( &mtx, 1, 1, dim, dim );
  fmf_createAlloc( &imtx, 1, 1, dim, dim );
  fmf_createAlloc( &xint, 1, 1, 1, dim );

  fmf_fillC( out, 0.0 );
  ii = 0;
  while (ii < i_max) {
    rezidual( res, out, coors, e_coors, bf, xint );
    err = 0.0;
    for (id = 0; id < dim; id++) {
      err += res->val[id] * res->val[id];
    }
    err = sqrt( err );
/*     fmf_print( bf, stdout, 0 ); */
/*     printf( "%d %f\n", ii, err ); */
    if (err < eps) break;

    matrix( mtx, out, e_coors, bfg );
/*     fmf_print( bfg, stdout, 0 ); */
/*     fmf_print( mtx, stdout, 0 ); */
    geme_invert3x3( imtx, mtx );
    ERR_CheckGo( ret );

    fmf_mulAB_nn( xint, res, imtx );
    fmf_addAB_nn( out, out, xint );
    ii += 1;
  }

 end_label:
  fmf_freeDestroy( &bf );
  fmf_freeDestroy( &bfg );
  fmf_freeDestroy( &res );
  fmf_freeDestroy( &mtx );
  fmf_freeDestroy( &imtx );
  fmf_freeDestroy( &xint );

  return( ret );
}
