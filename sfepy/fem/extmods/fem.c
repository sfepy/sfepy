#include "fem.h"
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
  int32 ii, iel, ir, irg, stride;
  int32 *pconn;
  float64 aux_r, aux_i;
  float64 *val_r, *val_i;

  stride = vec_r->offset;
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
      aux_r = (sign_r * vecInEls_r->val[stride*ir])
	- (sign_i * vecInEls_i->val[stride*ir]);
      aux_i = (sign_r * vecInEls_i->val[stride*ir])
	+ (sign_i * vecInEls_r->val[stride*ir]);
      val_r[stride*irg] += aux_r - aux_i;
      val_i[stride*irg] += aux_r + aux_i;
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
  int32 ii, iel, ir, ic, irg, icg, is, iloc, found, stride;
  int32 *pconnR, *pconnC;
  float64 aux_r, aux_i;
  float64 *val_r, *val_i;

  stride = mtx_r->offset;
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
	iloc = stride * (nEPC * ir + ic);

/* 	output( "%d %d %d %d %d %d\n", iel, ir, ic, irg, icg, iloc ); */
	/* Try bsearch instead here... */
	found = 0;
	for (is = prows[irg]; is < prows[irg+1]; is++) {
	  if (cols[is] == icg) {
	    aux_r = (sign_r * mtxInEls_r->val[iloc])
	      - (sign_i * mtxInEls_i->val[iloc]);
	    aux_i = (sign_r * mtxInEls_i->val[iloc])
	      + (sign_i * mtxInEls_r->val[iloc]);
	    val_r[stride*is] += aux_r;
	    val_i[stride*is] += aux_i;
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
