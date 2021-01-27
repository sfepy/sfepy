#include "geommech.h"

#define Pi 3.14159265358979

#undef __FUNC__
#define __FUNC__ "geme_invert3x3"
/*!
  @par Revision history:
  - 09.02.2001, c
  - 22.02.2001
  - 06.03.2003, adopted from rcfem2
  - 21.08.2003
*/
int32 geme_invert3x3( FMField *mtxI, FMField *mtx )
{
  int32 il, dim;
  int32 tmp1;
  float64 eps = 1e-55;
  float64 idet, det, *jI, *j;

  dim = mtx->nRow;
  for (il = 0; il < mtx->nLev; il++) {
    tmp1 = dim*dim*il;
    jI = mtxI->val + tmp1;
    j = mtx->val + tmp1;
    switch (dim) {
    case 1:
      idet = j[0];
      if (fabs(idet) < eps) {
        det = 0.0;
      } else {
        det = 1.0 / idet;
      }
      jI[0] = 1.0 * det;
      break;
    case 2:
      idet = j[0] * j[3] - j[1] * j[2];
      if (fabs(idet) < eps) {
        det = 0.0;
      } else {
        det = 1.0 / idet;
      }
      jI[0] = j[3] * det;
      jI[1] = -j[1] * det;
      jI[2] = -j[2] * det;
      jI[3] = j[0] * det;
      break;
    case 3:
      jI[0] = (j[4]*j[8] - j[7]*j[5]);
      jI[1] = -(j[1]*j[8] - j[2]*j[7]);
      jI[2] = (j[1]*j[5] - j[2]*j[4]);
      jI[3] = -(j[3]*j[8] - j[5]*j[6]);
      jI[4] = (j[0]*j[8] - j[2]*j[6]);
      jI[5] = -(j[0]*j[5] - j[2]*j[3]);
      jI[6] = (j[3]*j[7] - j[4]*j[6]);
      jI[7] = -(j[0]*j[7] - j[1]*j[6]);
      jI[8] = (j[0]*j[4] - j[1]*j[3]);
      idet = j[0] * jI[0] + j[1] * jI[3] + j[2] * jI[6];
      if (fabs(idet) < eps) {
        det = 0.0;
      } else {
        det = 1.0 / idet;
      }
      jI[0] *= det; jI[1] *= det; jI[2] *= det;
      jI[3] *= det; jI[4] *= det; jI[5] *= det;
      jI[6] *= det; jI[7] *= det; jI[8] *= det;
      break;
    default:
      errput( ErrHead "ERR_Switch\n" );
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "geme_invert4x4"
int32 geme_invert4x4( FMField *mtxI, FMField *mtx )
{
  int32 ii, il;
  float64 det, buf[16];
  float64 *pm, *pi;

  for (il = 0; il < mtx->nLev; il++) {
    pm = FMF_PtrLevel( mtx, il );
    pi = FMF_PtrLevel( mtxI, il );

    buf[0] = pm[5]*pm[10]*pm[15] - pm[5]*pm[11]*pm[14] - pm[9]*pm[6]*pm[15]
      + pm[9]*pm[7]*pm[14] + pm[13]*pm[6]*pm[11] - pm[13]*pm[7]*pm[10];
    buf[4] = -pm[4]*pm[10]*pm[15] + pm[4]*pm[11]*pm[14] + pm[8]*pm[6]*pm[15]
      - pm[8]*pm[7]*pm[14] - pm[12]*pm[6]*pm[11] + pm[12]*pm[7]*pm[10];
    buf[8] = pm[4]*pm[9]*pm[15] - pm[4]*pm[11]*pm[13] - pm[8]*pm[5]*pm[15]
      + pm[8]*pm[7]*pm[13] + pm[12]*pm[5]*pm[11] - pm[12]*pm[7]*pm[9];
    buf[12] = -pm[4]*pm[9]*pm[14] + pm[4]*pm[10]*pm[13] + pm[8]*pm[5]*pm[14]
      - pm[8]*pm[6]*pm[13] - pm[12]*pm[5]*pm[10] + pm[12]*pm[6]*pm[9];
    buf[1] = -pm[1]*pm[10]*pm[15] + pm[1]*pm[11]*pm[14] + pm[9]*pm[2]*pm[15]
      - pm[9]*pm[3]*pm[14] - pm[13]*pm[2]*pm[11] + pm[13]*pm[3]*pm[10];
    buf[5] = pm[0]*pm[10]*pm[15] - pm[0]*pm[11]*pm[14] - pm[8]*pm[2]*pm[15]
      + pm[8]*pm[3]*pm[14] + pm[12]*pm[2]*pm[11] - pm[12]*pm[3]*pm[10];
    buf[9] = -pm[0]*pm[9]*pm[15] + pm[0]*pm[11]*pm[13] + pm[8]*pm[1]*pm[15]
      - pm[8]*pm[3]*pm[13] - pm[12]*pm[1]*pm[11] + pm[12]*pm[3]*pm[9];
    buf[13] = pm[0]*pm[9]*pm[14] - pm[0]*pm[10]*pm[13] - pm[8]*pm[1]*pm[14]
      + pm[8]*pm[2]*pm[13] + pm[12]*pm[1]*pm[10] - pm[12]*pm[2]*pm[9];
    buf[2] = pm[1]*pm[6]*pm[15] - pm[1]*pm[7]*pm[14] - pm[5]*pm[2]*pm[15]
      + pm[5]*pm[3]*pm[14] + pm[13]*pm[2]*pm[7] - pm[13]*pm[3]*pm[6];
    buf[6] = -pm[0]*pm[6]*pm[15] + pm[0]*pm[7]*pm[14] + pm[4]*pm[2]*pm[15]
      - pm[4]*pm[3]*pm[14] - pm[12]*pm[2]*pm[7] + pm[12]*pm[3]*pm[6];
    buf[10] = pm[0]*pm[5]*pm[15] - pm[0]*pm[7]*pm[13] - pm[4]*pm[1]*pm[15]
      + pm[4]*pm[3]*pm[13] + pm[12]*pm[1]*pm[7] - pm[12]*pm[3]*pm[5];
    buf[14] = -pm[0]*pm[5]*pm[14] + pm[0]*pm[6]*pm[13] + pm[4]*pm[1]*pm[14]
      - pm[4]*pm[2]*pm[13] - pm[12]*pm[1]*pm[6] + pm[12]*pm[2]*pm[5];
    buf[3] = -pm[1]*pm[6]*pm[11] + pm[1]*pm[7]*pm[10] + pm[5]*pm[2]*pm[11]
      - pm[5]*pm[3]*pm[10] - pm[9]*pm[2]*pm[7] + pm[9]*pm[3]*pm[6];
    buf[7] = pm[0]*pm[6]*pm[11] - pm[0]*pm[7]*pm[10] - pm[4]*pm[2]*pm[11]
      + pm[4]*pm[3]*pm[10] + pm[8]*pm[2]*pm[7] - pm[8]*pm[3]*pm[6];
    buf[11] = -pm[0]*pm[5]*pm[11] + pm[0]*pm[7]*pm[9] + pm[4]*pm[1]*pm[11]
      - pm[4]*pm[3]*pm[9] - pm[8]*pm[1]*pm[7] + pm[8]*pm[3]*pm[5];
    buf[15] = pm[0]*pm[5]*pm[10] - pm[0]*pm[6]*pm[9] - pm[4]*pm[1]*pm[10]
      + pm[4]*pm[2]*pm[9] + pm[8]*pm[1]*pm[6] - pm[8]*pm[2]*pm[5];

    det = pm[0]*buf[0] + pm[1]*buf[4] + pm[2]*buf[8] + pm[3]*buf[12];
    if (fabs(det) < 1e-55) {
      output("possibly singular matrix!\n");
    }
    det = 1.0 / det;

    for (ii = 0; ii < 16; ii++) {
      pi[ii] = buf[ii] * det;
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "geme_tensor2vectorS3"
/*!
  @par Revision history:
  - 13.11.2001, c
  - 06.03.2003, adopted from rcfem2
*/
int32 geme_tensor2vectorS3( FMField *vec, FMField *mtx )
{
  int32 il, dim;
  float64 *pmtx, *pvec;

  dim = mtx->nRow;
  for (il = 0; il < mtx->nLev; il++) {
    pvec = FMF_PtrLevel( vec, il );
    pmtx = FMF_PtrLevel( mtx, il );
    switch (dim) {
    case 1:
      pvec[0] = pmtx[0];
      break;
    case 2:
      pvec[0] = pmtx[0];
      pvec[1] = pmtx[3];
      pvec[2] = pmtx[1];
      break;
    case 3:
      pvec[0] = pmtx[0];
      pvec[1] = pmtx[4];
      pvec[2] = pmtx[8];
      pvec[3] = pmtx[1];
      pvec[4] = pmtx[2];
      pvec[5] = pmtx[5];
      break;
    default:
      errput( ErrHead "ERR_Switch\n" );
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "geme_det3x3"
/*!
  @par Revision history:
  - 09.02.2001, c
  - 06.03.2003, adopted from rcfem2
*/
int32 geme_det3x3( float64 *det, FMField *mtx )
{
  int32 il, dim;
  float64 *j;

  dim = mtx->nRow;
  for (il = 0; il < mtx->nLev; il++) {
    j = mtx->val + dim*dim*il;
    switch (dim) {
    case 1:
      det[il] = j[0];
      break;
    case 2:
      det[il] = j[0] * j[3] - j[1] * j[2];
      break;
    case 3:
      det[il] = j[0]*j[4]*j[8] + j[3]*j[7]*j[2] + j[1]*j[5]*j[6]
	- j[2]*j[4]*j[6] - j[5]*j[7]*j[0] - j[1]*j[3]*j[8];
      break;
    default:
      errput( ErrHead "ERR_Switch\n" );
    }
  }
  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "geme_trace3x3"
/*!
  @par Revision history:
  - 13.11.2001, c
  - 06.03.2003, adopted from rcfem2
*/
int32 geme_trace3x3( float64 *tr, FMField *mtx )
{
  int32 il, dim;
  float64 *j;

  dim = mtx->nRow;
  for (il = 0; il < mtx->nLev; il++) {
    j = mtx->val + dim*dim*il;
    switch (dim) {
    case 1:
      tr[il] = j[0];
      break;
    case 2:
      tr[il] = j[0] + j[3];
      break;
    case 3:
      tr[il] = j[0] +j[4] + j[8];
      break;
    default:
      errput( ErrHead "ERR_Switch\n" );
    }
  }
  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "geme_invar1"
int32 geme_invar1( float64 *invar, FMField *mtx )
{
  int32 il, dim;
  float64 *j;

  dim = mtx->nRow;
  for (il = 0; il < mtx->nLev; il++) {
    j = mtx->val + dim*dim*il;
    switch (dim) {
    case 1:
      invar[il] = j[0];
      break;
    case 2: /* plain strain */
      invar[il] = 1.0 + j[0] + j[3];
      break;
    case 3:
      invar[il] = j[0] +j[4] + j[8];
      break;
    default:
      errput( ErrHead "ERR_Switch\n" );
    }
  }
  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "geme_invar2"
int32 geme_invar2( float64 *invar, FMField *mtx )
{
  int32 il, dim;
  float64 *j;

  dim = mtx->nRow;
  for (il = 0; il < mtx->nLev; il++) {
    j = mtx->val + dim*dim*il;
    switch (dim) {
    case 1: /* no sense in 1D */
      invar[il] = 0.0;
      break;
    case 2: /* plain strain */
      invar[il] = j[0]*j[3] + j[0] + j[3] - j[1]*j[1];
      break;
    case 3:
      invar[il] = j[0]*j[4] + j[0]*j[8] + j[4]*j[8]
	- j[1]*j[1] - j[2]*j[2] - j[5]*j[5];
      break;
    default:
      errput( ErrHead "ERR_Switch\n" );
    }
  }
  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "geme_norm3"
/*!
  @par Revision history:
  - 23.02.2001, c
  - 06.03.2003, adopted from rcfem2
*/
int32 geme_norm3( float64 *out, FMField *mtx )
{
  int32 il, dim;
  float64 *j;

  dim = mtx->nRow;
  for (il = 0; il < mtx->nLev; il++) {
    j = mtx->val + dim*il;
    switch (dim) {
    case 1:
      out[il] = fabs(j[0]);
      break;
    case 2:
      out[il] = sqrt( j[0] * j[0] + j[1] * j[1] );
      break;
    case 3:
      out[il] = sqrt( j[0] * j[0] + j[1] * j[1] + j[2] * j[2] );
      break;
    default:
      errput( ErrHead "ERR_Switch\n" );
    }
  }
  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "geme_eig3x3"
/*!
  @par Revision history:
  - 18.12.2003, c
*/
int32 geme_eig3x3( float64 *out, FMField *mtx )
{
  int32 il, dim;
  float64 *j, *val;

  dim = mtx->nRow;
  for (il = 0; il < mtx->nLev; il++) {
    j = mtx->val + dim*dim*il;
    val = out + dim*il;
    switch (dim) {
    case 1:
      val[0] = j[0];
      break;
    case 2: {
      // See Numerical Recipes.
      float64 b, c, q;

      b = -j[0] - j[2];
      c = j[0] * j[2] - j[1] * j[3];
      q = - 0.5 * (b + Sgn(b) * sqrt( b * b - 4.0 * c ));

      val[0] = q;
      val[1] = c / q;
    } break;
    case 3: {
      // See Numerical Recipes.
      float64 a, b, c, q, r, t;

      a = -(j[0] + j[4] + j[8]);
      b = j[0] * j[4] + j[0] * j[8] + j[4] * j[8]
	- j[3] * j[1] - j[6] * j[2] - j[7] * j[5];
      c = j[0] * j[5] * j[7] + j[4] * j[6] * j[2] + j[8] * j[1] * j[3]
	- j[6] * j[1] * j[5] - j[0] * j[4] * j[8] - j[3] * j[2] * j[7];

      q = (a * a - 3.0 * b) / 9.0;
      r = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c) / 54.0;

      if (((q * q * q) - (r * r)) > MachEps) {
	t = acos( r / sqrt( q * q * q ) );
      } else {
	t = Pi;
      }
/*       output( "%e %e %e\n", r * r, q * q * q, t ); */
      val[0] = -2.0 * sqrt( q ) * cos( (t) / 3.0 ) - a / 3.0;
      val[1] = -2.0 * sqrt( q ) * cos( (t+2.0*Pi) / 3.0 ) - a / 3.0;
      val[2] = -2.0 * sqrt( q ) * cos( (t-2.0*Pi) / 3.0 ) - a / 3.0;
    } break;
    default:
      errput( ErrHead "ERR_Switch\n" );
    }
  }
  return( RET_OK );
}


#undef __FUNC__
#define __FUNC__ "geme_mulAVSB3"
/*!
  R = A B, A is upper triangle stored as vector, for dim <= 3.
  e.g. \sigma n == \sigma^T n, \sigma stored as vector (sym).

  @par Revision history:
  - 06.09.2006, c
*/
int32 geme_mulAVSB3( FMField *out, FMField *vs, FMField *in )
{
  int32 iqp, ir, ic, ii, dim, nQP, nc;
  int32 _is[] = {0, 0, 0, 0, 0, 0, 0, 0, 0,
		 0, 2, 2, 1, 0, 0, 0, 0, 0,
		 0, 3, 4, 3, 1, 5, 4, 5, 2};
  int32 *is;
  float64 *pout, *pvs, *pin;

  nQP = vs->nLev;
  dim = in->nRow;
  nc = out->nCol;

#ifdef DEBUG_FMF
  if ((vs->nRow != (dim * (dim + 1) / 2)) || (out->nLev != in->nLev)
      || (vs->nCol != 1) || (out->nLev != vs->nLev) || (nc != in->nCol)
      || (dim != out->nRow)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d), (%d %d %d), (%d %d %d)\n",
	    out->nLev, out->nRow, out->nCol,
	    vs->nLev, vs->nRow, vs->nCol,
	    in->nLev, in->nRow, in->nCol );
  }
#endif

  is = _is + 9 * (dim - 1);

  for (iqp = 0; iqp < nQP; iqp++) {
    pvs = FMF_PtrLevel( vs, iqp );
    pout = FMF_PtrLevel( out, iqp );
    pin = FMF_PtrLevel( in, iqp );

    for (ir = 0; ir < dim; ir++) {
      for (ic = 0; ic < nc; ic++) {
	pout[nc*ir+ic] = 0.0;
	for (ii = 0; ii < dim; ii++) {
	  pout[nc*ir+ic] += pvs[is[dim*ir+ii]] * pin[nc*ii+ic];
	}
      }
    }
  }

  return( RET_OK );
}

/*!
  @par Revision history:
  - 14.11.2001, c
  - 15.11.2001
*/
int32 t2i3D[] = {0, 1, 2, 0, 0, 1};
int32 t2j3D[] = {0, 1, 2, 1, 2, 2};
int32 t4s3D[] = {0, 3, 4,
		 3, 1, 5,
		 4, 5, 2};

int32 t2i2D[] = {0, 1, 0};
int32 t2j2D[] = {0, 1, 1};
int32 t4s2D[] = {0, 2,
		 2, 1};

int32 t2i1D[] = {0};
int32 t2j1D[] = {0};
int32 t4s1D[] = {0};

#undef __FUNC__
#define __FUNC__ "geme_mulT2ST2S_T4S_ikjl"
/*!
  @par Revision history:
  - 14.11.2001, c
  - 06.03.2003, adopted from rcfem2
*/
int32 geme_mulT2ST2S_T4S_ikjl( FMField *t4, FMField *t21, FMField *t22 )
{
  int32 iqp, ir, ic, ii, ij, ik, il, s1, s2;
  int32 sym = t4->nRow, dim;
  float64 *pt4, *pt21, *pt22;
  int32 *t2i = 0, *t2j = 0, *t4s = 0;

  dim = sym2dim( sym );

  switch (dim) {
  case 1:
    t2i = t2i1D;
    t2j = t2j1D;
    t4s = t4s1D;
    break;
  case 2:
    t2i = t2i2D;
    t2j = t2j2D;
    t4s = t4s2D;
    break;
  case 3:
    t2i = t2i3D;
    t2j = t2j3D;
    t4s = t4s3D;
    break;
  default:
    errput( ErrHead "ERR_Switch\n" );
  }

  for (iqp = 0; iqp < t4->nLev; iqp++) {
    pt4 = FMF_PtrLevel( t4, iqp );
    pt21 = FMF_PtrLevel( t21, iqp );
    pt22 = FMF_PtrLevel( t22, iqp );
    for (ir = 0; ir < sym; ir++) {
      for (ic = 0; ic < sym; ic++) {
	ii = t2i[ir];
	ij = t2j[ir];
	ik = t2i[ic];
	il = t2j[ic];
	s1 = t4s[dim*ii+ik];
	s2 = t4s[dim*ij+il];
	pt4[sym*ir+ic] = pt21[s1] * pt22[s2];
      }
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "geme_mulT2ST2S_T4S_iljk"
/*!
  @par Revision history:
  - 14.11.2001, c
  - 06.03.2003, adopted from rcfem2
*/
int32 geme_mulT2ST2S_T4S_iljk( FMField *t4, FMField *t21, FMField *t22 )
{
  int32 iqp, ir, ic, ii, ij, ik, il, s1, s2;
  int32 sym = t4->nRow, dim;
  float64 *pt4, *pt21, *pt22;
  int32 *t2i = 0, *t2j = 0, *t4s = 0;

  dim = sym2dim( sym );

  switch (dim) {
  case 1:
    t2i = t2i1D;
    t2j = t2j1D;
    t4s = t4s1D;
    break;
  case 2:
    t2i = t2i2D;
    t2j = t2j2D;
    t4s = t4s2D;
    break;
  case 3:
    t2i = t2i3D;
    t2j = t2j3D;
    t4s = t4s3D;
    break;
  default:
    errput( ErrHead "ERR_Switch\n" );
  }

  for (iqp = 0; iqp < t4->nLev; iqp++) {
    pt4 = FMF_PtrLevel( t4, iqp );
    pt21 = FMF_PtrLevel( t21, iqp );
    pt22 = FMF_PtrLevel( t22, iqp );
    for (ir = 0; ir < sym; ir++) {
      for (ic = 0; ic < sym; ic++) {
	ii = t2i[ir];
	ij = t2j[ir];
	ik = t2i[ic];
	il = t2j[ic];
	s1 = t4s[dim*ii+il];
	s2 = t4s[dim*ij+ik];
	pt4[sym*ir+ic] = pt21[s1] * pt22[s2];
      }
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "geme_mulT2S_AA"
/*! @f$ R_{ij} = A_{ik} * A_{kj} */
int32 geme_mulT2S_AA( FMField *R, FMField *A )
{
  int32 il, sym;
  float64 *pr, *pa;

  sym = R->nRow;

  pr = R->val;
  pa = A->val;

  for ( il = 0; il < R->nLev; il++ ) {
    switch( sym ) {
    case 1:
      pr[0] = pa[0]*pa[0];
      break;
    case 3:
      pr[0] = pa[0]*pa[0] + pa[2]*pa[2];
      pr[1] = pa[2]*pa[2] + pa[1]*pa[1];
      pr[2] = pa[2]*pa[0] + pa[1]*pa[2];
      break;
    case 6:
      pr[0] = pa[0]*pa[0] + pa[5]*pa[5] + pa[4]*pa[4];
      pr[1] = pa[5]*pa[5] + pa[1]*pa[1] + pa[3]*pa[3];
      pr[2] = pa[4]*pa[4] + pa[3]*pa[3] + pa[2]*pa[2];
      pr[3] = pa[4]*pa[5] + pa[1]*pa[3] + pa[2]*pa[3];
      pr[4] = pa[0]*pa[4] + pa[5]*pa[3] + pa[4]*pa[2];
      pr[5] = pa[0]*pa[5] + pa[5]*pa[1] + pa[4]*pa[3];
      break;
    default:
      errput( ErrHead "ERR_Switch\n" );
    } /* switch */
    pr += sym;
    pa += sym;
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "geme_elementVolume"
/*!
  @par Revision history:
  - 02.05.2001, c
  - 06.03.2003, adopted from rcfem2
*/
int32 geme_elementVolume( float64 *volume, float64 *jacobian, int32 nQP )
{
  int32 qp;

  *volume = 0.0;
  for (qp = 0; qp < nQP; qp++) {
    *volume += jacobian[qp];
  }

  return( RET_OK );
}

static int32 aux1[] = {0, 1, 2, 0, 0, 1};
static int32 aux2[] = {0, 1, 2, 1, 2, 2};

#undef __FUNC__
#define __FUNC__ "geme_buildOpOmega_VS3"
/*!
  @par Revision history:
  - 18.10.2001, c
  - 06.03.2003, adopted from rcfem2
*/
int32 geme_buildOpOmega_VS3( float64 *pomega, float64 *pdir,
			     int32 nItem, int32 dim, int32 sym )
{
  int32 ii, ir;
  for (ii = 0; ii < nItem; ii++) {
    for (ir = 0; ir < sym; ir++) {
      pomega[ir] = pdir[aux1[ir]] * pdir[aux2[ir]];
/*        printf( "%d %d, %f %f %f\n", ii, ir, */
/*  	      pomega[ir], pdir[aux1[ir]], pdir[aux2[ir]] ); */
    }
    pomega += sym;
    pdir += dim;
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "geme_projectToDir"
/*!
  @par Revision history:
  - 18.10.2001, c
  - 06.03.2003, adopted from rcfem2
*/
int32 geme_projectToDir( float64 *pdef, float64 *pomega,
			 float64 *pstrain, int32 nItem, int32 size )
{
  int32 ii, ir;
  for (ii = 0; ii < nItem; ii++) {
    pdef[ii] = 0.0;
    for (ir = 0; ir < size; ir++) {
      pdef[ii] += pomega[ir] * pstrain[ir];
/*        printf( "%d %d, %f %f %f\n", ii, ir, */
/*  	      pdef[ii], pomega[ir], pstrain[ir] ); */
    }
    pomega += size;
    pstrain += size;
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "bf_act"
/*!
  @a in is a state vector stored as matrix (dim, nEP) - in->nLev == 1.

  @par Revision history:
  - 12.12.2005, c
  - 14.12.2005
*/
int32 bf_act( FMField *out, FMField *bf, FMField *in )
{
  int32 iqp, ir, ic, nEP, nQP, dim;
  float64 *pout, *pbf;

  nEP = bf->nCol;
  nQP = bf->nLev;
  dim = in->nRow;

#ifdef DEBUG_FMF
  if ((out->nRow != dim) || (out->nCol != 1)
      || (out->nLev != bf->nLev) || (in->nCol != nEP)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d), (%d %d %d), (%d %d %d)\n",
	    out->nLev, out->nRow, out->nCol,
	    bf->nLev, bf->nRow, bf->nCol,
	    in->nLev, in->nRow, in->nCol );
  }
#endif

  fmf_fillC( out, 0.0 );
  for (iqp = 0; iqp < nQP; iqp++) {
    pbf = FMF_PtrLevel( bf, iqp );
    pout = FMF_PtrLevel( out, iqp );

    for (ic = 0; ic < dim; ic++ ) {
      for (ir = 0; ir < nEP; ir++) {
	pout[ic] += pbf[ir] * in->val[nEP*ic+ir];
      }
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "bf_ract"
/*!
  Act from right.

  @par Revision history:
  - 12.12.2005, c
*/
int32 bf_ract( FMField *out, FMField *bf, FMField *in )
{
  int32 iqp, ir, ic, ii, nEP, nQP, dim;
  float64 *pout, *pbf, *pin;

  nEP = bf->nCol;
  nQP = bf->nLev;
  dim = in->nCol;

#ifdef DEBUG_FMF
  if ((out->nRow != in->nRow) || (out->nCol != (dim * nEP))
      || (out->nLev != bf->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d), (%d %d %d), (%d %d %d)\n",
	    out->nLev, out->nRow, out->nCol,
	    bf->nLev, bf->nRow, bf->nCol,
	    in->nLev, in->nRow, in->nCol );
  }
#endif

  fmf_fillC( out, 0.0 );
  for (iqp = 0; iqp < nQP; iqp++) {
    pbf = FMF_PtrLevel( bf, iqp );
    pout = FMF_PtrLevel( out, iqp );
    pin = FMF_PtrLevel( in, iqp );

    for (ir = 0; ir < out->nRow; ir++) {
      for (ic = 0; ic < dim; ic++) {
	for (ii = 0; ii < nEP; ii++) {
	  (*pout) = (*pin) * pbf[ii];
	  pout++;
	}
	pin++;
      }
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "bf_actt"
/*!
  Act transposed.

  @par Revision history:
  - 12.12.2005, c
  - 20.12.2005
*/
int32 bf_actt( FMField *out, FMField *bf, FMField *in )
{
  int32 iqp, ir, ic, ii, nEP, nQP, dim;
  float64 *pout, *pbf, *pin;

  nEP = bf->nCol;
  nQP = bf->nLev;
  dim = in->nRow;

#ifdef DEBUG_FMF
  if ((out->nCol != in->nCol) || (out->nRow != (dim * nEP))
      || (out->nLev != bf->nLev)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d), (%d %d %d), (%d %d %d)\n",
	    out->nLev, out->nRow, out->nCol,
	    bf->nLev, bf->nRow, bf->nCol,
	    in->nLev, in->nRow, in->nCol );
  }
#endif

  fmf_fillC( out, 0.0 );
  for (iqp = 0; iqp < nQP; iqp++) {
    pbf = FMF_PtrLevel( bf, iqp );
    pout = FMF_PtrLevel( out, iqp );
    pin = FMF_PtrLevel( in, iqp );

    for (ir = 0; ir < dim; ir++) {
      for (ic = 0; ic < out->nCol; ic++) {
	for (ii = 0; ii < nEP; ii++) {
	  pout[out->nCol*ii+ic] = pbf[ii] * (*pin);
	}
	pin++;
      }
      pout += out->nCol * nEP;
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "bf_actt_c1"
/*!
  @a in is a vector (dim, 1).

  @par Revision history:
  - 06.09.2006, c
*/
int32 bf_actt_c1( FMField *out, FMField *bf, FMField *in )
{
  int32 iqp, ir, ic, nEP, nQP, dim;
  float64 *pout, *pbf, *pin;

  nEP = bf->nCol;
  nQP = bf->nLev;
  dim = in->nRow;

#ifdef DEBUG_FMF
  if ((out->nRow != dim * nEP) || (out->nCol != 1) || (out->nLev != in->nLev)
      || (out->nLev != bf->nLev) || (in->nRow != dim) || (in->nCol != 1)) {
    errput( ErrHead "ERR_BadMatch: (%d %d %d), (%d %d %d), (%d %d %d)\n",
	    out->nLev, out->nRow, out->nCol,
	    bf->nLev, bf->nRow, bf->nCol,
	    in->nLev, in->nRow, in->nCol );
  }
#endif

  for (iqp = 0; iqp < nQP; iqp++) {
    pbf = FMF_PtrLevel( bf, iqp );
    pout = FMF_PtrLevel( out, iqp );
    pin = FMF_PtrLevel( in, iqp );

    for (ic = 0; ic < dim; ic++ ) {
      for (ir = 0; ir < nEP; ir++) {
	pout[nEP*ic+ir] = pbf[ir] * pin[ic];
      }
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "bf_buildFTF"
/*!
  @par Revision history:
  - 21.11.2006, c
*/
int32 bf_buildFTF( FMField *ftf, FMField *ftf1 )
{
  int32 iqp, ir, ic, nEPR, nEPC, nQP, dim;
  float64 *pftf, *pftf1;
  float64 val;

  fmf_fillC( ftf, 0.0 );

  nEPR = ftf1->nRow;
  nEPC = ftf1->nCol;
  nQP = ftf1->nLev;
  dim = ftf->nRow / nEPR;

  for (iqp = 0; iqp < nQP; iqp++) {
    pftf1 = FMF_PtrLevel( ftf1, iqp );
    pftf = FMF_PtrLevel( ftf, iqp );
    for (ir = 0; ir < nEPR; ir++) {
      for (ic = 0; ic < nEPC; ic++) {
	val = pftf1[nEPC*ir+ic];

	pftf[dim*nEPC*ir+ic] = val;
	if (dim == 1) continue;
	pftf[dim*nEPC*(nEPR+ir)+ic+nEPC] = val;
	if (dim == 2) continue;
	pftf[dim*nEPC*(2*nEPR+ir)+ic+2*nEPC] = val;
      }
    }
  }

  return( RET_OK );
}


/*!
  @par Revision history:
  - 04.08.2006, c
  - 24.04.2007
*/
void debug_printConn( int32 *conn, int32 nEP )
{
  int32 ii;

  for (ii = 0; ii < nEP; ii++) {
    printf( FI32" ", conn[ii] );
  }
  printf( "\n" );
}

#undef __FUNC__
#define __FUNC__ "ele_extractNodalValuesNBN"
/*!
  Extract values from element nodes, for example coordinates of nodes.
  Each node must have @ dim values!

  Input vector order is node-by-node.
  The extraction order is node-by-node.

  Use when local shape is (nEP, dpn).

  @par Revision history:
  - 26.10.2005, adapted from mafest1
*/
int32 ele_extractNodalValuesNBN( FMField *out, FMField *in,
				 int32 *conn )
{
  int32 inod, idof;

  for (inod = 0; inod < out->nRow; inod++) {
    for (idof = 0; idof < out->nCol; idof++ ) {
      out->val[out->nCol*inod+idof] = in->val[out->nCol*conn[inod]+idof];
    }
  }

  return( RET_OK );
}

#undef __FUNC__
#define __FUNC__ "ele_extractNodalValuesDBD"
/*!
  Extract values from element nodes, for example coordinates of nodes.
  Each node must have @ dim values!

  Input vector order is node-by-node.
  The extraction order is dof-by-dof.

  Use when local shape is (dpn, nEP).

  @par Revision history:
  - 14.12.2005, c
*/
int32 ele_extractNodalValuesDBD( FMField *out, FMField *in,
				 int32 *conn )
{
  int32 inod, idof;

  for (idof = 0; idof < out->nRow; idof++ ) {
    for (inod = 0; inod < out->nCol; inod++) {
      out->val[out->nCol*idof+inod] = in->val[out->nRow*conn[inod]+idof];
    }
  }

  return( RET_OK );
}
