#include "nurbs.h"

inline void ravel_multi_index(uint32 *index, uint32 *indices,
                              uint32 *shape, uint32 num)
{
  uint32 ii, stride = 1;
  uint32 raveled = 0;

  for (ii = num - 1; ii >= 1; ii--) {
    raveled += stride * indices[ii];
    stride *= shape[ii];
  }
  raveled += stride * indices[0];

  *index = raveled;
}

inline void unravel_index(uint32 *indices, uint32 index,
                          uint32 *shape, uint32 num)
{
  uint32 ii;

  for (ii = num - 1; ii > 0; ii--) {
    indices[ii] = index % shape[ii];
    index /= shape[ii];
  }
  indices[0] = index % shape[0];
}

int32 eval_bernstein_basis(FMField *funs, FMField *ders,
                           float64 x, uint32 degree)
{
  uint32 ip, ifun;
  uint32 n_fun = degree + 1;
  float64 prev, tmp;

  fmf_fillC(funs, 0.0);
  fmf_fillC(ders, 0.0);

  funs->val[0] = 1.0;

  if (degree == 0) {
    return(RET_OK);
  }

  for (ip = 1; ip < n_fun - 1; ip++) {
    prev = 0.0;
    for (ifun = 0; ifun < ip + 1; ifun++) {
      tmp = x * funs->val[ifun];
      funs->val[ifun] = (1.0 - x) * funs->val[ifun] + prev;
      prev = tmp;
    }
  }

  ders->val[0] = degree * (- funs->val[0]);
  for (ifun = 1; ifun < n_fun; ifun++) {
    ders->val[ifun] = degree * (funs->val[ifun - 1] - funs->val[ifun]);
  }

  prev = 0.0;
  for (ifun = 0; ifun < n_fun; ifun++) {
    tmp = x * funs->val[ifun];
    funs->val[ifun] = (1.0 - x) * funs->val[ifun] + prev;
    prev = tmp;
  }

  return(RET_OK);
}

#undef __FUNC__
#define __FUNC__ "eval_nurbs_basis_tp"
/*
  dR_dx has shape (dim, n_efun), transposed w.r.t. the Python version!
*/
int32 eval_nurbs_basis_tp(FMField *R, FMField *dR_dx, FMField *det,
                          FMField *dR_dxi,
                          FMField *dx_dxi, FMField *dxi_dx,
                          FMField *B, FMField *dB_dxi,
                          FMField *N, FMField *dN_dxi,
                          FMField *qp, uint32 ie, FMField *control_points,
                          FMField *weights, int32 *degrees, int32 dim,
                          FMField *cs,
                          int32 *conn, int32 n_el, int32 n_ep)
{
  int32 ret = RET_OK;
  uint32 ii, jj, a, i0, i1, i2;
  uint32 n_efuns[3];
  uint32 n_efun = 1;
  uint32 n_els[3];
  uint32 ic[3];
  int32 *ec;
  FMField *C, *N0, *N1, *N2, *dN0_dxi, *dN1_dxi, *dN2_dxi;
  float64 w, W, P;
  float64 dw_dxi[3];

#ifdef DEBUG_FMF
  if (!((dim == qp->nCol) && (dim <= 3))) {
    errput(ErrHead "inconsistent dimension! (%d == $d <= 3)", dim, qp->nCol);
    ERR_CheckGo(ret);
  }
#endif

  for (ii = 0; ii < dim; ii++) {
    n_efuns[ii] = degrees[ii] + 1;
    n_efun *= n_efuns[ii];
  }

#ifdef DEBUG_FMF
  if (n_efun != n_ep) {
    errput(ErrHead "inconsistent number of bases! (%d == $d)", n_efun, n_ep);
    ERR_CheckGo(ret);
  }
#endif

  // Element connectivity.
  ec = conn + n_ep * ie;

  // 1D Bernstein basis B, dB/dxi.
  for (ii = 0; ii < dim; ii++) {
    eval_bernstein_basis(B + ii, dB_dxi + ii, qp->val[ii], degrees[ii]);
  }

  // 1D B-spline basis N = CB, dN/dxi = C dB/dxi.
  for (ii = 0; ii < dim; ii++) {
    n_els[ii] = (cs + ii)->nCell;
  }

  unravel_index(ic, ie, n_els, dim);

  for (ii = 0; ii < dim; ii++) {
    C = cs + ii;
    FMF_SetCell(C, ic[ii]);

    fmf_mulAB_nn(N + ii, C, B + ii);
    fmf_mulAB_nn(dN_dxi + ii, C, dB_dxi + ii);
  }

  ERR_CheckGo(ret);

  // Numerators and denominator for tensor-product NURBS basis R, dR/dxi.
  w = 0; // w_b
  for (ii = 0; ii < dim; ii++) {
    dw_dxi[ii] = 0.0; // dw_b/dxi
  }
  a = 0; // Basis function index.
  if (dim == 3) {
    N0 = N + 0;
    N1 = N + 1;
    N2 = N + 2;
    dN0_dxi = dN_dxi + 0;
    dN1_dxi = dN_dxi + 1;
    dN2_dxi = dN_dxi + 2;
    for (i0 = 0; i0 < n_efuns[0]; i0++) {
      for (i1 = 0; i1 < n_efuns[1]; i1++) {
        for (i2 = 0; i2 < n_efuns[2]; i2++) {
          W = weights->val[ec[a]];

          R->val[a] = N0->val[i0] * N1->val[i1] * N2->val[i2] * W;
          w += R->val[a];

          dR_dxi->val[a+n_ep*0] = dN0_dxi->val[i0] * N1->val[i1] * N2->val[i2] * W;
          dw_dxi[0] += dR_dxi->val[a+n_ep*0];

          dR_dxi->val[a+n_ep*1] = N0->val[i0] * dN1_dxi->val[i1] * N2->val[i2] * W;
          dw_dxi[1] += dR_dxi->val[a+n_ep*1];

          dR_dxi->val[a+n_ep*2] = N0->val[i0] * N1->val[i1] * dN2_dxi->val[i2] * W;
          dw_dxi[2] += dR_dxi->val[a+n_ep*2];

          a += 1;
        }
      }
    }
  } else if (dim == 2) {
    N0 = N + 0;
    N1 = N + 1;
    dN0_dxi = dN_dxi + 0;
    dN1_dxi = dN_dxi + 1;
    for (i0 = 0; i0 < n_efuns[0]; i0++) {
      for (i1 = 0; i1 < n_efuns[1]; i1++) {
        W = weights->val[ec[a]];

        R->val[a] = N0->val[i0] * N1->val[i1] * W;
        w += R->val[a];

        dR_dxi->val[a+n_ep*0] = dN0_dxi->val[i0] * N1->val[i1] * W;
        dw_dxi[0] += dR_dxi->val[a+n_ep*0];

        dR_dxi->val[a+n_ep*1] = N0->val[i0] * dN1_dxi->val[i1] * W;
        dw_dxi[1] += dR_dxi->val[a+n_ep*1];

        a += 1;
      }
    }
  } else {
    N0 = N + 0;
    dN0_dxi = dN_dxi + 0;
    for (i0 = 0; i0 < n_efuns[0]; i0++) {
        W = weights->val[ec[a]];

        R->val[a] = N0->val[i0] * W;
        w += R->val[a];

        dR_dxi->val[a+n_ep*0] = dN0_dxi->val[i0] * W;
        dw_dxi[0] += dR_dxi->val[a+n_ep*0];

        a += 1;
    }
  }

  // Finish R <- R / w_b.
  fmf_mulC(R, 1.0 / w);

  // Finish dR/dxi. D == W C dB/dxi, dR/dxi = (D - R dw_b/dxi) / w_b.
  for (a = 0; a < dR_dxi->nCol; a++) {
    for (ii = 0; ii < dim; ii++) {
      dR_dxi->val[a+n_ep*ii] = (dR_dxi->val[a+n_ep*ii]
                                - R->val[a] * dw_dxi[ii]) / w;
    }
  }

  // Mapping reference -> physical domain dxi/dx.
  // x = sum P_a R_a, dx/dxi = sum P_a dR_a/dxi, invert.
  for (ii = 0; ii < dim; ii++) {
    for (jj = 0; jj < dim; jj++) {
      dx_dxi->val[dim*ii+jj] = 0.0;
      for (a = 0; a < dR_dxi->nCol; a++) {
        P = control_points->val[dim*ec[a]+ii];

        dx_dxi->val[dim*ii+jj] += P * dR_dxi->val[a+n_ep*jj];
      }
    }
  }
  geme_det3x3(det->val, dx_dxi);

  geme_invert3x3(dxi_dx, dx_dxi);

  // dR/dx.
  fmf_mulATB_nn(dR_dx, dxi_dx, dR_dxi);

 end_label:
  return(ret);
}
