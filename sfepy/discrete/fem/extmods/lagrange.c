#include "geommech.h"

#include "lagrange.h"

void print_context_lagrange(void *_ctx)
{
  LagrangeContext *ctx = (LagrangeContext *) _ctx;

  int32 ir, ic;

  output("order: %d\n", ctx->order);
  output("is_bubble: %d\n", ctx->is_bubble);
  output("tdim: %d\n", ctx->tdim);
  output("nodes:\n");
  for (ir = 0; ir < ctx->n_nod; ir ++) {
    for (ic = 0; ic < ctx->n_col; ic ++) {
      output(" %d", ctx->nodes[ctx->n_col*ir+ic]);
    }
    output("\n");
  }
  output("n_nod: %d\n", ctx->n_nod);
  output("n_col: %d\n", ctx->n_col);

  output("ref_coors:\n");
  fmf_print(ctx->ref_coors, stdout, 0);
  output("vmin: %.4e\n", ctx->vmin);
  output("vmax: %.4e\n", ctx->vmax);

  output("mtx_i:\n");
  fmf_print(ctx->mtx_i, stdout, 0);
  output("bc:\n");
  fmf_print(ctx->bc, stdout, 0);
  output("base1d:\n");
  fmf_print(ctx->base1d, stdout, 0);

  output("eps: %.4e\n", ctx->eps);
  output("check_errors: %d\n", ctx->check_errors);
  output("i_max: %d\n", ctx->i_max);
  output("newton_eps: %.4e\n", ctx->newton_eps);
}

/*
  Get barycentric (area in 2D, volume in 3D) coordinates of points.

  Parameters
  ----------
  bc : FMField
      The barycentric coordinates, shape `(n_coor, dim + 1)`. Then
      reference element coordinates `xi = dot(bc, ref_coors)`.
  coors : FMField
      The coordinates of the points, shape `(n_coor, dim)`.
  ctx : LagrangeContext
      The context data.

  Notes
  -----
  Uses the following items of LagrangeContext:

  mtx_i : FMField
      The inverse of simplex coordinates matrix, shape `(dim + 1, dim + 1)`.
  eps : float
      The tolerance for snapping out-of-simplex point back to the simplex.
  check_errors : bool
      If True, raise ValueError if a barycentric coordinate is outside
      the snap interval `[-eps, 1 + eps]`.
*/
int32 get_barycentric_coors(FMField *bc, FMField *coors, void *_ctx)
{
  LagrangeContext *ctx = (LagrangeContext *) _ctx;
  FMField *mtx_i = ctx->mtx_i;
  float64 eps = ctx->eps;
  int32 check_errors = ctx->check_errors;

  int32 ii, ir, ic, error, ret = RET_OK;
  int32 n_coor = coors->nRow;
  int32 nc = coors->nCol;
  int32 n_v = mtx_i->nRow;
  int32 dim = n_v - 1;
  float64 val;

  for (ir = 0; ir < n_coor; ir++) {
    for (ic = 0; ic < n_v; ic++) {
      val = 0.0;
      for (ii = 0; ii < dim; ii++) {
        val += mtx_i->val[n_v*ic+ii] * coors->val[nc*ir+ii];
      }
      val += mtx_i->val[n_v*ic+dim];
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

      if ((check_errors) && (error)) {
        errput("quadrature point %d outside of element! (%.e)\n", ir, val);
        errset("quadrature point outside of element (see above)!");
      }
      bc->val[n_v*ir+ic] = val;

      ERR_CheckGo(ret);
    }
  }

 end_label:
  return(ret);
}

#undef __FUNC__
#define __FUNC__ "get_xi_dist"
// Returns: ok = 1: success, ok = 0: failure (only for tensor product cells).
int32 get_xi_dist(float64 *pdist, FMField *xi,
                  FMField *point, FMField *e_coors,
                  void *_ctx)
{
  LagrangeContext *ctx = (LagrangeContext *) _ctx;
  FMField *bc = ctx->bc;
  int32 n_v = e_coors->nRow;
  int32 dim = e_coors->nCol;
  int32 tdim = ctx->tdim;
  float64 vmin = ctx->vmin;
  float64 vmax = ctx->vmax;

  int32 ii, ok;
  float64 dist = 0.0, val;
  float64 buf4[4];
  float64 buf8[8];

  if (n_v == (dim + 1)) {
    fmf_pretend_nc(bc, 1, 1, 1, tdim + 1, buf4);
    get_xi_simplex(xi, point, e_coors, _ctx);

    // dist == 0 for 0 <= bc <= 1.
    for (ii = 0; ii < n_v; ii++) {
      val = Min(Max(bc->val[ii] - 1.0, 0.0), 100.0);
      dist += val * val;
      val = Min(Max(0.0 - bc->val[ii], 0.0), 100.0);
      dist += val * val;
    }
    ok = 1;

  } else {
    fmf_pretend_nc(ctx->base1d, 1, 1, 1, n_v, buf8);

    ok = get_xi_tensor(xi, point, e_coors, _ctx);

    // dist == 0 for vmin <= xi <= vmax and ok == 0.
    if (ok == 0) {
      ok = 1;
      for (ii = 0; ii < dim; ii++) {
        val = Min(Max(xi->val[ii] - vmax, 0.0), 100.0);
        dist += val * val;
        val = Min(Max(vmin - xi->val[ii], 0.0), 100.0);
        dist += val * val;
      }
    } else {
      ok = 0;
      dist = 1e10;
    }
  }
  *pdist = dist;

  return(ok);
}

/*
  Get reference simplex coordinates `xi` of `dest_point` given spatial
  element coordinates `e_coors` and coordinates of reference simplex
  vertices `ref_coors` (in `ctx`). Output also the corresponding barycentric
  coordinates `bc` (in `ctx`).
*/
int32 get_xi_simplex(FMField *xi, FMField *dest_point, FMField *e_coors,
                     void *_ctx)
{
  LagrangeContext *ctx = (LagrangeContext *) _ctx;

  int32 idim, ii;
  int32 n_v = e_coors->nRow;
  int32 dim = e_coors->nCol;
  FMField mtx[1], mtx_i[1], rhs[1];
  float64 buf16[16], buf16_2[16], buf4[4];

  fmf_pretend_nc(mtx, 1, 1, n_v, n_v, buf16);
  fmf_pretend_nc(mtx_i, 1, 1, n_v, n_v, buf16_2);
  fmf_pretend_nc(rhs, 1, 1, 1, n_v, buf4);

  for (idim = 0; idim < dim; idim++) {
    for (ii = 0; ii < n_v; ii++) {
      mtx->val[n_v*idim+ii] = e_coors->val[dim*ii+idim];
      rhs->val[idim] = dest_point->val[idim];
    }
  }

  for (ii = 0; ii < n_v; ii++) {
    mtx->val[n_v*dim+ii] = 1.0;
    rhs->val[dim] = 1.0;
  }
  if (dim == 3) {
    geme_invert4x4(mtx_i, mtx);
  } else {
    geme_invert3x3(mtx_i, mtx);
  }

  fmf_mulABT_nn(ctx->bc, rhs, mtx_i);
  fmf_mulAB_nn(xi, ctx->bc, ctx->ref_coors);

  return(RET_OK);
}

/*
  Get reference tensor product element coordinates using Newton method.

  Uses linear 1D base functions.
*/
int32 get_xi_tensor(FMField *xi, FMField *dest_point, FMField *e_coors,
                    void *_ctx)
{
  LagrangeContext *ctx = (LagrangeContext *) _ctx;
  FMField *bc = ctx->bc;

  float64 vmin = ctx->vmin;
  float64 vmax = ctx->vmax;
  int32 i_max = ctx->i_max;
  float64 newton_eps = ctx->newton_eps;

  int32 idim, ii, ok = 0;
  int32 dim = e_coors->nCol;
  int32 powd = 1 << dim;
  FMField bf[1], bfg[1], xint[1], res[1], mtx[1], imtx[1];
  float64 err;
  float64 buf6[6], buf8[8], buf24[24];
  float64 buf3_1[3], buf3_2[3], buf9_1[9], buf9_2[9];

  fmf_pretend_nc(bc, dim, 1, 1, 2, buf6);
  fmf_pretend_nc(bf, 1, 1, 1, powd, buf8);
  fmf_pretend_nc(bfg, 1, 1, dim, powd, buf24);
  fmf_pretend_nc(res, 1, 1, 1, dim, buf3_1);
  fmf_pretend_nc(xint, 1, 1, 1, dim, buf3_2);
  fmf_pretend_nc(mtx, 1, 1, dim, dim, buf9_1);
  fmf_pretend_nc(imtx, 1, 1, dim, dim, buf9_2);

  ii = 0;
  fmf_fillC(xi, 0.5 * (vmin + vmax));
  while (ii < i_max) {
    // Base(xi).
    for (idim = 0; idim < dim; idim++) {
      FMF_SetCell(bc, idim);
      // slice [:,idim:idim+1]
      bc->val[1] = (xi->val[idim] - vmin) / (vmax - vmin);
      bc->val[0] = 1.0 - bc->val[1];
    }

    eval_lagrange_tensor_product(bf, 1, 0, _ctx);

    // X(xi).
    fmf_mulAB_n1(xint, bf, e_coors);
    // Rezidual.
    fmf_subAB_nn(res, dest_point, xint);

    err = 0.0;
    for (idim = 0; idim < dim; idim++) {
      err += res->val[idim] * res->val[idim];
    }
    err = sqrt(err);

    if (err < newton_eps) {
      break;
    }

    // grad Base(xi).
    eval_lagrange_tensor_product(bfg, 1, 1, _ctx);
    // - Matrix.
    fmf_mulAB_n1(mtx, bfg, e_coors);

    geme_invert3x3(imtx, mtx);

    fmf_mulAB_nn(xint, res, imtx);
    fmf_addAB_nn(xi, xi, xint);
    ii += 1;
  }

  if (ii == i_max) {
    ok = 2;
    // Newton did not converge.
    // Use centre if nothing better is available.
    fmf_fillC(xi, 0.5 * (vmin + vmax));
  }

  return(ok);
}

/*
  Evaluate Lagrange base polynomials in given points on simplex domain.
*/
int32 eval_lagrange_simplex(FMField *out, int32 order, int32 diff,
                            void *_ctx)
{
  LagrangeContext *ctx = (LagrangeContext *) _ctx;
  FMField *bc = ctx->bc;
  FMField *mtx_i = ctx->mtx_i;
  int32 *nodes = ctx->nodes;
  int32 n_col = ctx->n_col;

  int32 ret = RET_OK;
  int32 n_coor = bc->nRow;
  int32 n_v = bc->nCol;
  int32 dim = n_v - 1;
  int32 n_nod = out->nCol;
  int32 ii, ir, ic, i1, i2, inod, n_i1, n_ii;
  float64 dval, dd, vv, bci1, bcii;
  float64 *pout;

  if (n_coor != out->nLev) {
    errset("coordinates size mismatch!");
    ERR_CheckGo(ret);
  }

  if (!diff) {
    for (ic = 0; ic < n_coor; ic++) {
      pout = FMF_PtrLevel(out, ic);

      for (inod = 0; inod < n_nod; inod++) {
        pout[inod] = 1.0;

        for (i1 = 0; i1 < n_v; i1++) {
          n_i1 = nodes[n_col*inod+i1];
          bci1 = bc->val[n_v*ic+i1];

          for (i2 = 0; i2 < n_i1; i2++) {
            pout[inod] *= (order * bci1 - i2) / (i2 + 1.0);
          }
        }
      }
    }
  } else {
    fmf_fillC(out, 0.0);

    for (ic = 0; ic < n_coor; ic++) {
      pout = FMF_PtrLevel(out, ic);

      for (inod = 0; inod < n_nod; inod++) {
        for (ii = 0; ii < n_v; ii++) {
          vv = 1.0;
          bcii = bc->val[n_v*ic+ii];

          for (i1 = 0; i1 < n_v; i1++) {
            if (i1 == ii) continue;
            n_i1 = nodes[n_col*inod+i1];
            bci1 = bc->val[n_v*ic+i1];

            for (i2 = 0; i2 < n_i1; i2++) {
              vv *= (order * bci1 - i2) / (i2 + 1.0);
            }

          }

          dval = 0.0;
          n_ii = nodes[n_col*inod+ii];
          for (i1 = 0; i1 < n_ii; i1++) {
            dd = 1.0;

            for (i2 = 0; i2 < n_ii; i2++) {
              if (i1 == i2) continue;

              dd *= (order * bcii - i2) / (i2 + 1.0);
            }
            dval += dd * order / (i1 + 1.0);
          }

          for (ir = 0; ir < dim; ir++) {
            pout[n_nod*ir+inod] += vv * dval * mtx_i->val[n_v*ii+ir];
          }
        }
      }
    }
  }

 end_label:
  return(ret);
}

/*
  Evaluate Lagrange base polynomials in given points on tensor product
  domain.
*/
int32 eval_lagrange_tensor_product(FMField *out, int32 order, int32 diff,
                                   void *_ctx)
{
  LagrangeContext *ctx = (LagrangeContext *) _ctx;
  FMField *bc = ctx->bc;
  FMField *base1d = ctx->base1d;
  int32 *nodes = ctx->nodes;

  int32 ret = RET_OK;
  int32 ii, idim, im, ic;
  int32 nr = out->nRow;
  int32 n_nod = out->nCol;
  int32 dim = bc->nCell;

  fmf_fillC(out, 1.0);

  if (!diff) {
    for (ii = 0; ii < dim; ii++) {
      // slice [:,2*ii:2*ii+2]
      ctx->nodes = nodes + 2 * ii;
      // slice [:,ii:ii+1]
      FMF_SetCell(bc, ii);

      eval_lagrange_simplex(base1d, order, diff, _ctx);

      for (im = 0; im < out->cellSize; im++) {
        out->val[im] *= base1d->val[im];
      }

      ERR_CheckGo(ret);
    }

  } else {
    for (ii = 0; ii < dim; ii++) {
      // slice [:,2*ii:2*ii+2]
      ctx->nodes = nodes + 2 * ii;
      // slice [:,ii:ii+1]
      FMF_SetCell(bc, ii);

      for (idim = 0; idim < dim; idim++) {
        if (ii == idim) {
          eval_lagrange_simplex(base1d, order, diff, _ctx);
        } else {
          eval_lagrange_simplex(base1d, order, 0, _ctx);
        }

        // slice [:,idim:idim+1,:]
        for (im = 0; im < out->nLev; im++) {
          for (ic = 0; ic < n_nod; ic++) {
            out->val[nr*n_nod*im+n_nod*idim+ic] *= base1d->val[n_nod*im+ic];
          }
        }
      }

      ERR_CheckGo(ret);
    }
  }

 end_label:
  // Restore nodes.
  ctx->nodes = nodes;

  return(ret);
}
