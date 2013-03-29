#include "lobatto.h"

/*
  Evaluate 1D Lobatto functions of the given order in given points.
*/
int32 eval_lobatto1d(FMField *out, FMField *coors, int32 order)
{
  int32 ret = RET_OK;
  int32 ii;
  fun eval_fun;
  int32 n_coor = coors->nRow;

  if ((order < 0) || (order > max_order)) {
    errput("order must be in [0, %d]! (was %d)", max_order, order);
    errset("wrong order value (see above!");
    ERR_CheckGo(ret);
  }

  eval_fun = lobatto[order];
  for (ii = 0; ii < n_coor; ii++) {
    out->val[ii] = eval_fun(coors->val[ii]);
  }

 end_label:
  return(ret);
}

/*
  Evaluate tensor product Lobatto functions of the given order in given points.
*/
int32 eval_lobatto_tensor_product(FMField *out, FMField *coors,
                                  int32 *nodes,
                                  float64 cmin, float64 cmax,
                                  int32 diff)
{
  int32 ret = RET_OK;
  float64 dx = cmax - cmin;
  int32 ii, ifun, ic, ir, io;
  int32 n_coor = coors->nRow;
  int32 dim = coors->nCol;
  int32 n_fun = out->nCol;
  FMField xis[1];
  fun eval_fun;

  fmf_alloc(xis, 1, 1, n_coor, dim);

  for (ii = 0; ii < (n_fun * dim); ii++) {
    if (nodes[ii] > max_order) {
      errput("order must be in [0, %d]! (was %d)", max_order, nodes[ii]);
      errset("wrong order value (see above!");
      ERR_CheckGo(ret);
    }
  }

  // Transform coordinates via affine coordinates lambda to be in [-1, 1].
  for (ii = 0; ii < (n_coor * dim); ii++) {
    xis->val[ii] = 2.0 * (coors->val[ii] - cmin) / dx - 1.0;
  }


  fmf_fillC(out, 1.0);

  if (!diff) {
    for (ii = 0; ii < dim; ii++) {
      for (ifun = 0; ifun < n_fun; ifun++) {
        eval_fun = lobatto[nodes[dim * ifun + ii]];
        for (ic = 0; ic < n_coor; ic++) {
          out->val[n_fun * ic + ifun] *= eval_fun(xis->val[dim * ic + ii]);
        }
      }
    }
  } else {
    for (ii = 0; ii < dim; ii++) {
      for (ifun = 0; ifun < n_fun; ifun++) {
        for (ir = 0; ir < dim; ir++) {
          if (ir == ii) {
            eval_fun = d_lobatto[nodes[dim * ifun + ii]];
          } else {
            eval_fun = lobatto[nodes[dim * ifun + ii]];
          }
          for (ic = 0; ic < n_coor; ic++) {
            io = n_fun * (dim * ic + ir) + ifun;
            out->val[io] *= eval_fun(xis->val[dim * ic + ii]);
          }
        }
      }
    }
    // Multiply by 2 due to the transformation of coordinates.
    fmf_mulC(out, 2.0);
  }

 end_label:
  fmf_free(xis);

  return(ret);
}
