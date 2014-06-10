#include "nurbs.h"

inline void ravel_multi_index(uint32 *index, uint32 *indices,
                              uint32 *shape, uint32 num)
{
  uint32 ii, stride = 1;
  uint32 raveled = 0;

  for (ii = num - 1; ii >= 1; ii--) {
    raveled += stride * indices[ii];
    stride *= shape[ii - 1];
  }
  raveled += stride * indices[0];

  *index = raveled;
}

inline void unravel_index(uint32 *indices, uint32 index,
                          uint32 *shape, uint32 num)
{
  int32 ii; // To iterate to zero...

  for (ii = num - 1; ii >= 0; ii--) {
    indices[ii] = index % shape[ii];
    index /= shape[ii];
  }
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
