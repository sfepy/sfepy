#include "fem.h"
#include "geommech.h"

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
	errset("quadrature point outside of element!");
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

#undef __FUNC__
#define __FUNC__ "evaluate_at"
int32 evaluate_at( FMField *out,
		   int32 *cells, int32 n_cells, int32 n_cells_col,
		   int32 *status, int32 n_status,
		   FMField *dest_coors, FMField *source_vals,
		   int32 *ics, int32 n_ics,
		   int32 *offsets, int32 n_offsets,
		   int32 *iconn0, int32 n_iconn0,
		   FMField *mesh_coors,
		   int32 *nEls0, int32 *nEPs0, int32 **conns0,
		   int32 *nEls, int32 *nEPs, int32 **conns,
		   int32 n_ref_coorss, FMField *ref_coorss,
		   int32 *nNod, int32 *nCol, int32 **nodess,
		   int32 *orders, int32 n_orders,
		   int32 n_mtx_is, FMField *mtx_is,
		   int32 allow_extrapolation,
		   float64 close_limit, float64 qp_eps,
		   int32 i_max, float64 newton_eps )
{
  int32 ii, ie, ie_min, ig, iel, ip, ic, id, ik, dim, nEl, nEP, nGr, dpn;
  int32 nEP_max, n_v_max, n_max;
  int32 order = 0, n_v = 0, ok, ret = RET_OK;
  int32 *conn, *iconn, *nodes = 0;
  float64 aux, err, dist, d_min, vmin, vmax;
  float64 buf16[16], buf16_2[16], buf4[4];
  FMField bc_mtx[1], bc_mtx_i[1], bc_rhs[1];
  FMField e_coors[1], base1d[1], bc[1], dest_point[1], src[1];
  FMField *ref_coors = 0, *mtx_i = 0;
  FMField *bc_max = 0, *b1d_max = 0, *ec_max = 0, *src_max = 0;
  FMField *bf = 0, *bfg = 0, *res = 0, *mtx = 0, *imtx = 0, *xi = 0, *xint = 0;
  FMField *bfs_max, *bfgs_max;
  FMField *bfs, *bfgs;

  dim = mesh_coors->nCol;
  dpn = out->nCol;
  nGr = n_mtx_is;

  e_coors->nAlloc = -1;
  bc->nAlloc = -1;
  base1d->nAlloc = -1;
  dest_point->nAlloc = -1;
  src->nAlloc = -1;
  bc_mtx->nAlloc = -1;
  bc_mtx_i->nAlloc = -1;
  bc_rhs->nAlloc = -1;

  n_max = 0;
  for (ii = 0; ii < (n_offsets - 1); ii++) {
    n_max = Max(n_max, offsets[ii+1] - offsets[ii]);
  }

  /* output("AAA %d %d %d n_max: %d\n", dim, dpn, nGr, n_max); */

  nEP_max = 0;
  n_v_max = 0;
  for (ig = 0; ig < nGr; ig++) {
    nEP_max = Max(nEP_max, nEPs[ig]);
    n_v_max = Max(n_v_max, ref_coorss[ig].nRow);
  }

  bfs = alloc_mem( FMField, n_max );
  bfgs = alloc_mem( FMField, n_max );
  bfs_max = alloc_mem( FMField, n_max );
  bfgs_max = alloc_mem( FMField, n_max );
  for (ie = 0; ie < n_max; ie++) {
    bfs[ie].nAlloc = -1;
    bfgs[ie].nAlloc = -1;
    fmf_alloc( bfs_max + ie, 1, 1, 1, nEP_max );
    fmf_alloc( bfgs_max + ie, 1, 1, dim, nEP_max );
  }

  fmf_createAlloc( &res, 1, 1, 1, dim );
  fmf_createAlloc( &mtx, 1, 1, dim, dim );
  fmf_createAlloc( &imtx, 1, 1, dim, dim );
  fmf_createAlloc( &xint, 1, 1, 1, dim );
  fmf_createAlloc( &xi, n_max, 1, 1, dim );
  fmf_createAlloc( &ec_max, 1, 1, nEP_max, dim );
  fmf_createAlloc( &bc_max, 1, 1, n_v_max, 1 );
  fmf_createAlloc( &b1d_max, 1, 1, 1, nEP_max );
  fmf_createAlloc( &src_max, 1, 1, dpn, nEP_max );

  fmf_fillC( out, 0.0 );

  fmf_pretend( dest_point, dest_coors->nRow, 1, 1, dim, dest_coors->val );
  for (ip = 0; ip < dest_coors->nRow; ip++) {
    ic = ics[ip];

    FMF_SetCell( dest_point, ip );

    ok = 0;
    d_min = 100.0 * 100.0;
    ie_min = -1;
    iconn = iconn0 + 2 * offsets[ic];

    nEl = offsets[ic+1] - offsets[ic];
    /* output("AA %d %d nEl: %d\n", ip, ic, nEl); */
    if (nEl == 0) {
      status[ip] = 3;
      continue;
    }
    
    for (ie = 0; ie < nEl; ie++) {
      ig  = iconn[0];
      iel = iconn[1];
      iconn += 2;

      /* output("BB %d %d %d\n", ie, ig, iel); */

      if (nNod[ig] != nEPs[ig]) {
	errput("incompatible elements!");
      }

      FMF_SetCell( xi, ie );
      
      bf = bfs + ie;
      fmf_pretend( bf, 1, 1, 1, nEPs[ig], bfs_max[ie].val );

      bfg = bfgs + ie;
      fmf_pretend( bfg, 1, 1, dim, nEPs[ig], bfgs_max[ie].val );

      ref_coors = ref_coorss + ig;
      nodes = nodess[ig];
      order = orders[ig];
      mtx_i = mtx_is + ig;
      conn = conns[ig];

      n_v = ref_coors->nRow;

      vmin = ref_coors->val[0];
      vmax = ref_coors->val[dim];

      /* output("BBB %d %d, %d %d %d, %f %f\n", */
      /* 	     n_v, nEPs[ig], */
      /* 	     nNod[ig], nCol[ig], order, vmin, vmax); */

      fmf_pretend( e_coors, 1, 1, nEPs[ig], dim, ec_max->val );

      ele_extractNodalValuesNBN( e_coors, mesh_coors,
				 conn + nEPs[ig] * iel );

      /* fmf_print( e_coors, stdout, 0 ); */
      /* fmf_print( dest_point, stdout, 0 ); */

      if (n_v == (dim + 1)) {
	// Barycentric coordinates.
	fmf_pretend( bc_mtx, 1, 1, n_v, n_v, buf16 );
	fmf_pretend( bc_mtx_i, 1, 1, n_v, n_v, buf16_2 );
	fmf_pretend( bc_rhs, 1, 1, n_v, 1, buf4 );
	for (id = 0; id < dim; id++) {
	  for (ii = 0; ii < n_v; ii++) {
	    bc_mtx->val[n_v*id+ii] = e_coors->val[dim*ii+id];
	  }
	  bc_rhs->val[id] = dest_point->val[id];
	}
	for (ii = 0; ii < n_v; ii++) {
	  bc_mtx->val[n_v*dim+ii] = 1.0;
	}
	bc_rhs->val[dim] = 1.0;

	if (dim == 3) {
	  geme_invert4x4( bc_mtx_i, bc_mtx );
	} else {
	  geme_invert3x3( bc_mtx_i, bc_mtx );
	}

	fmf_pretend( bc, 1, 1, n_v, 1, bc_max->val );
	fmf_mulAB_nn( bc, bc_mtx_i, bc_rhs );
	/* fmf_print( bc, stdout, 0 ); */

	fmf_mulATB_nn( xi, bc, ref_coors );
	
      } else {
	
	fmf_pretend( bc, 1, 1, 2, 1, bc_max->val );
	fmf_pretend( base1d, 1, 1, 1, nEPs[ig], b1d_max->val );

	fmf_fillC( xi, 0.5 * (vmin + vmax) );

	// Newton method
	ii = 0;
	while (ii < i_max) {
	  // Base(xi).
	  eval_lagrange_tensor_product( bf, xi,
					nodes, nNod[ig], nCol[ig],
					Max( order, 1 ), 0,
					mtx_i, bc, base1d,
					1, qp_eps );
	  /* fmf_print( xi, stdout, 0 ); */
	  /* fmf_print( bf, stdout, 0 ); */
	  // X(xi).
	  fmf_mulAB_n1( xint, bf, e_coors );
	  // Rezidual.
	  fmf_subAB_nn( res, dest_point, xint );

	  /* fmf_print( xint, stdout, 0 ); */
	  /* fmf_print( dest_point, stdout, 0 ); */

	  err = 0.0;
	  for (id = 0; id < dim; id++) {
	    err += res->val[id] * res->val[id];
	  }
	  err = sqrt( err );

	  /* output("%d %f\n", ii, err ); */
	  if (err < newton_eps) break;

	  // grad Base(xi).
	  eval_lagrange_tensor_product( bfg, xi,
					nodes, nNod[ig], nCol[ig],
					Max( order, 1 ), 1,
					mtx_i, bc, base1d,
					1, qp_eps );
	  // - Matrix.
	  fmf_mulAB_n1( mtx, bfg, e_coors );

	  geme_invert3x3( imtx, mtx );
	  ERR_CheckGo( ret );
	  
	  fmf_mulAB_nn( xint, res, imtx );
	  fmf_addAB_nn( xi, xi, xint );
	  ii += 1;
	}
        if (ii == i_max) {
          ok = 2;
          // Newton did not converge.
          // Use centre if nothing better is available, but do not spoil
          // possible d_min (see below).
          fmf_fillC( xi, 0.5 * (vmin + vmax) );
        }
      }
      /* fmf_print( xi, stdout, 0 ); */
      /* fmf_print( bf, stdout, 0 ); */

      if (n_v == (dim + 1)) {
	// dist == 0 for 0 <= bc <= 1.
	dist = 0.0;
	for (ii = 0; ii < n_v; ii++) {
	  aux = Min( Max( bc->val[ii] - 1.0, 0.0 ), 100.0 );
	  dist += aux * aux;
	  aux = Min( Max( 0.0 - bc->val[ii], 0.0 ), 100.0 );
	  dist += aux * aux;
	}

      } else {
        if (ok != 2) {
          // dist == 0 for vmin <= xi <= vmax.
          dist = 0.0;
          for (id = 0; id < dim; id++) {
            aux = Min( Max( xi->val[id] - vmax, 0.0 ), 100.0 );
            dist += aux * aux;
            aux = Min( Max( vmin - xi->val[id], 0.0 ), 100.0 );
            dist += aux * aux;
          }
        } else {
          dist = d_min + 1.0;
          ok = 0;
        }
      }
      /* output("CCC %d, %d, %f\n", ie, ii, dist ); */

      if (dist < qp_eps) {
	ok = 1;
	ie_min = ie;
	break;
      } else {
	if (dist < d_min) {
	  d_min = dist;
	  ie_min = ie;
	}
      }
    }

    // Restore ig, iel.
    iconn = iconn0 + 2 * offsets[ic];
    ig  = iconn[2*ie_min+0];
    iel = iconn[2*ie_min+1];

    cells[2*ip+0] = ig;
    cells[2*ip+1] = iel;

    /* output("DDD %d: %d, %f, %d %d\n", ok, ie_min, d_min, ig, iel ); */

    if (!ok) {
      if (allow_extrapolation) {
	// Try using minimum distance xi.
	if (sqrt(d_min) < close_limit) {
	  status[ip] = 1;
	} else {
	  status[ip] = 2;
	}
	FMF_SetCell( xi, ie_min );

	nodes = nodess[ig];
	order = orders[ig];

	if (order > 0) {
	  mtx_i = mtx_is + ig;

	  if (n_v == (dim + 1)) {
	    fmf_pretend( bc, 1, 1, n_v, 1, bc_max->val );
	    eval_lagrange_simplex( bf, xi,
				   nodes, nNod[ig], nCol[ig],
				   order, 0,
				   mtx_i, bc,
				   1, qp_eps );
	  } else {
	    fmf_pretend( bc, 1, 1, 2, 1, bc_max->val );
	    fmf_pretend( base1d, 1, 1, 1, nEPs[ig], b1d_max->val );

	    eval_lagrange_tensor_product( bf, xi,
					  nodes, nNod[ig], nCol[ig],
					  order, 0,
					  mtx_i, bc, base1d,
					  1, qp_eps );
	  }
	}
      } else {
	status[ip] = 3;
      }
    } else {
      status[ip] = 0;

      if ((n_v == (dim + 1)) && (order > 0)) {
	fmf_pretend( bc, 1, 1, n_v, 1, bc_max->val );
	eval_lagrange_simplex( bf, xi,
			       nodes, nNod[ig], nCol[ig],
			       order, 0,
			       mtx_i, bc,
			       1, qp_eps );
      }
    }
    /* output("EEE %d\n", status[ip] ); */

    if (status[ip] <= 1) {
      if (order > 0) {
	conn = conns[ig];
	nEP = nEPs[ig];
      } else {
	conn = conns0[ig];
	nEP = 1;

	fmf_pretend( bf, 1, 1, 1, 1, bfs_max[ie].val );
	fmf_fillC( bf, 1.0 );
      }

      // Interpolate source_vals using bf.
      fmf_pretend( src, 1, 1, dpn, nEP, src_max->val );
      ele_extractNodalValuesDBD( src, source_vals,
				 conn + nEP * iel );

      for (ic = 0; ic < dpn; ic++) {
	aux = 0.0;
	for (ik = 0; ik < nEP; ik++) {
	  aux += bf->val[ik] * src->val[nEP*ic+ik];
	}
	out->val[dpn*ip+ic] = aux;
      }
    }

    /* sys_pause(); */
  }

 end_label:

  free_mem( bfs );
  free_mem( bfgs );

  for (ie = 0; ie < n_max; ie++) {
    fmf_free( bfs_max + ie );
    fmf_free( bfgs_max + ie );
  }
  free_mem( bfs_max );
  free_mem( bfgs_max );

  fmf_freeDestroy( &res );
  fmf_freeDestroy( &mtx );
  fmf_freeDestroy( &imtx );
  fmf_freeDestroy( &xint );
  fmf_freeDestroy( &xi );
  fmf_freeDestroy( &ec_max );
  fmf_freeDestroy( &bc_max );
  fmf_freeDestroy( &b1d_max );
  fmf_freeDestroy( &src_max );

  return( ret );
}
