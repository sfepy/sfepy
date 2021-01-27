/**
    \file contres.cpp
    CONTact RESidual
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
//#include <algorithm>

#include "contres.h"
#include "common.h"

/*! Evaluate shape functions and their 1st partial derivatives of 4-node bilinear element

  \param r - 1st isoparametric (parent, reference) coordinate
  \param s - 2nd isoparametric coordinate

  \return H 1d array (4x1) of shape functions values
  \return dH 2d array (4x2) of 1st partial derivatives of shape functions with respect to r (1st column) and s (2nd column)

*/
void sfd4(double* H, double* dH, double r, double s) {
  const double h1 = 0.25*(1-r)*(1-s);
  const double h2 = 0.25*(1+r)*(1-s);
  const double h3 = 0.25*(1+r)*(1+s);
  const double h4 = 0.25*(1-r)*(1+s);
  /***********************************************/
  const double h1r = -0.25*(1-s);
  const double h2r =  0.25*(1-s);
  const double h3r =  0.25*(s+1);
  const double h4r = -0.25*(1+s);
  /***********************************************/
  const double h1s =  0.25*(r-1);
  const double h2s = -0.25*(1+r);
  const double h3s =  0.25*(r+1);
  const double h4s =  0.25*(1-r);
  /***********************************************/
  H[0] = h1;
  H[1] = h2;
  H[2] = h3;
  H[3] = h4;

  dH[0] = h1r;
  dH[1] = h2r;
  dH[2] = h3r;
  dH[3] = h4r;

  dH[4] = h1s;
  dH[5] = h2s;
  dH[6] = h3s;
  dH[7] = h4s;
}

/*! Evaluate shape functions and their 1st partial derivatives of 8-node (serendipity) quadrilateral element

  \param r - 1st isoparametric (parent, reference) coordinate
  \param s - 2nd isoparametric coordinate

  \return H 1d array (8x1) of shape functions values
  \return dH 2d array (8x2) of 1st partial derivatives of shape functions with respect to r (1st column) and s (2nd column)

*/
void sfd8(float64* H, float64* dH, float64 r, float64 s) {
  const float64 h5 = 0.5*(1-r*r)*(1-s);
  const float64 h6 = 0.5*(1+r)*(1-s*s);
  const float64 h7 = 0.5*(1-r*r)*(1+s);
  const float64 h8 = 0.5*(1-r)*(1-s*s);

  const float64 h1 = 0.25*(1-r)*(1-s) - 0.5*h5 - 0.5*h8;
  const float64 h2 = 0.25*(1+r)*(1-s) - 0.5*h5 - 0.5*h6;
  const float64 h3 = 0.25*(1+r)*(1+s) - 0.5*h6 - 0.5*h7;
  const float64 h4 = 0.25*(1-r)*(1+s) - 0.5*h7 - 0.5*h8;

  /***********************************************/
  const float64 h5r = -r*(1-s);
  const float64 h6r = 0.5*(1-s*s);
  const float64 h7r = -r*(1+s);
  const float64 h8r = -0.5*(1-s*s);

  const float64 h1r = -0.25*(1-s) - 0.5*h5r - 0.5*h8r;
  const float64 h2r =  0.25*(1-s) - 0.5*h5r - 0.5*h6r;
  const float64 h3r =  0.25*(s+1) - 0.5*h6r - 0.5*h7r;
  const float64 h4r = -0.25*(1+s) - 0.5*h7r - 0.5*h8r;
  /***********************************************/
  const float64 h5s = -0.5*(1-r*r);
  const float64 h6s = -(1+r)*s;
  const float64 h7s = 0.5*(1-r*r);
  const float64 h8s = -(1-r)*s;

  const float64 h1s =  0.25*(r-1) - 0.5*h5s - 0.5*h8s;
  const float64 h2s = -0.25*(1+r) - 0.5*h5s - 0.5*h6s;
  const float64 h3s =  0.25*(r+1) - 0.5*h6s - 0.5*h7s;
  const float64 h4s =  0.25*(1-r) - 0.5*h7s - 0.5*h8s;
  /***********************************************/
  H[0] = h1;
  H[1] = h2;
  H[2] = h3;
  H[3] = h4;
  H[4] = h5;
  H[5] = h6;
  H[6] = h7;
  H[7] = h8;

  dH[0] = h1r;
  dH[1] = h2r;
  dH[2] = h3r;
  dH[3] = h4r;
  dH[4] = h5r;
  dH[5] = h6r;
  dH[6] = h7s;
  dH[7] = h8s;

  dH[8]  = h1s;
  dH[9]  = h2s;
  dH[10] = h3s;
  dH[11] = h4s;
  dH[12] = h5s;
  dH[13] = h6s;
  dH[14] = h7s;
  dH[15] = h8s;
}

/*! Evaluate shape functions and their 1st partial derivatives of 6-node (serendipity) triangular element

  \param r - 1st isoparametric coordinate
  \param s - 2nd isoparametric coordinate

  \return H 1d array (6x1) of shape functions values
  \return dH 2d array (6x2) of 1st partial derivatives of shape functions with respect to r (1st column) and s (2nd column)

*/
void sfd6(float64* H, float64* dH, float64 r, float64 s) {
  const float64 h4 = 4 * r*(1 - r - s);
  const float64 h5 = 4 * r*s;
  const float64 h6 = 4 * s*(1 - r - s);

  const float64 h1 = 1 - r - s - 0.5*h4 - 0.5*h6;
  const float64 h2 = r - 0.5*h4 - 0.5*h5;
  const float64 h3 = s - 0.5*h5 - 0.5*h6;
  /***********************************************/
  const float64 h4r = 4 * (1 - r - s) - 4 * r;
  const float64 h5r = 4 * s;
  const float64 h6r = -4 * s;

  const float64 h1r = -1 - 0.5*h4r - 0.5*h6r;
  const float64 h2r = 1 - 0.5*h4r - 0.5*h5r;
  const float64 h3r = 0 - 0.5*h5r - 0.5*h6r;
  /***********************************************/
  const float64 h4s = -4 * r;
  const float64 h5s = 4 * r;
  const float64 h6s = 4 * (1 - r - s) - 4 * s;

  const float64 h1s = -1 - 0.5*h4s - 0.5*h6s;
  const float64 h2s = 0 - 0.5*h4s - 0.5*h5s;
  const float64 h3s = 1 - 0.5*h5s - 0.5*h6s;
  /***********************************************/
  H[0] = h1;
  H[1] = h2;
  H[2] = h3;
  H[3] = h4;
  H[4] = h5;
  H[5] = h6;

  dH[0] = h1r;
  dH[1] = h2r;
  dH[2] = h3r;
  dH[3] = h4r;
  dH[4] = h5r;
  dH[5] = h6r;

  dH[6] = h1s;
  dH[7] = h2s;
  dH[8] = h3s;
  dH[9] = h4s;
  dH[10] = h5s;
  dH[11] = h6s;
}

/*! Evaluate shape functions and their 1st derivatives of 2-node bar element

  \param r - 1st isoparametric coordinate

  \return H 1d array (2x1) of shape functions values
  \return dH 1d array (2x1) of 1st derivatives of shape functions with respect to r

*/
void sfd2(float64* H, float64* dH, float64 r) {
  const float64 h1 = 0.5*(1 - r);
  const float64 h2 = 0.5*(1 + r);
  /***********************************************/
  const float64 h1r = -0.5;
  const float64 h2r = 0.5;
  /***********************************************/
  H[0] = h1;
  H[1] = h2;

  dH[0] = h1r;
  dH[1] = h2r;
}

/*! Calculate contact residual term (gradient) and contact tangent term (Hessian)

  \param len - maximal length of 1d arrays rows,cols, and vals
  \param GPs - 2d array (GPs_len x ??? cols)
  \param ISN - 2d array (nsn*)
  \param IEN -
  \param X -
  \param U -
  \param H -
  \param dH -
  \param gw -
  \param activeGPsOld -
  \param neq - Number of Equations
  \param nsd - Number of Space Dimensions (usually 2 or 3)
  \param npd - Number of Parametric Dimensions (usually 1 or 2) (parametric=reference=parent coordinates)
  \param ngp - Number of GaussPoints on contact segment (segment means face of element)
  \param nes - Number of Element Segments
  \param nsn - Number of Segment Nodes
  \param GPs_len - length of GPs array
  \param epss - penalty parameter (usualy 100*Young's mudulus)
  \param keyContactDetection - if is true, ...
  \param keyAssembleKc - if is true, contact tangent term is assembled only

  \return Gc - 1d array
  \return rows - 1d array
  \return cols - 1d array
  \return valc - 1d array

*/
#undef __FUNC__
#define __FUNC__ "assembleContactResidualAndStiffness"
void assembleContactResidualAndStiffness(float64* Gc, float64* vals, int32* rows, int32* cols, int* len, float64* GPs, int32* ISN, int32* IEN, float64* X, float64* U, float64* H, float64* dH, float64* gw, float64* activeGPsOld, int neq, int nsd, int npd, int ngp, int nes, int nsn, int nen, int GPs_len, float64 epss, int keyContactDetection, int keyAssembleKc)
{
  int ii, i, j, k, g, sdf, pdf;
  int col;
  int* segmentNodesIDs = alloc_mem(int, nsn);
  int* segmentNodesIDm = alloc_mem(int, nsn);
  float64* Xs = alloc_mem(float64, nsn*nsd);
  float64* Xm = alloc_mem(float64, nsn*nsd);
  float64* Us = alloc_mem(float64, nsn*nsd);
  float64* Um = alloc_mem(float64, nsn*nsd);
  float64* dXs = alloc_mem(float64, nsn*npd);
  float64* dxs = alloc_mem(float64, nsn*npd);
  float64* GAPs = alloc_mem(float64, ngp);
  bool* activeGPs = alloc_mem(bool, ngp);
  float64* C_s = alloc_mem(float64, nsn*nsd);
  float64* C_m = alloc_mem(float64, nsn*nsd);
  float64* Hm = alloc_mem(float64, nsn);
  float64* dHm = alloc_mem(float64, nsn*npd);
  float64 Xp[3];
  float64 Xg[3];
  int els, elm, sgs, sgm, IENrows, IENrowm;
  int len_guess = *len;
  float64 r, s, dh;
  float64 Normal[3];
  float64 normal[3];
  float64 dXs_dr1, dXs_dr2, dXs_dr3, dXs_ds1, dXs_ds2, dXs_ds3;
  float64 dxs_dr1, dxs_dr2, dxs_dr3, dxs_ds1, dxs_ds2, dxs_ds3;
  float64 jacobian, normal_length, signGAPs, hs, hm;
  int jdof, jnode, kdof, knode;

  *len = 0;

  // Fill Gc_s array by zeros:
  for (i = 0; i < neq; ++i) {
    Gc[i] = 0.0;
  }

  for (i = 0; i < GPs_len; i += ngp) {

    // Fill C_s and C_m arrays by zeros:
    for (j = 0; j < nsn*nsd; ++j) {
      C_s[j] = 0.0;
      C_m[j] = 0.0;
    }

    // slave element index:
    col = nsd*GPs_len;
    els = (int)GPs[col + i]; // Python numbering starts with 0

    // slave segment index:
    col = (nsd + 1)*GPs_len;
    sgs = (int)GPs[col + i]; // Python numbering starts with 0

    // slave segment coords Xs and displacements Us:
    for (j = 0; j < nsn; ++j) {
      col = nes*j;
      IENrows = ISN[col + sgs]; // Python numbering starts with 0
      col = nen*els;
      segmentNodesIDs[j] = IEN[col + IENrows]; // Python numbering starts with 0
      for (k = 0; k < nsd; ++k) {
	col = k*(int)(neq / nsd);
	Xs[k*nsn + j] = X[col + segmentNodesIDs[j]];
	Us[k*nsn + j] = U[col + segmentNodesIDs[j]];
      }
    }

    for (g = 0; g < ngp; ++g) {
      // Gausspoint gap values and activeGPs:
      activeGPs[g] = (bool) activeGPsOld[i + g];

      if (!activeGPs[g]) {
	continue;
      }

      // master element index:
      col = (2 * nsd + 4)*GPs_len;
      elm = (int)GPs[col + i + g]; // Python numbering starts with 0
      // master segment index:
      col = (2 * nsd + 5)*GPs_len;
      sgm = (int)GPs[col + i + g]; // Python numbering starts with 0

      // master segment coords Xm and displacements Um:
      for (j = 0; j < nsn; ++j) {
	col = nes*j;
	IENrowm = ISN[col + sgm]; // Python numbering starts with 0
	col = nen*elm;
	segmentNodesIDm[j] = IEN[col + IENrowm]; // Python numbering starts with 0
	for (k = 0; k < nsd; ++k) {
	  col = k*(int)(neq / nsd);
	  Xm[k*nsn + j] = X[col + segmentNodesIDm[j]];
	  Um[k*nsn + j] = U[col + segmentNodesIDm[j]];
	}
      }

      // Shape function and its derivatives of gausspoint's master segment

      r = GPs[(nsd + 3)*GPs_len + i + g];
      s = 0.0;

      if (npd == 2) {
	s = GPs[(nsd + 4)*GPs_len + i + g];
      }


      switch (nsn) {
      case 2:
	sfd2(Hm, dHm, r);
	break;
      case 4:
	sfd4(Hm, dHm, r, s);
	break;
      case 6:
	sfd6(Hm, dHm, r, s);
	break;
      case 8:
	sfd8(Hm, dHm, r, s);
      }

      // evaluate gausspoint coords:
      for (sdf = 0; sdf < nsd; ++sdf) {
	Xp[sdf] = 0.0;
	Xg[sdf] = 0.0;
	for (k = 0; k < nsn; ++k) {
	  Xp[sdf] += Hm[k] * (Xm[sdf*nsn + k] + Um[sdf*nsn + k]);
	  Xg[sdf] += H[k*ngp + g] * (Xs[sdf*nsn + k] + Us[sdf*nsn + k]);
	}
      }

      if (keyContactDetection) {
	col = (nsd + 2)*GPs_len;
	GAPs[g] = GPs[col + i + g];
      }
      else {
	GAPs[g] = sqrt((Xp[0] - Xg[0])*(Xp[0] - Xg[0]) + (Xp[1] - Xg[1])*(Xp[1] - Xg[1]) + (Xp[2] - Xg[2])*(Xp[2] - Xg[2]));
      }

      // Fill dXs arrays by zeros:
      for (j = 0; j < nsn*npd; ++j) {
	dXs[j] = 0.0;
	dxs[j] = 0.0;
      }

      // Evaluate tangent vectors:
      for (pdf = 0; pdf < npd; ++pdf) {
	for (sdf = 0; sdf < nsd; ++sdf) {
	  for (j = 0; j < nsn; ++j) {
	    col = j*npd*ngp;
	    dh = dH[col + g*npd + pdf];
	    dXs[npd*sdf + pdf] += dh*Xs[sdf*nsn + j];
	    dxs[npd*sdf + pdf] += dh*(Xs[sdf*nsn + j] + Us[sdf*nsn + j]);
	  }
	}
      }

      // Evaluate normal vector:
      for (ii = 0; ii < 3; ii++) {
        Normal[ii] = 0.0;
        normal[ii] = 0.0;
      }
      if (nsd == 2) {
	Normal[0] = dXs[1];
	Normal[1] = -dXs[0];
	Normal[2] = 0.0;

	normal[0] = dxs[1];
	normal[1] = -dxs[0];
	normal[2] = 0.0;
      }
      else if (nsd == 3) {
	dXs_dr1 = dXs[0];
	dXs_dr2 = dXs[npd + 0];
	dXs_dr3 = dXs[2 * npd + 0];
	dXs_ds1 = dXs[1];
	dXs_ds2 = dXs[npd + 1];
	dXs_ds3 = dXs[2 * npd + 1];
	Normal[0] = dXs_dr2*dXs_ds3 - dXs_dr3*dXs_ds2;
	Normal[1] = dXs_dr3*dXs_ds1 - dXs_dr1*dXs_ds3;
	Normal[2] = dXs_dr1*dXs_ds2 - dXs_dr2*dXs_ds1;

	dxs_dr1 = dxs[0];
	dxs_dr2 = dxs[npd + 0];
	dxs_dr3 = dxs[2 * npd + 0];
	dxs_ds1 = dxs[1];
	dxs_ds2 = dxs[npd + 1];
	dxs_ds3 = dxs[2 * npd + 1];

	normal[0] = (dxs_dr2*dxs_ds3 - dxs_dr3*dxs_ds2);
	normal[1] = (dxs_dr3*dxs_ds1 - dxs_dr1*dxs_ds3);
	normal[2] = (dxs_dr1*dxs_ds2 - dxs_dr2*dxs_ds1);
      }

      jacobian      = sqrt(Normal[0] * Normal[0] + Normal[1] * Normal[1] + Normal[2] * Normal[2]);
      normal_length = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
      normal[0] /= normal_length;
      normal[1] /= normal_length;
      normal[2] /= normal_length;

      if (!keyContactDetection) {
	signGAPs = (Xp[0] - Xg[0])*normal[0] + (Xp[1] - Xg[1])*normal[1] + (Xp[2] - Xg[2])*normal[2];
	if (signGAPs <= 0) {
	  GAPs[g] = -GAPs[g];
	}
	GPs[(nsd + 2)*GPs_len + i + g] = GAPs[g];
      }

      if (GAPs[g] > 0.0) {
	GPs[(2*nsd + 3)*GPs_len + i + g] = 0;
	continue;
      } else {
	GPs[(2*nsd + 3)*GPs_len + i + g] = 1;
      }

      // evaluate shape functions and contact residual vectors:
      for (j = 0; j < nsn; ++j) {
	hs = H[j*ngp + g];
	hm = Hm[j];

	for (sdf = 0; sdf < nsd; ++sdf) {

	  C_s[j*nsd + sdf] += hs*normal[sdf];
	  C_m[j*nsd + sdf] -= hm*normal[sdf];
	  Gc[segmentNodesIDs[j] * nsd + sdf] -= epss*GAPs[g] * hs*normal[sdf] * gw[g] * jacobian;
	}
      }

      if(keyAssembleKc) {
	//if (i < 0.5*GPs_len) {
	for (j = 0; j < nsn*nsd; ++j) { // loop over cols
	  for (k = 0; k < nsn*nsd; ++k) { // loop over rows

	    jdof = j%nsd;
	    jnode = (j - jdof) / nsd; // row node

	    kdof = k%nsd;
	    knode = (k - kdof) / nsd; // col node

	    if(fabs(C_m[j] * C_m[k]) > 1e-50) {
	      if (*len >= len_guess) printf("Error, len is too small: len = %i.\n", len_guess);
	      cols[*len] = segmentNodesIDm[jnode] * nsd + jdof;
	      rows[*len] = segmentNodesIDm[knode] * nsd + kdof;
	      vals[*len] = 0.5 * epss * C_m[j] * C_m[k] * gw[g] * jacobian;
	      (*len)++;
	    }

	    if(fabs(C_s[j] * C_m[k]) > 1e-50) {
	      if (*len >= len_guess) printf("Error, len is too small: len = %i.\n", len_guess);
	      cols[*len] = segmentNodesIDs[jnode] * nsd + jdof;
	      rows[*len] = segmentNodesIDm[knode] * nsd + kdof;
	      vals[*len] = 0.5 * epss * C_s[j] * C_m[k] * gw[g] * jacobian;
	      (*len)++;
	    }

	    if(fabs(C_m[j] * C_s[k]) > 1e-50) {
	      if (*len >= len_guess) printf("Error, len is too small: len = %i.\n", len_guess);
	      cols[*len] = segmentNodesIDm[jnode] * nsd + jdof;
	      rows[*len] = segmentNodesIDs[knode] * nsd + kdof;
	      vals[*len] = 0.5 * epss * C_m[j] * C_s[k] * gw[g] * jacobian;
	      (*len)++;
	    }

	    if(fabs(C_s[j] * C_s[k]) > 1e-50) {
	      if (*len >= len_guess) printf("Error, len is too small: len = %i.\n", len_guess);
	      cols[*len] = segmentNodesIDs[jnode] * nsd + jdof;
	      rows[*len] = segmentNodesIDs[knode] * nsd + kdof;
	      vals[*len] = 0.5 * epss * C_s[j] * C_s[k] * gw[g] * jacobian;
	      (*len)++;
	    }
	  }
	}
      }
      // Fill C_m array by zeros:
      for (j = 0; j < nsn*nsd; ++j) {
	C_s[j] = 0.0;
	C_m[j] = 0.0;
      }
      //}
    } // loop over gausspoints
  } // loop over GPs rows

  free_mem(segmentNodesIDs);
  free_mem(segmentNodesIDm);
  free_mem(Xs);
  free_mem(Xm);
  free_mem(Us);
  free_mem(Um);
  free_mem(dXs);
  free_mem(dxs);
  free_mem(GAPs);
  free_mem(activeGPs);
  free_mem(C_s);
  free_mem(C_m);
  free_mem(Hm);
  free_mem(dHm);
}

#undef __FUNC__
#define __FUNC__ "getLongestEdgeAndGPs"
void getLongestEdgeAndGPs(float64* longestEdge, float64* GPs, int n, int nsd, int ngp, int neq, int nsn, int nes, int nen, uint32* elementID, uint32* segmentID, int32* ISN, int32* IEN, float64* H, float64* X) {
  // GPs legend:              Xg    els  sgs   gap  Xm    isActive elm sgm
  //float64* GPs = alloc_mem(float64, n*(nsd + 1 + 1  + 1  + nsd + 1       + 1 + 1));
  int e, i, j, sdf, el, sg, IENrow;
  int* segmentNodesID = alloc_mem(int, nsn);
  float64 lengthOfEdge;
  float64* Xs = alloc_mem(float64, nsn*nsd);
  float64* Xg = alloc_mem(float64, ngp*nsd);
  int g = 0;
  *longestEdge = 0.0;

  for (e = 0; e < n; ++e) {
    el = elementID[e];
    sg = segmentID[e];

    // segment coords Xs:
    for (i = 0; i < nsn; ++i) {
      IENrow = ISN[nes*i + sg]; // Python numbering starts with 0
      segmentNodesID[i] = IEN[nen*el + IENrow]; // Python numbering starts with 0
      for (j = 0; j < nsd; ++j) {
	Xs[j*nsn + i] = X[j*(int)(neq / nsd) + segmentNodesID[i]];
      }
    }

    // evaluate gausspoint coords:
    for (i = 0; i < ngp; ++i) {
      for (sdf = 0; sdf < nsd; ++sdf) {
	Xg[i*nsd + sdf] = 0.0;
	for (j = 0; j < nsn; ++j) {
	  Xg[i*nsd + sdf] += H[j*ngp + i] * Xs[sdf*nsn + j];
	}

	GPs[sdf*n*ngp + g] = Xg[i*nsd + sdf]; // slave gausspoint coords
	GPs[(nsd + 3 + sdf)*n*ngp + g] = 0.0; // init baricentric coords
      }
      GPs[nsd*n*ngp + g] = el;                // slave element
      GPs[(nsd + 1)*n*ngp + g] = sg;          // slave segment
      GPs[(nsd + 2)*n*ngp + g] = FLT_MAX;     // init gap
      GPs[(2 * nsd + 3)*n*ngp + g] = 0;       // init is NO active
      GPs[(2 * nsd + 4)*n*ngp + g] = 0;       // master element
      GPs[(2 * nsd + 5)*n*ngp + g] = 0;       // master segment

      g++;
    }

    for (i = 0; i < nsn; ++i) {
      for (j = i+1; j < nsn; ++j) {
	lengthOfEdge = 0.0;
	for (sdf = 0; sdf < nsd; ++sdf) {
	  lengthOfEdge += pow(Xs[sdf*nsn + i] - Xs[sdf*nsn + j ], 2);
	}
	*longestEdge = Max(*longestEdge, sqrt(lengthOfEdge) );
      }
    }

  } // loop over elements

  free_mem(segmentNodesID);
  free_mem(Xs);
  free_mem(Xg);

}

#undef __FUNC__
#define __FUNC__ "getAABB"
void getAABB(float64* AABBmin, float64* AABBmax, int nsd, int nnod, float64* X, float64 longestEdge, int32* IEN, int32* ISN, uint32* elementID, uint32* segmentID, int n, int nsn, int nes, int nen, int neq) {
  int e, i, sdf, el, sg;
  int* segmentNodesID = alloc_mem(int, nsn);
  int IENrow;
  float64 x;

  for (sdf = 0; sdf < nsd; ++sdf) {
    AABBmin[sdf] = FLT_MAX;
    AABBmax[sdf] = -FLT_MAX;

    for (e = 0; e < n; ++e) {
      el = elementID[e]; // Python numbering starts with 0
      sg = segmentID[e]; // Python numbering starts with 0

      // segment coords Xs:
      for (i = 0; i < nsn; ++i) {
	IENrow = ISN[nes*i + sg]; // Python numbering starts with 0
	segmentNodesID[i] = IEN[nen*el + IENrow]; // Python numbering starts with 0
	x = X[sdf*(int)(neq / nsd) + segmentNodesID[i]];
	AABBmin[sdf] = Min(AABBmin[sdf], x);
	AABBmax[sdf] = Max(AABBmax[sdf], x);
      }
    }

    if ((AABBmax[sdf] - AABBmin[sdf]) < longestEdge) {
      AABBmax[sdf] += 0.5*longestEdge;
      AABBmin[sdf] -= 0.5*longestEdge;
    }
  }
  free_mem(segmentNodesID);
}

#undef __FUNC__
#define __FUNC__ "evaluateContactConstraints"
void evaluateContactConstraints(float64* GPs, int32* ISN, int32* IEN, int32* N, float64* AABBmin, float64* AABBmax, int32* head, int32* next, float64* X, uint32* elementID, uint32* segmentID, int n, int nsn, int nsd, int npd, int ngp, int nen, int nes, int neq, float64 longestEdge) {
  int ii, i, e, k, j, it, sdf, i2, i1, i0;
  int* segmentNodesID = alloc_mem(int, nsn);
  float64* Xm = alloc_mem(float64, nsn*nsd);
  float64* Xmin = alloc_mem(float64, nsd);
  float64* Xmax = alloc_mem(float64, nsd);
  float64* Hm = alloc_mem(float64, nsn);
  float64* dHm = alloc_mem(float64, nsn*npd);
  float64 Xt[9];
  float64 Xc[3];
  int el, sg, IENrow, Ic, v, els, sgs, niter, max_niter;
  int Imin[3];
  int Imax[3];
  float64 normal[3];
  float64 t1[3];
  float64 t2[3];
  float64 t3[3];
  float64 Xg[3];
  float64 Xp[3];
  float64 ra[9];
  float64 Q1[3];
  float64 Q2[3];
  float64 Q3[3];
  float64 normalLength, d_tmp, d, t1_norm, sign, Q1n, Q2n, Q3n;
  float64 b1, b2, A11, A22, A12;
  float64 recDetA, invA11, invA22, invA12;
  float64 r, s, r_len, s_len, dr, ds, x, dx_dr, dx_ds, dr_norm;
  bool isInside;

  // If segment element is quad then it is divided to 4 triangles:
  int ntr = 1;
  if(nsn == 4 || nsn == 8) {
    ntr = 4;
  }

  for (i = 0; i < n*ngp; ++i) {
    GPs[(nsd + 2)*n*ngp + i] = FLT_MAX;
  }

  for (e = 0; e < n; ++e) {
    el = elementID[e];
    sg = segmentID[e];

    // segment coords Xm:
    for (k = 0; k < nsd; ++k) {

      Xmin[k] = FLT_MAX;
      Xmax[k] = -FLT_MAX;
      Xc[k] = 0;

      for (j = 0; j < nsn; ++j) {
	IENrow = ISN[nes*j + sg]; // Python numbering starts with 0
	segmentNodesID[j] = IEN[nen*el + IENrow]; // Python numbering starts with 0
	Xm[k*nsn + j] = X[k*(int)(neq / nsd) + segmentNodesID[j]];
	Xmin[k] = Min(Xmin[k], Xm[k*nsn + j]);
	Xmax[k] = Max(Xmax[k], Xm[k*nsn + j]);
	Xc[k] += Xm[k*nsn + j];
      }
      Xmin[k] -= 0.5*longestEdge;
      Xmax[k] += 0.5*longestEdge;
      Xc[k] /= nsn;
    }

    // Loop over segment triangles:
    for (it = 0; it < ntr; ++it) {

      if(nsd == 3) {
	// triangle coords Xt:
	if(ntr == 1) {
	  Xt[0] = Xm[0];
	  Xt[1] = Xm[1];
	  Xt[2] = Xm[2];

	  Xt[3] = Xm[nsn];
	  Xt[4] = Xm[nsn+1];
	  Xt[5] = Xm[nsn+2];

	  Xt[6] = Xm[2*nsn];
	  Xt[7] = Xm[2*nsn+1];
	  Xt[8] = Xm[2*nsn+2];
	}
	else if(ntr == 4) {
	  Xt[0] = Xm[it];
	  Xt[1] = Xm[(it+1)%ntr];
	  Xt[2] = Xc[0];

	  Xt[3] = Xm[nsn+it];
	  Xt[4] = Xm[nsn+(it+1)%ntr];
	  Xt[5] = Xc[1];

	  Xt[6] = Xm[2*nsn+it];
	  Xt[7] = Xm[2*nsn+(it+1)%ntr];
	  Xt[8] = Xc[2];
	}

	// Min and Max of the triangle coords:
	for (k = 0; k < nsd; ++k) {
	  Xmin[k] = FLT_MAX;
	  Xmax[k] = -FLT_MAX;
	  for (j = 0; j < 3; ++j) {
	    Xmin[k] = Min(Xmin[k], Xt[k*3 + j]);
	    Xmax[k] = Max(Xmax[k], Xt[k*3 + j]);
	  }
	  Xmin[k] -= 0.5*longestEdge;
	  Xmax[k] += 0.5*longestEdge;
	}
      }

      for (ii = 0; ii < 3; ii++) {
        Imin[ii] = 0;
        Imax[ii] = 0;
        normal[ii] = 0.0;
        t1[ii] = 0.0;
        t2[ii] = 0.0;
        t3[ii] = 0.0;
        Xg[ii] = 0.0;
        Xp[ii] = 0.0;
      }

      Imin[0] = (int)(N[0] * (Xmin[0] - AABBmin[0]) / (AABBmax[0] - AABBmin[0]));
      Imin[1] = (int)(N[1] * (Xmin[1] - AABBmin[1]) / (AABBmax[1] - AABBmin[1]));
      Imax[0] = (int)(N[0] * (Xmax[0] - AABBmin[0]) / (AABBmax[0] - AABBmin[0]));
      Imax[1] = (int)(N[1] * (Xmax[1] - AABBmin[1]) / (AABBmax[1] - AABBmin[1]));

      if (nsd == 2) {

	Imin[2] = 0;
	Imax[2] = 0;

	t1[0] = Xm[1] - Xm[0];
	t1[1] = Xm[3] - Xm[2];
	t1[2] = 0.0;

	normal[0] = t1[1];
	normal[1] = -t1[0];
	normal[2] = 0.0;
      }
      else if (nsd == 3) {

	Imin[2] = (int)(N[2] * (Xmin[2] - AABBmin[2]) / (AABBmax[2] - AABBmin[2]));
	Imax[2] = (int)(N[2] * (Xmax[2] - AABBmin[2]) / (AABBmax[2] - AABBmin[2]));

	// Tangent vectors parallel with element edges 1 and 2:
	// Component: X     Y      Z
	// Vertex 1:  Xt[0] Xt[3]  Xt[6]
	// Vertex 2:  Xt[1] Xt[4]  Xt[7]
	// Vertex 3:  Xt[2] Xt[5]  Xt[8]

	t1[0] = Xt[1] - Xt[0];
	t1[1] = Xt[4] - Xt[3];
	t1[2] = Xt[7] - Xt[6];

	t2[0] = Xt[2] - Xt[1];
	t2[1] = Xt[5] - Xt[4];
	t2[2] = Xt[8] - Xt[7];

	t3[0] = Xt[0] - Xt[2];
	t3[1] = Xt[3] - Xt[5];
	t3[2] = Xt[6] - Xt[8];

	// Normal vector:
	normal[0] = t1[1] * t2[2] - t1[2] * t2[1];
	normal[1] = t1[2] * t2[0] - t1[0] * t2[2];
	normal[2] = t1[0] * t2[1] - t1[1] * t2[0];
      }

      normalLength = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
      normal[0] /= normalLength;
      normal[1] /= normalLength;
      normal[2] /= normalLength;

      for (sdf = 0; sdf < nsd; ++sdf) {
	if (Imin[sdf] < 0) {
	  Imin[sdf] = 0;
	}
	if (Imax[sdf] < 0) {
	  Imax[sdf] = 0;
	}
	if (Imin[sdf] >= N[sdf]) {
	  Imin[sdf] = N[sdf] - 1;
	}
	if (Imax[sdf] >= N[sdf]) {
	  Imax[sdf] = N[sdf] - 1;
	}
      }

      for (i2 = Imin[2]; i2 <= Imax[2]; ++i2) {
	for (i1 = Imin[1]; i1 <= Imax[1]; ++i1) {
	  for (i0 = Imin[0]; i0 <= Imax[0]; ++i0) {
	    Ic = i2*N[0] * N[1] + i1*N[0] + i0;
	    v = head[Ic];

	    while (v != -1) {

	      // Jump if Gausspoit segment is equal to master segment
	      els = GPs[nsd*n*ngp + v];            // slave element
	      sgs = GPs[(nsd + 1)*n*ngp + v];      // slave segment
	      if (el == els && sg == sgs) {
		v = next[v];
		continue;
	      }

	      // Inside-outside algorithm:
	      isInside = 0;
	      if (nsd == 2) {
		d = 0.0;
		t1_norm = 0.0;
		for (i = 0; i < nsd; ++i) {
		  Xg[i] = GPs[i*n*ngp + v];
		  ra[i] = Xg[i] - Xm[i*nsn + 0];
		  d += ra[i] * t1[i];
		  t1_norm += t1[i] * t1[i];
		}
		t1_norm = sqrt(t1_norm);
		d = d / t1_norm;
		// Check if inside edge1:
		if (d >= 0.0 && d <= t1_norm) {
		  isInside = 1;
		  d_tmp = d;
		  d = 0.0;
		  sign = 0.0;
		  for (i = 0; i < nsd; ++i) {
		    Xp[i] = Xm[i*nsn + 0] + d_tmp * t1[i] / t1_norm;
		    sign += (Xg[i] - Xp[i])*normal[i];
		    d += pow(Xg[i] - Xp[i], 2);
		  }
		  d = sqrt(d);
		  if (sign < 0) {
		    d *= -1;
		  }
		}
	      }
	      else if (nsd == 3) {
		d = 0.0;
		for (i = 0; i < nsd; ++i) {
		  Xg[i] = GPs[i*n*ngp + v];
		  ra[i * 3 + 0] = Xg[i] - Xt[i*3 + 0];
		  ra[i * 3 + 1] = Xg[i] - Xt[i*3 + 1];
		  ra[i * 3 + 2] = Xg[i] - Xt[i*3 + 2];
		}

		// component:  X    Y    Z
		// r1:       ra[0] ra[3] ra[6]
		// r2:       ra[1] ra[4] ra[7]
		// r3:       ra[2] ra[5] ra[8]
		Q1[0] = ra[3] * t1[2] - ra[6] * t1[1];
		Q1[1] = ra[6] * t1[0] - ra[0] * t1[2];
		Q1[2] = ra[0] * t1[1] - ra[3] * t1[0];

		Q2[0] = ra[4] * t2[2] - ra[7] * t2[1];
		Q2[1] = ra[7] * t2[0] - ra[1] * t2[2];
		Q2[2] = ra[1] * t2[1] - ra[4] * t2[0];

		Q3[0] = ra[5] * t3[2] - ra[8] * t3[1];
		Q3[1] = ra[8] * t3[0] - ra[2] * t3[2];
		Q3[2] = ra[2] * t3[1] - ra[5] * t3[0];

		Q1n = Q1[0] * normal[0] + Q1[1] * normal[1] + Q1[2] * normal[2];
                Q2n = Q2[0] * normal[0] + Q2[1] * normal[1] + Q2[2] * normal[2];
                Q3n = Q3[0] * normal[0] + Q3[1] * normal[1] + Q3[2] * normal[2];

		if (Q1n*Q2n >= 0) {
		  if (Q1n*Q3n >= 0) {
		    isInside = 1;
		    d = ra[0] * normal[0] + ra[3] * normal[1] + ra[6] * normal[2];
		    for (i = 0; i < nsd; ++i) {
		      Xp[i] = Xg[i] - d*normal[i];
		    }
		  }
		}
	      } // else if (nsd = 3)

	      /* printf("Xg = (%f, %f, %f)\n", Xg[0], Xg[1], Xg[2]); */
	      /* printf("Xp = (%f, %f, %f)\n", Xp[0], Xp[1], Xp[2]); */
	      /* printf("normal = (%f, %f, %f)\n", normal[0], normal[1], normal[2]); */

	      if (isInside) {
		// If distance is less then current closest distance:
		if (d < GPs[(nsd + 2)*n*ngp + v]) {

		  // Initial guess of the parametric coordinates on the triangle:
		  r_len = 0, r = 0;
		  s_len = 0, s = 0;
		  switch (nsn) {
		  case 2:
		    // Tangent vectors parallel with element edges 1 and 2:
		    // Component: X     Y      Z
		    // Node 1:  Xm[0] Xm[2]  Xm[4]
		    // Node 2:  Xm[1] Xm[3]  Xm[5]

		    r_len = pow(Xm[1] - Xm[0], 2.0) +
		      pow(Xm[3] - Xm[2], 2.0) +
		      pow(Xm[5] - Xm[4], 2.0);

		    r = ((Xp[0] - Xm[0])  * (Xm[1] - Xm[0]) +
			 (Xp[1] - Xm[2])  * (Xm[3] - Xm[2]) +
			 (Xp[2] - Xm[4])  * (Xm[5] - Xm[4])) / r_len;
		    r = 2*r-1;
		    break;
		  case 6:
		    r_len = pow(Xm[1] - Xm[0], 2.0) +
		      pow(Xm[7] - Xm[6], 2.0) +
		      pow(Xm[13] - Xm[12], 2.0);

		    s_len = pow(Xm[2] - Xm[0], 2.0) +
		      pow(Xm[8] - Xm[6], 2.0) +
		      pow(Xm[14] - Xm[12], 2.0);

		    r = ((Xp[0] - Xm[0])  * (Xm[1] - Xm[0]) +
			 (Xp[1] - Xm[6])  * (Xm[7] - Xm[6]) +
			 (Xp[2] - Xm[12]) * (Xm[13] - Xm[12])) / r_len;

		    s = ((Xp[0] - Xm[0])  * (Xm[2] - Xm[0]) +
			 (Xp[1] - Xm[6])  * (Xm[8] - Xm[6]) +
			 (Xp[2] - Xm[12]) * (Xm[14] - Xm[12])) / s_len;
		    break;
		  case 8:
		    r_len = pow(Xm[1]  - Xm[0],  2.0) +
		      pow(Xm[9]  - Xm[8],  2.0) +
		      pow(Xm[17] - Xm[16], 2.0);

		    s_len = pow(Xm[3]  - Xm[0],  2.0) +
		      pow(Xm[11] - Xm[8],  2.0) +
		      pow(Xm[19] - Xm[16], 2.0);

		    r = ((Xp[0] - Xm[0])  * (Xm[1] - Xm[0]) +
			 (Xp[1] - Xm[8])  * (Xm[9] - Xm[8]) +
			 (Xp[2] - Xm[16]) * (Xm[17] - Xm[16])) / r_len;

		    r = 2*r-1;

		    s = ((Xp[0] - Xm[0])  * (Xm[3]  - Xm[0]) +
			 (Xp[1] - Xm[8])  * (Xm[11] - Xm[8]) +
			 (Xp[2] - Xm[16]) * (Xm[19] - Xm[16])) / s_len;

		    s = 2*s-1;
		  }


		  // Local contact search by Least-square projection method:
		  dr_norm = FLT_MAX;
		  niter = 0;
		  max_niter = 1000;
		  do {
		    switch (nsn) {
		    case 2:
		      sfd2(Hm, dHm, r);
		      break;
                    case 4:
		      sfd4(Hm, dHm, r, s);
		      break;
		    case 6:
		      sfd6(Hm, dHm, r, s);
		      break;
		    case 8:
		      sfd8(Hm, dHm, r, s);
		    }

		    A11 = 0.0;
		    A22 = 0.0;
		    A12 = 0.0;
		    b1 = 0.0;
		    b2 = 0.0;
		    d_tmp = 0.0;

		    for (sdf = 0; sdf < nsd; ++sdf) {
		      x = 0.0;
		      dx_dr = 0.0;
		      dx_ds = 0.0;

		      for (k = 0; k < nsn; ++k) {
			x += Hm[k] * Xm[sdf*nsn + k];
			dx_dr += dHm[k] * Xm[sdf*nsn + k];
			if (npd == 2) {
			  dx_ds += dHm[nsn + k] * Xm[sdf*nsn + k];
			}
		      }

		      b1 += dx_dr*(Xg[sdf] - x);
		      A11 += dx_dr*dx_dr;

		      if (npd == 2) {
			b2 += dx_ds*(Xg[sdf] - x);
			A22 += dx_ds*dx_ds;
			A12 += dx_dr*dx_ds;
		      }
		      d_tmp += (Xg[sdf] - x) * (Xg[sdf] - x);
		    }

		    d = (d<0) ? -sqrt(d_tmp) : sqrt(d_tmp);

		    if (npd == 1) {
		      invA11 = 1 / A11;
		      dr = invA11*b1;
		      r += dr;
		      dr_norm = fabs(dr);
		    }

		    if (npd == 2) {
		      recDetA = 1 / (A11*A22 - A12*A12);
		      invA11 = recDetA * A22;
		      invA22 = recDetA * A11;
		      invA12 = -recDetA * A12;
		      dr = invA11*b1 + invA12*b2;
		      ds = invA12*b1 + invA22*b2;

		      r += dr;
		      s += ds;
		      dr_norm = sqrt(dr*dr + ds*ds);
		    }

		    niter++;
		  } while (dr_norm > 1e-3 && niter < max_niter);

		  if (niter >= max_niter) {
		    printf("Fatal error: Local contact search do NOT converge.\n");
		  }

		  GPs[(nsd       + 2)*n*ngp + v] = d;      // store penetration
		  if(d < 0) {
		    GPs[(2 * nsd + 3)*n*ngp + v] = 1.0;    // set gausspoint to active state
		    GPs[(2 * nsd + 4)*n*ngp + v] = el; // set master element
		    GPs[(2 * nsd + 5)*n*ngp + v] = sg; // set master segment
		    GPs[(    nsd + 3)*n*ngp + v] = r;
		    if (npd == 2) {
		      GPs[(  nsd + 4)*n*ngp + v] = s;
		    }
		  }
		} // if d is less then
	      } // if inside

	      v = next[v];
	    } // while
	  } // i0
	} // i1
      }	// i2
    } // loop over triangles
  } // loop over elements

  /*
  printf("GPs = \n");
  for(int row = 0; row<(n*ngp); ++row) {
    for(int col = 0; col<(2*nsd+6); ++col) {
      printf("%f\t", GPs[col*n*ngp + row]);
    }
    printf("\n");
  }
  */

  free_mem(segmentNodesID);
  free_mem(Xm);
  free_mem(Xmin);
  free_mem(Xmax);
  free_mem(Hm);
  free_mem(dHm);
}
