/**
    \file contres.cpp
    CONTact RESidual
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cfloat>
#include <algorithm>

#include "contres.h"

/*! Evaluate shape functions and their 1st partial derivatives of 8-node (serendipity) quadrilateral element

  \param r - 1st isoparametric (parent, reference) coordinate
  \param s - 2nd isoparametric coordinate

  \return H 1d array (8x1) of shape functions values
  \return dH 2d array (8x2) of 1st partial derivatives of shape functions with respect to r (1st column) and s (2nd column)

*/
void sfd8(double* H, double* dH, double r, double s) {
  const double h5 = 0.5*(1-r*r)*(1-s);
  const double h6 = 0.5*(1+r)*(1-s*s);
  const double h7 = 0.5*(1-r*r)*(1+s);
  const double h8 = 0.5*(1-r)*(1-s*s);

  const double h1 = 0.25*(1-r)*(1-s) - 0.5*h5 - 0.5*h8;
  const double h2 = 0.25*(1+r)*(1-s) - 0.5*h5 - 0.5*h6;
  const double h3 = 0.25*(1+r)*(1+s) - 0.5*h6 - 0.5*h7;
  const double h4 = 0.25*(1-r)*(1+s) - 0.5*h7 - 0.5*h8;

  /***********************************************/
  const double h5r = -r*(1-s);
  const double h6r = 0.5*(1-s*s);
  const double h7r = -r*(1+s);
  const double h8r = -0.5*(1-s*s);

  const double h1r = -0.25*(1-s) - 0.5*h5r - 0.5*h8r;
  const double h2r =  0.25*(1-s) - 0.5*h5r - 0.5*h6r;
  const double h3r =  0.25*(s+1) - 0.5*h6r - 0.5*h7r;
  const double h4r = -0.25*(1+s) - 0.5*h7r - 0.5*h8r;
  /***********************************************/
  const double h5s = -0.5*(1-r*r);
  const double h6s = -(1+r)*s;
  const double h7s = 0.5*(1-r*r);
  const double h8s = -(1-r)*s;

  const double h1s =  0.25*(r-1) - 0.5*h5s - 0.5*h8s;
  const double h2s = -0.25*(1+r) - 0.5*h5s - 0.5*h6s;
  const double h3s =  0.25*(r+1) - 0.5*h6s - 0.5*h7s;
  const double h4s =  0.25*(1-r) - 0.5*h7s - 0.5*h8s;
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
void sfd6(double* H, double* dH, double r, double s) {
  const double h4 = 4 * r*(1 - r - s);
  const double h5 = 4 * r*s;
  const double h6 = 4 * s*(1 - r - s);

  const double h1 = 1 - r - s - 0.5*h4 - 0.5*h6;
  const double h2 = r - 0.5*h4 - 0.5*h5;
  const double h3 = s - 0.5*h5 - 0.5*h6;
  /***********************************************/
  const double h4r = 4 * (1 - r - s) - 4 * r;
  const double h5r = 4 * s;
  const double h6r = -4 * s;

  const double h1r = -1 - 0.5*h4r - 0.5*h6r;
  const double h2r = 1 - 0.5*h4r - 0.5*h5r;
  const double h3r = 0 - 0.5*h5r - 0.5*h6r;
  /***********************************************/
  const double h4s = -4 * r;
  const double h5s = 4 * r;
  const double h6s = 4 * (1 - r - s) - 4 * s;

  const double h1s = -1 - 0.5*h4s - 0.5*h6s;
  const double h2s = 0 - 0.5*h4s - 0.5*h5s;
  const double h3s = 1 - 0.5*h5s - 0.5*h6s;
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
void sfd2(double* H, double* dH, double r) {
  const double h1 = 0.5*(1 - r);
  const double h2 = 0.5*(1 + r);
  /***********************************************/
  const double h1r = -0.5;
  const double h2r = 0.5;
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
void assembleContactResidualAndStiffness(double* Gc, double* vals, double* rows, double* cols, int* len, double* GPs, int* ISN, int* IEN, double* X, double* U, double* H, double* dH, double* gw, double* activeGPsOld, int neq, int nsd, int npd, int ngp, int nes, int nsn, int nen, int GPs_len, double epss, bool keyContactDetection, bool keyAssembleKc)
{
  int col;
  int* segmentNodesIDs = new int[nsn];
  int* segmentNodesIDm = new int[nsn];
  double* Xs = new double[nsn*nsd];
  double* Xm = new double[nsn*nsd];
  double* Us = new double[nsn*nsd];
  double* Um = new double[nsn*nsd];
  double* dXs = new double[nsn*npd];
  double* dxs = new double[nsn*npd];
  double* GAPs = new double[ngp];
  bool* activeGPs = new bool[ngp];
  double* C_s = new double[nsn*nsd];
  double* C_m = new double[nsn*nsd];
  double* Hm = new double[nsn];
  double* dHm = new double[nsn*npd];

  int len_guess = *len;
  *len = 0;

  double Xp[3];
  double Xg[3];

  // Fill Gc_s array by zeros:
  for (int i = 0; i < neq; ++i) {
    Gc[i] = 0.0;
  }

  for (int i = 0; i < GPs_len; i += ngp) {

    // Fill C_s and C_m arrays by zeros:
    for (int j = 0; j < nsn*nsd; ++j) {
      C_s[j] = 0.0;
      C_m[j] = 0.0;
    }

    // slave element index:
    col = nsd*GPs_len;
    const int els = (int)GPs[col + i] - 1; // Matlab numbering starts with 1

    // slave segment index:
    col = (nsd + 1)*GPs_len;
    const int sgs = (int)GPs[col + i] - 1; // Matlab numbering starts with 1

    // slave segment coords Xs and displacements Us:
    for (int j = 0; j < nsn; ++j) {
      col = nes*j;
      int IENrows = ISN[col + sgs] - 1; // Matlab numbering starts with 1
      col = nen*els;
      segmentNodesIDs[j] = IEN[col + IENrows] - 1; // Matlab numbering starts with 1
      for (int k = 0; k < nsd; ++k) {
	col = k*(int)(neq / nsd);
	Xs[k*nsn + j] = X[col + segmentNodesIDs[j]];
	Us[k*nsn + j] = U[col + segmentNodesIDs[j]];
      }
    }

    for (int g = 0; g < ngp; ++g) {

      // Gausspoint gap values and activeGPs:
      activeGPs[g] = (bool) activeGPsOld[i + g];

      if (!activeGPs[g]) {
	continue;
      }

      // master element index:
      col = (2 * nsd + 4)*GPs_len;
      const int elm = (int)GPs[col + i + g] - 1; // Matlab numbering starts with 1
      // master segment index:
      col = (2 * nsd + 5)*GPs_len;
      const int sgm = (int)GPs[col + i + g] - 1; // Matlab numbering starts with 1

      // master segment coords Xm and displacements Um:
      for (int j = 0; j < nsn; ++j) {
	col = nes*j;
	int IENrowm = ISN[col + sgm] - 1; // Matlab numbering starts with 1
	col = nen*elm;
	segmentNodesIDm[j] = IEN[col + IENrowm] - 1; // Matlab numbering starts with 1
	for (int k = 0; k < nsd; ++k) {
	  col = k*(int)(neq / nsd);
	  Xm[k*nsn + j] = X[col + segmentNodesIDm[j]];
	  Um[k*nsn + j] = U[col + segmentNodesIDm[j]];
	}
      }

      // Shape function and its derivatives of gausspoint's master segment

      double r = GPs[(nsd + 3)*GPs_len + i + g];
      double s = 0.0;

      if (npd == 2) {
	s = GPs[(nsd + 4)*GPs_len + i + g];
      }


      switch (nsn) {
      case 2:
	sfd2(Hm, dHm, r);
	break;
      case 6:
	sfd6(Hm, dHm, r, s);
	break;
      case 8:
	sfd8(Hm, dHm, r, s);
      }

      // evaluate gausspoint coords:
      for (int sdf = 0; sdf < nsd; ++sdf) {
	Xp[sdf] = 0.0;
	Xg[sdf] = 0.0;
	for (int k = 0; k < nsn; ++k) {
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
      for (int j = 0; j < nsn*npd; ++j) {
	dXs[j] = 0.0;
	dxs[j] = 0.0;
      }

      // Evaluate tangent vectors:
      for (int pdf = 0; pdf < npd; ++pdf) {
	for (int sdf = 0; sdf < nsd; ++sdf) {
	  for (int j = 0; j < nsn; ++j) {
	    col = j*npd*ngp;
	    const double dh = dH[col + g*npd + pdf];
	    dXs[npd*sdf + pdf] += dh*Xs[sdf*nsn + j];
	    dxs[npd*sdf + pdf] += dh*(Xs[sdf*nsn + j] + Us[sdf*nsn + j]);
	  }
	}
      }

      // Evaluate normal vector:
      double Normal[3];
      double normal[3];
      if (nsd == 2) {
	Normal[0] = dXs[1];
	Normal[1] = -dXs[0];
	Normal[2] = 0.0;

	normal[0] = dxs[1];
	normal[1] = -dxs[0];
	normal[2] = 0.0;
      }
      else if (nsd == 3) {
	const double dXs_dr1 = dXs[0];
	const double dXs_dr2 = dXs[npd + 0];
	const double dXs_dr3 = dXs[2 * npd + 0];
	const double dXs_ds1 = dXs[1];
	const double dXs_ds2 = dXs[npd + 1];
	const double dXs_ds3 = dXs[2 * npd + 1];
	Normal[0] = dXs_dr2*dXs_ds3 - dXs_dr3*dXs_ds2;
	Normal[1] = dXs_dr3*dXs_ds1 - dXs_dr1*dXs_ds3;
	Normal[2] = dXs_dr1*dXs_ds2 - dXs_dr2*dXs_ds1;

	const double dxs_dr1 = dxs[0];
	const double dxs_dr2 = dxs[npd + 0];
	const double dxs_dr3 = dxs[2 * npd + 0];
	const double dxs_ds1 = dxs[1];
	const double dxs_ds2 = dxs[npd + 1];
	const double dxs_ds3 = dxs[2 * npd + 1];
	normal[0] = dxs_dr2*dxs_ds3 - dxs_dr3*dxs_ds2;
	normal[1] = dxs_dr3*dxs_ds1 - dxs_dr1*dxs_ds3;
	normal[2] = dxs_dr1*dxs_ds2 - dxs_dr2*dxs_ds1;
      }

      const double jacobian      = sqrt(Normal[0] * Normal[0] + Normal[1] * Normal[1] + Normal[2] * Normal[2]);
      const double normal_length = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
      normal[0] /= normal_length;
      normal[1] /= normal_length;
      normal[2] /= normal_length;

      if (!keyContactDetection) {
	const double signGAPs = (Xp[0] - Xg[0])*normal[0] + (Xp[1] - Xg[1])*normal[1] + (Xp[2] - Xg[2])*normal[2];
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
      for (int j = 0; j < nsn; ++j) {
	const double hs = H[j*ngp + g];
	const double hm = Hm[j];

	for (int sdf = 0; sdf < nsd; ++sdf) {

	  C_s[j*nsd + sdf] += hs*normal[sdf];
	  C_m[j*nsd + sdf] -= hm*normal[sdf];
	  Gc[segmentNodesIDs[j] * nsd + sdf] -= epss*GAPs[g] * hs*normal[sdf] * gw[g] * jacobian;
	}
      }

      if(keyAssembleKc) {
	//if (i < 0.5*GPs_len) {
	for (int j = 0; j < nsn*nsd; ++j) { // loop over cols
	  for (int k = 0; k < nsn*nsd; ++k) { // loop over rows

	    const int jdof = j%nsd;
	    const int jnode = (j - jdof) / nsd; // row node

	    const int kdof = k%nsd;
	    const int knode = (k - kdof) / nsd; // col node

	    if(fabs(C_m[j] * C_m[k]) > 1e-50) {
	      cols[*len] = segmentNodesIDm[jnode] * nsd + jdof + 1;
	      rows[*len] = segmentNodesIDm[knode] * nsd + kdof + 1;
	      vals[*len] = 0.5 * C_m[j] * C_m[k] * gw[g] * jacobian;
	      (*len)++;
	      if (*len >= len_guess) printf("Error, len is too small: len = %i.\n", len_guess);
	    }

	    if(fabs(C_s[j] * C_m[k]) > 1e-50) {
	      cols[*len] = segmentNodesIDs[jnode] * nsd + jdof + 1;
	      rows[*len] = segmentNodesIDm[knode] * nsd + kdof + 1;
	      vals[*len] = 0.5 * C_s[j] * C_m[k] * gw[g] * jacobian;
	      (*len)++;
	      if (*len >= len_guess) printf("Error, len is too small: len = %i.\n", len_guess);
	    }

	    if(fabs(C_m[j] * C_s[k]) > 1e-50) {
	      cols[*len] = segmentNodesIDm[jnode] * nsd + jdof + 1;
	      rows[*len] = segmentNodesIDs[knode] * nsd + kdof + 1;
	      vals[*len] = 0.5 * C_m[j] * C_s[k] * gw[g] * jacobian;
	      (*len)++;
	      if (*len >= len_guess) printf("Error, len is too small: len = %i.\n", len_guess);
	    }

	    if(fabs(C_s[j] * C_s[k]) > 1e-50) {
	      cols[*len] = segmentNodesIDs[jnode] * nsd + jdof + 1;
	      rows[*len] = segmentNodesIDs[knode] * nsd + kdof + 1;
	      vals[*len] = 0.5 * C_s[j] * C_s[k] * gw[g] * jacobian;
	      (*len)++;
	      if (*len >= len_guess) printf("Error, len is too small: len = %i.\n", len_guess);
	    }
	  }
	}
      }
      // Fill C_m array by zeros:
      for (int j = 0; j < nsn*nsd; ++j) {
	C_s[j] = 0.0;
	C_m[j] = 0.0;
      }
      //}
    } // loop over gausspoints
  } // loop over GPs rows

  delete[] segmentNodesIDs;
  delete[] segmentNodesIDm;
  delete[] Xs;
  delete[] Xm;
  delete[] Us;
  delete[] Um;
  delete[] dXs;
  delete[] dxs;
  delete[] GAPs;
  delete[] activeGPs;
  delete[] C_s;
  delete[] C_m;
  delete[] Hm;
  delete[] dHm;
}

void getLongestEdgeAndGPs(double* longestEdge, double* GPs, int n, int nsd, int ngp, int neq, int nsn, int nes, int nen, int* elementID, int* segmentID, int* ISN, int* IEN, double* H, double* X) {
  // GPs legend:              Xg    els  sgs   gap  Xm    isActive elm sgm
  //double* GPs = new double[n*(nsd + 1 + 1  + 1  + nsd + 1       + 1 + 1)];
  int* segmentNodesID = new int[nsn];
  double* Xs = new double[nsn*nsd];
  double* Xg = new double[ngp*nsd];
  int g = 0;
  *longestEdge = 0.0;

  for (int e = 0; e < n; ++e) {
    int el = elementID[e] - 1;
    int sg = segmentID[e] - 1;

    // segment coords Xs:
    for (int i = 0; i < nsn; ++i) {
      const int IENrow = ISN[nes*i + sg] - 1; // Matlab numbering starts with 1
      segmentNodesID[i] = IEN[nen*el + IENrow] - 1; // Matlab numbering starts with 1
      for (int j = 0; j < nsd; ++j) {
	Xs[j*nsn + i] = X[j*(int)(neq / nsd) + segmentNodesID[i]];
      }
    }

    // evaluate gausspoint coords:
    for (int i = 0; i < ngp; ++i) {
      for (int sdf = 0; sdf < nsd; ++sdf) {
	Xg[i*nsd + sdf] = 0.0;
	for (int j = 0; j < nsn; ++j) {
	  Xg[i*nsd + sdf] += H[j*ngp + i] * Xs[sdf*nsn + j];
	}

	GPs[sdf*n*ngp + g] = Xg[i*nsd + sdf]; // slave gausspoint coords
	GPs[(nsd + 3 + sdf)*n*ngp + g] = 0.0; // init baricentric coords
      }
      GPs[nsd*n*ngp + g] = el + 1;            // slave element
      GPs[(nsd + 1)*n*ngp + g] = sg + 1;      // slave segment
      GPs[(nsd + 2)*n*ngp + g] = FLT_MAX;     // init gap
      GPs[(2 * nsd + 3)*n*ngp + g] = 0;       // init is NO active
      GPs[(2 * nsd + 4)*n*ngp + g] = 0;       // master element
      GPs[(2 * nsd + 5)*n*ngp + g] = 0;       // master segment

      g++;
    }

    // element edge length / diagonal:
    if (nsd == 2) {
      const double lengthOfEdge = sqrt(pow(Xs[0] - Xs[1], 2) + pow(Xs[2] - Xs[3], 2));
      *longestEdge = std::max(*longestEdge, lengthOfEdge);
    }
    else if (nsd == 3) {
      if (nsn == 6) {
	const double lengthOfEdge1 = sqrt(pow(Xs[0] - Xs[1], 2) + pow(Xs[3] - Xs[4], 2) + pow(Xs[6] - Xs[7], 2));
	const double lengthOfEdge2 = sqrt(pow(Xs[2] - Xs[1], 2) + pow(Xs[5] - Xs[4], 2) + pow(Xs[8] - Xs[7], 2));
	const double lengthOfEdge3 = sqrt(pow(Xs[0] - Xs[2], 2) + pow(Xs[3] - Xs[5], 2) + pow(Xs[6] - Xs[8], 2));
	*longestEdge = std::max(std::max(std::max(lengthOfEdge1, lengthOfEdge2), lengthOfEdge3), *longestEdge);
      }
      if (nsn == 8) {
	const double diagonal1 = sqrt(pow(Xs[2] - Xs[0], 2) + pow(Xs[10] - Xs[8], 2) + pow(Xs[18] - Xs[16], 2));
	const double diagonal2 = sqrt(pow(Xs[3] - Xs[1], 2) + pow(Xs[11] - Xs[9], 2) + pow(Xs[19] - Xs[17], 2));
	*longestEdge = std::max(std::max(diagonal1, diagonal2), *longestEdge);
      }
    }

  } // loop over elements

  delete[] segmentNodesID;
  delete[] Xs;
  delete[] Xg;

}

void getAABB(double* AABBmin, double* AABBmax, int nsd, int nnod, double* X, double longestEdge, int* IEN, int* ISN, int* elementID, int* segmentID, int n, int nsn, int nes, int nen, int neq) {

  int* segmentNodesID = new int[nsn];

  for (int sdf = 0; sdf < nsd; ++sdf) {
    AABBmin[sdf] = FLT_MAX;
    AABBmax[sdf] = -FLT_MAX;

    for (int e = 0; e < n; ++e) {
      int el = elementID[e] - 1;
      int sg = segmentID[e] - 1;

      // segment coords Xs:
      for (int i = 0; i < nsn; ++i) {
	const int IENrow = ISN[nes*i + sg] - 1; // Matlab numbering starts with 1
	segmentNodesID[i] = IEN[nen*el + IENrow] - 1; // Matlab numbering starts with 1
	for (int sdf = 0; sdf < nsd; ++sdf) {
	  const double x = X[sdf*(int)(neq / nsd) + segmentNodesID[i]];
	  AABBmin[sdf] = std::min(AABBmin[sdf], x);
	  AABBmax[sdf] = std::max(AABBmax[sdf], x);
	}
      }
    }

    for (int sdf = 0; sdf < nsd; ++sdf) {
      if ((AABBmax[sdf] - AABBmin[sdf]) < longestEdge) {
	AABBmax[sdf] += 0.5*longestEdge;
	AABBmin[sdf] -= 0.5*longestEdge;
      }
    }
  }
  delete[] segmentNodesID;
}

void evaluateContactConstraints(double* GPs, int* ISN, int* IEN, int* N, double* AABBmin, double* AABBmax, int* head, int* next, double* X, int* elementID, int* segmentID, int n, int nsn, int nsd, int npd, int ngp, int nen, int nes, int neq, double longestEdge) {

  int* segmentNodesID = new int[nsn];
  double* Xm = new double[nsn*nsd];
  double* Xmin = new double[nsd];
  double* Xmax = new double[nsd];
  double* Hm = new double[nsn];
  double* dHm = new double[nsn*npd];

  double Xt[9];
  double Xc[3];
  // If segment element is quad then it is divided to 4 triangles:
  int ntr = 1;
  if(nsn == 8) {
    ntr = 4;
  }

  for (int i = 0; i < n*ngp; ++i) {
    GPs[(nsd + 2)*n*ngp + i] = FLT_MAX;
  }

  for (int e = 0; e < n; ++e) {
    int el = elementID[e] - 1;
    int sg = segmentID[e] - 1;

    // segment coords Xm:
    for (int k = 0; k < nsd; ++k) {

      Xmin[k] = FLT_MAX;
      Xmax[k] = -FLT_MAX;
      Xc[k] = 0;

      for (int j = 0; j < nsn; ++j) {
	const int IENrow = ISN[nes*j + sg] - 1; // Matlab numbering starts with 1
	segmentNodesID[j] = IEN[nen*el + IENrow] - 1; // Matlab numbering starts with 1
	Xm[k*nsn + j] = X[k*(int)(neq / nsd) + segmentNodesID[j]];
	Xmin[k] = std::min(Xmin[k], Xm[k*nsn + j]);
	Xmax[k] = std::max(Xmax[k], Xm[k*nsn + j]);
	Xc[k] += Xm[k*nsn + j];
      }
      Xmin[k] -= 0.5*longestEdge;
      Xmax[k] += 0.5*longestEdge;
      Xc[k] /= nsn;
    }

    // Loop over segment triangles:
    for (int it = 0; it < ntr; ++it) {

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
	  Xt[1] = Xm[it+1];
	  Xt[2] = Xc[0];

	  Xt[3] = Xm[nsn+it];
	  Xt[4] = Xm[nsn+it+1];
	  Xt[5] = Xc[1];

	  Xt[6] = Xm[2*nsn+it];
	  Xt[7] = Xm[2*nsn+it+1];
	  Xt[8] = Xc[2];
	}

	// Min and Max of the triangle coords:
	for (int k = 0; k < nsd; ++k) {
	  Xmin[k] = FLT_MAX;
	  Xmax[k] = -FLT_MAX;
	  for (int j = 0; j < 3; ++j) {
	    Xmin[k] = std::min(Xmin[k], Xt[k*3 + j]);
	    Xmax[k] = std::max(Xmax[k], Xt[k*3 + j]);
	  }
	  Xmin[k] -= 0.5*longestEdge;
	  Xmax[k] += 0.5*longestEdge;
	}
      }

      int Imin[3];
      int Imax[3];
      double normal[3];
      double t1[3];
      double t2[3];
      double t3[3];
      double Xg[3];
      double Xp[3];

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

      const double normalLength = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
      normal[0] /= normalLength;
      normal[1] /= normalLength;
      normal[2] /= normalLength;

      for (int sdf = 0; sdf < nsd; ++sdf) {
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

      for (int i2 = Imin[2]; i2 <= Imax[2]; ++i2) {
	for (int i1 = Imin[1]; i1 <= Imax[1]; ++i1) {
	  for (int i0 = Imin[0]; i0 <= Imax[0]; ++i0) {
	    const int Ic = i2*N[0] * N[1] + i1*N[0] + i0;
	    int v = head[Ic];

	    while (v != -1) {
	      double d;
	      // Inside-outside algorithm:
	      bool isInside = false;
	      if (nsd == 2) {
		double r[2];
		d = 0.0;
		double t1_norm = 0.0;
		for (int i = 0; i < nsd; ++i) {
		  Xg[i] = GPs[i*n*ngp + v];
		  r[i] = Xg[i] - Xm[i*nsn + 0];
		  d += r[i] * t1[i];
		  t1_norm += t1[i] * t1[i];
		}
		d = sqrt(d);
		t1_norm = sqrt(t1_norm);
		// Check if inside edge1:
		if (d > 0.0 && d < t1_norm) {
		  isInside = true;
		  d = 0.0;
		  double sign = 0.0;
		  for (int i = 0; i < nsd; ++i) {
		    Xp[i] = Xm[i*nsn + 0] + t1[i] / t1_norm*r[i] * t1[i] / t1_norm;
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
		double r[9];
		double Q1[3];
		double Q2[3];
		double Q3[3];
		d = 0.0;
		for (int i = 0; i < nsd; ++i) {
		  Xg[i] = GPs[i*n*ngp + v];
		  r[i * 3 + 0] = Xg[i] - Xt[i*3 + 0];
		  r[i * 3 + 1] = Xg[i] - Xt[i*3 + 1];
		  r[i * 3 + 2] = Xg[i] - Xt[i*3 + 2];
		}

		// component:  X    Y    Z
		// r1:       r[0] r[3] r[6]
		// r2:       r[1] r[4] r[7]
		// r3:       r[2] r[5] r[8]
		Q1[0] = r[3] * t1[2] - r[6] * t1[1];
		Q1[1] = r[6] * t1[0] - r[0] * t1[2];
		Q1[2] = r[0] * t1[1] - r[3] * t1[0];

		Q2[0] = r[4] * t2[2] - r[7] * t2[1];
		Q2[1] = r[7] * t2[0] - r[1] * t2[2];
		Q2[2] = r[1] * t2[1] - r[4] * t2[0];

		Q3[0] = r[5] * t3[2] - r[8] * t3[1];
		Q3[1] = r[8] * t3[0] - r[2] * t3[2];
		Q3[2] = r[2] * t3[1] - r[5] * t3[0];

		const double Q1n = Q1[0] * normal[0] + Q1[1] * normal[1] + Q1[2] * normal[2];
		const double Q2n = Q2[0] * normal[0] + Q2[1] * normal[1] + Q2[2] * normal[2];
		const double Q3n = Q3[0] * normal[0] + Q3[1] * normal[1] + Q3[2] * normal[2];

		if (Q1n*Q2n > 0) {
		  if (Q1n*Q3n > 0) {
		    isInside = true;
		    d = r[0] * normal[0] + r[3] * normal[1] + r[6] * normal[2];
		    for (int i = 0; i < nsd; ++i) {
		      Xp[i] = Xg[i] - d*normal[i];
		    }
		  }
		}
	      } // else if (nsd = 3)

	      if (isInside) {
		// If distance is less then current closest distance:
		if (d < GPs[(nsd + 2)*n*ngp + v]) {

		  // Initial guess of the parametric coordinates on the triangle:
		  double r_len, r;
		  double s_len, s;
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
		    //printf("r0 = %f, s0 = %f\n", r, s);
		  }


		  // Local contact search by Least-square projection method:
		  double dr_norm;
		  int niter = 0;
		  int max_niter = 1000;
		  do {
		    switch (nsn) {
		    case 2:
		      sfd2(Hm, dHm, r);
		      break;
		    case 6:
		      sfd6(Hm, dHm, r, s);
		      break;
		    case 8:
		      sfd8(Hm, dHm, r, s);
		    }

		    double b1, b2, A11, A22, A12;
		    A11 = 0.0;
		    A22 = 0.0;
		    A12 = 0.0;
		    b1 = 0.0;
		    b2 = 0.0;
		    for (int sdf = 0; sdf < nsd; ++sdf) {
		      double x = 0.0;
		      double dx_dr = 0.0;
		      double dx_ds = 0.0;

		      for (int k = 0; k < nsn; ++k) {
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
		    }

		    double recDetA;
		    double invA11;
		    double invA22;
		    double invA12;
		    double dr;
		    double ds;

		    if (npd == 1) {
		      invA11 = 1 / A11;
		      dr = invA11*b1;
		      r += dr;
		      dr_norm = dr;
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
		      //printf("r = %f, s = %f, dr_norm = %f\n", r, s, dr_norm);
		    }

		    niter++;
		  } while (dr_norm > 1e-3 && niter < max_niter);

		  if (niter >= max_niter) {
		    printf("Fatal error: Local contact search do NOT converge.\n");
		  }

		  GPs[(nsd       + 2)*n*ngp + v] = d;      // store penetration
		  if(d < 0) {
		    GPs[(2 * nsd + 3)*n*ngp + v] = 1.0;    // set gausspoint to active state
		    GPs[(2 * nsd + 4)*n*ngp + v] = el + 1; // set master element
		    GPs[(2 * nsd + 5)*n*ngp + v] = sg + 1; // set master segment
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

  delete[] segmentNodesID;
  delete[] Xm;
  delete[] Xmin;
  delete[] Xmax;
  delete[] Hm;
  delete[] dHm;
}
