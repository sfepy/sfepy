/* dft.f -- translated by f2c (version 20050501).
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b2 = .33333333333333331;

/* Subroutine */ int getvxc_(real *n, real *vxc, integer *relat)
{
    /* Initialized data */

    static real thrd = .3333333333333333f;

    /* System generated locals */
    real r__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), atan(
	    doublereal), log(doublereal);

    /* Local variables */
    static real a, b, c__, q, r__, y, t1, t2, t3, y0, y2, t30, t31, t32, pi, 
	    mu, rs, aq2, yy0, dth, yyy, q2yb, ecca, fcca, beta, vcca, vxlda, 
	    dadbeta, drdbeta;

/* calculates V_LDA, from "n". */
/* relat=0,1, for relat = 0, it returns V_LDA, */
/* for relat=1 it returns V_RLDA */
/* adapted from VSCXC(D), by J. Vackar */
/* f2py intent(out) VXC */
/* f2py real*8 VXC */
/* f2py real*8 n */
    if (*n == 0.f) {
	*vxc = 0.f;
	return 0;
    }
    d__1 = (doublereal) dmax(0.f,*n);
    d__2 = (doublereal) thrd;
    dth = pow_dd(&d__1, &d__2);
/* CA C-term */
    a = .0621814f;
    b = 3.72744f;
    c__ = 12.93532f;
    y0 = -.10498f;
/*       Q=SQRT(4.*C-B*B) */
    q = 6.1520298314f;
/* -- */
    vxlda = dth * -.98474502184f;
    rs = .6203504909f / dth;
    y = sqrt(rs);
    y2 = rs;
    yyy = y2 + b * y + c__;
    yy0 = y0 * y0 + b * y0 + c__;
    q2yb = q / (y * 2.f + b);
    aq2 = atan(q2yb);
/*       AQ2=1./TAN(Q2YB) */
    t1 = log(y2 / yyy);
    t2 = b * 2.f / q * aq2;
    t30 = -b * y0 / yy0;
/* Computing 2nd power */
    r__1 = y - y0;
    t31 = log(r__1 * r__1 / yyy);
    t32 = (b + y0 * 2.f) * 2.f / q * aq2;
    t3 = t30 * (t31 + t32);
    ecca = a * .5f * (t1 + t2 + t3);
    fcca = a * .5f * ((c__ * (y - y0) - b * y0 * y) / ((y - y0) * yyy));
    vcca = ecca - thrd * fcca;
    if (*relat == 1) {
	c__ = 137.036f;
	pi = 3.1415926535897931f;
/* Computing 2nd power */
	r__1 = pi;
	d__1 = (doublereal) (r__1 * r__1 * 3 * *n);
	beta = 1.f / c__ * pow_dd(&d__1, &c_b2);
/* Computing 2nd power */
	r__1 = beta;
	mu = sqrt(r__1 * r__1 + 1);
/* Computing 2nd power */
	r__1 = beta;
	a = (beta * mu - log(beta + mu)) / (r__1 * r__1);
/* Computing 2nd power */
	r__1 = a;
	r__ = 1 - r__1 * r__1 * 1.5f;
	dadbeta = 2 / mu - a * 2 / beta;
	drdbeta = a * -3 * dadbeta;
	vxlda = vxlda * r__ + vxlda * drdbeta * beta / 4;
    }
    *vxc = vxlda + vcca;
    return 0;
} /* getvxc_ */

