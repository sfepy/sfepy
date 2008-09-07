#include <math.h>
#include "dft.h"

double max(double a, double b)
{
    if (a > b)
        return a;
    else
        return b;
}

double vxc(double n, int relat)
{
    double DTH,A,B,C,Y0,Q,VXLDA,RS,Y,Y2, YYY,YY0;
    double Q2YB;
    double AQ2,T1,T2,T30,T31,T32,T3,ECCA,FCCA,VCCA;
    double R,pi,beta,mu, dRdbeta, dAdbeta, c;
    double THRD=0.3333333333333333;

    if (n == 0) {
        return 0.0;
    }
    DTH=pow(max(0.,n),THRD);
    A=.0621814;
    B=3.72744;
    C=12.93532;
    Y0=-.10498;
    Q=6.1520298314;
    VXLDA=-.98474502184 * DTH;
    RS=.62035049090 / DTH;
    Y=sqrt(RS);
    Y2=RS;
    YYY=Y2+B*Y+C;
    YY0=Y0*Y0+B*Y0+C;
    Q2YB=Q/(2.*Y+B);
    AQ2=atan(Q2YB);
    T1=log(Y2/YYY);
    T2=(2.*B/Q) * AQ2;
    T30=-B*Y0/YY0;
    T31=log(pow((Y-Y0),2)/YYY);
    T32=(2.*(B+2.*Y0)/Q) * AQ2;
    T3=T30*(T31+T32);
    ECCA=.5*A*(T1+T2+T3);
    FCCA=.5*A*((C*(Y-Y0)-B*Y0*Y)/((Y-Y0)*YYY));
    VCCA=ECCA-THRD*FCCA;
    if (relat == 1) {
        c = 137.036;
        pi = 3.1415926535897931;
        beta = 1./c * pow((3*pow(pi,2)*n),(1./3));
        mu = sqrt(1+pow(beta,2));
        A = (beta*mu-log(beta+mu))/pow(beta,2);
        R = 1-3./2*pow(A,2);
        dAdbeta = 2/mu - 2*A/beta;
        dRdbeta = -3*A*dAdbeta;
        VXLDA = VXLDA*R+VXLDA*dRdbeta*beta/4;
    }
    return VXLDA+VCCA;
}
