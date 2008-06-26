      subroutine getvxc(n,VXC,relat)
      !calculates V_LDA, from "n". 
      !relat=0,1, for relat = 0, it returns V_LDA, 
      !for relat=1 it returns V_RLDA
      !adapted from VSCXC(D), by J. Vackar
Cf2py intent(out) VXC
Cf2py real*8 VXC
Cf2py real*8 n
      implicit none
      real n,VXC
      real DTH,A,B,C,Y0,Q,VXLDA,RS,Y,Y2, YYY,YY0
      real Q2YB
      real AQ2,T1,T2,T30,T31,T32,T3,ECCA,FCCA,VCCA
      real R,pi,beta,mu,THRD, dRdbeta, dAdbeta
      DATA THRD/0.3333333333333333/
      integer relat

      if (n.eq.0) then
          VXC=0
          return
      endif
      DTH=MAX(0.,n)**THRD
C CA C-term
        A=.0621814
        B=3.72744
        C=12.93532
        Y0=-.10498
c       Q=SQRT(4.*C-B*B)
        Q=6.1520298314
C --
        VXLDA=-.98474502184 * DTH
        RS=.62035049090 / DTH
        Y=SQRT(RS)
        Y2=RS
        YYY=Y2+B*Y+C
        YY0=Y0*Y0+B*Y0+C
        Q2YB=Q/(2.*Y+B)
        AQ2=ATAN(Q2YB)
c       AQ2=1./TAN(Q2YB)
        T1=LOG(Y2/YYY)
        T2=(2.*B/Q) * AQ2
        T30=-B*Y0/YY0
        T31=LOG((Y-Y0)**2/YYY)
        T32=(2.*(B+2.*Y0)/Q) * AQ2
        T3=T30*(T31+T32)
        ECCA=.5*A*(T1+T2+T3)
        FCCA=.5*A*((C*(Y-Y0)-B*Y0*Y)/((Y-Y0)*YYY))
        VCCA=ECCA-THRD*FCCA
        if (relat.eq.1) then
            c = 137.036
            pi = 3.1415926535897931
            beta = 1./c * (3*pi**2*n)**(1./3)
            mu = sqrt(1+beta**2)
            A = (beta*mu-log(beta+mu))/beta**2
            R = 1-3./2*A**2
            dAdbeta = 2/mu - 2*A/beta
            dRdbeta = -3*A*dAdbeta
            VXLDA = VXLDA*R+VXLDA*dRdbeta*beta/4
        endif
        VXC=VXLDA+VCCA
      END
