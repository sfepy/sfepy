        subroutine integrate_radial_dirac_r(kappa,Z,E,nr,R,V,f,g)
        !integrates the Dirac eq, returns r*R, where R is the radial
        !solution.
Cf2py intent(out) f
Cf2py intent(out) g
        implicit none
        integer, intent(in) :: kappa
        integer, intent(in) :: Z
        real(8), intent(in) :: E
        integer, intent(in) :: nr
        real(8), intent(in) :: R(nr)
        real(8), intent(in) :: V(nr)
        real(8), intent(out) :: f(nr)
        real(8), intent(out) :: g(nr)

        real(8), parameter :: c=137.03599911

        real(8) y(2)
        integer i
        external derivsr

        real(8) beta
        integer nmax
        parameter (nmax=20000)
        common /dataderivs/ c_kappa,c_E,c_V,c_R,c_c,c_nr,c_i
        real(8) c_kappa,c_E,c_c,c_V(nmax),c_R(nmax)
        integer c_nr,c_i
        !copy some parameters to common, so that they can be used from derivs
        c_kappa=kappa
        c_E=E
        c_c=c
        c_nr=nr
        do i=1,nr
            c_R(i)=R(i)
            c_V(i)=V(i)
        enddo

        beta=sqrt(kappa**2-(Z/c)**2)
        g(1)=R(1)**(beta-1+1)
        f(1)=R(1)**(beta-2+1)*(beta+kappa)/((E-V(1))/c+2*c)

        !nn=0
        y(1)=g(1)
        y(2)=f(1)
        do i=2,nr
            c_i=i
            call rk4step(y,R(i-1),R(i)-R(i-1),derivsr)
            g(i)=y(1)
            f(i)=y(2)
            if ((abs(f(i)).ge.1e3).or.(abs(g(i)).ge.1e3)) return
            !if (g(i-1)*g(i).lt.0) nn=nn+1
        enddo
        return
        end subroutine

        subroutine derivsr(x,y,dydx)
Cf2py intent(inout) dydx
        implicit none
        real(8) x,y(2),dydx(2)
        integer nmax
        parameter (nmax=20000)
        common /dataderivs/ c_kappa,c_E,c_V,c_R,c_c,c_nr,c_i
        real(8) c_kappa,c_E,c_c,c_V(nmax),c_R(nmax)
        integer c_nr,c_i
        real(8) V
        external get_val
        call get_val(c_V,x,c_R,c_nr,V,c_i)
        !y=[g,f]
        !these equations can be found in:
        !http://ondrej.certik.cz/cookbook/thesis.ps
        !page 19, the [8] 2a and 2b equations
        dydx(1)=-(c_kappa)/x*y(1)+((c_E-V)/c_c+2*c_c)*y(2)
        dydx(2)=+(c_kappa)/x*y(2)-((c_E-V)/c_c      )*y(1)
        end

        subroutine get_val(f,x,R,n,V,i)
        implicit none
        real(8) f(n),x,R(n),V
        integer i
        integer j1,j2,n1,n2,n
        j1=i-1
        j2=j1+1

        n1=j1-1
        n2=j2+1 
        if (n1<1) then
          n2=n2-n1+1
          n1=1
        endif
        if (n2>n) then
          n1=n1-(n2-n)
          n2=n
        endif
        call interp(x,r(n1:n2),n2-n1+1,f(n1:n2),V)
        end

        subroutine rk4step(y,x,h,derivs)
        !advances the solution of
        !y_i'(x)=f_i(x,y_1,y_2,...,y_N)       for i=1,2,...,N; f_i are known
        !from x to x+h. 
        !parameters:
        !  y ... y(x) on the input and y(x+h) on the output.
        !  derivs(x,y,dydx) .... returns a vector of f_i(x,y1,y2,...)
cf2py intent(out) yout
        implicit none
        integer, parameter :: N=2
        real(8) h,x,dydx(N),y(N)
        EXTERNAL derivs
        real(8) dym(N),dyt(N),yt(N)
cf2py intent(inout) dydx
        call derivs(x,y,dydx)
        yt(:) = y(:) + h/2. * dydx(:)
        call derivs(x+h/2,yt,dyt)
        yt(:) = y(:) + h/2. * dyt(:)
        call derivs(x+h/2,yt,dym)
        yt(:) = y(:) + h * dym(:)
        dym(:) = dyt(:) + dym(:)
        call derivs(x+h,yt,dyt)
        y(:) = y(:) + h/6. * (dydx(:) + dyt(:) + 2. * dym(:))
        end

      subroutine integrate(x,f,n, S)
        implicit none
Cf2py intent(out) S
Cf2py real*8 S
Cf2py real*8 x
Cf2py real*8 f
      integer n
      real S, x(n), f(n), gr(n),cf(3,n)
      call fderiv(-1,n,x,f,gr,cf)
      S=gr(n)
      end

      subroutine normalize(K0,Y,R)
C     nanormuje Y na gridu R na jednicku, tj. int |Y|^2 dr = 1
Cf2py intent(in,out,copy) Y
Cf2py real*8 Y
Cf2py real*8 R
        implicit none
      integer K0
      real Y(K0),R(K0)
      real f(K0)
      real S
      f(:) = y(:)**2
      call integrate(r,f,k0, S)
      S = sqrt(abs(S))
      if (S.gt.0.d0) then
         y(:)=y(:)/S
      else
         write(*,*)
         write(*,'("Error(rdirac): zero wavefunction")')
         write(*,*)
         stop
      endif
      end

      subroutine WAVE4H_original (L,E,K0,U,Y,JM,AP,R,Z)
!
! L  (in)     ... orbital quantum number
! E  (in)     ... energy (Hartree units)
! K0 (in,out) ... max. grid index
! U  (array,in)     ... potential (Hartree units)
! Y  (array,out)    ... wave function
! JM, AP (in) ... grid parameters (R(j) = AP * j/(JM-j))
! R  (array,in)     ... radial grid (atomic units)
!
      IMPLICIT REAL*8 (A-H,O-Z)
      integer Z
C      IMPLICIT double precision (A-H,O-Z)
      DIMENSION U(K0),Y(K0),R(K0)
      DATA YJLIM/1.E+70/
C      DATA YJLIM/1.E+3/
!Cf2py intent(out) Y
C
C  HARTREE UNITS
C
C  STATEMENT FUNCTIONS
C     R(JX)=AP*JX/(JM-JX)
      EX2(RLN)=RLN*RLN
C      F(JX)=(SL/EX2(R(JX))+2.E+0*(U(JX)-E))/EX2(EX2(REAL(JM-JX)))
      F(JX)=(SL/EX2(R(JX))+2.d+0*(U(JX)-E))/EX2(EX2(REAL(JM-JX)))
C      A(JX)=1.E+0-G*F(JX)
      A(JX)=1.d+0-G*F(JX)
C

c      write (*,*) L
c      write (*,*) E
c      write (*,*) K0
c      write (*,*) U
c      write (*,*) Y
c      write (*,*) JM
c      write (*,*) AP
c      write (*,*) R

C  CONSTANT VALUES
C      G=AP*AP*JM*JM/12.E+0
      G=AP*AP*JM*JM/12.d+0
      SL=L*(L+1)
      K=K0
      IF (E.GE.0) THEN
        LF=1
        DO 5 LV=1,2*L+1,2
    5   LF=LF*LV
        CNR=1.d+0/LF*SQRT(2.*E)**L
      ELSE
        CNR=1.d+0
      END IF
C
C  --------------------------------/2.
C
      CNR=1.d+0
      YJ2=CNR*(1.d+0/JM)*R(1)**L
      Y(1)=YJ2
      YJ1=CNR*(2.d+0/JM)*R(2)**L
      !added by Ondrej Certik to fix the asymptotic
      if (L.eq.0) YJ1=YJ1-Z*YJ1*(R(2)-R(1))
      Y(2)=YJ1
      AJ2=A(1)
      AJ1=A(2)
C
      J=3
   10	AJ=A(J)
	YJ=((12.d+0-10.d+0*AJ1)*YJ1-AJ2*YJ2)/AJ
C
      IF (J.LT.K .AND. ABS(YJ).LT.YJLIM) THEN
	Y(J)=YJ
	YJ2=YJ1
	YJ1=YJ
	AJ2=AJ1
	AJ1=AJ
	J=J+1
	GOTO 10
      END IF
C
      Y(J)=YJ
C
      KX=J
      DO 20 J=1,KX
   20	Y(J)=Y(J)*JM/J
      IF (KX.LT.K0) THEN
        DO 31 J=KX+1,K0
   31     Y(J)=0.d+0
      END IF
      K0=KX
      RETURN
C
      END

      subroutine parsefunction(k0,y,nodes,minidx,positive)
C parses the function y, returns:
C nodes: the number of intersection with the x-axis, not counting
C   the beginning and infinity
C minidx: the index of the last minimum, i.e. the place, where the
C   function was last closest to the x-axis (in absolute value), so for
C   indexes > minidx, the function goes to +- infinity, or possibly to
C   zero in infinity
C positive: true or false, depending if the y is approaching the x-axis
C   in the infinity from above (positive) or below (negative).
C This information is is used in the routine checke to determine, if the
C function lies below or above a certain energy.
      implicit none
      integer k0,nodes,minidx
      logical positive
      real y(k0)

      integer last_sign,last_j,maxj,j,isy,k
      nodes=0
      last_sign=sign(1.0,y(0+1))
      last_j=-1

      maxj=k0-1
      do j=1,k0-1
        isy=sign(1.0,y(j+1))
        if (isy .eq. -last_sign) then
            last_sign=isy
            last_j=j
            nodes=nodes+1
        endif
        if (abs(y(j+1)) .gt. 990) then
C        if (abs(y(j+1)) .gt. 1e50) then
            maxj=j
            goto 21
        endif
      enddo
21    continue

      k=maxj
      do while (abs(y(k-1+1)) .lt. abs(y(k+1))) 
        k=k-1
        if (k.eq.0) goto 22
      enddo
22    continue

C      write (*,*) nodes,k,last_j,k.le.last_j
      if (k.le.last_j) nodes=nodes-1
C      write (*,*) nodes,k,last_j,k.le.last_j

      minidx=k
      positive = y(maxj+1) .gt. 0
C      if (minidx.eq.0) minidx=maxj
C     this is because we are normalizing the wavefunction later on
C     so that we had something to normalize. (minidx=0 causes y=0 for
C     all r)
      if (minidx.eq.0) minidx=6
      end

      subroutine checke(n,l,E,k0,u,r,Z,relat, y,Eabove,minidx)
C Checks the energy E: calls integrate_radial_dirac with the energy E and
C returns Eabove=true (false) if E is above (below) the eigen value.
C Also returns minidx, see parsefunction.
      !according to the parameter relat, it will call the correct
      !integration routine.
      implicit none
      integer l,minidx,n,relat,Z,kappa,k0
      real u(k0),y(k0),f(k0),g(k0),r(k0),e
      logical Eabove

      integer nods_theory,nods_actual
      logical positive_for_big_r

      integer jm
      real ap
      integer k0o

      if (relat.eq.-1) then
          AP=1.0/(1/R(1)-2/R(2))
          JM=int((R(1)-R(2))/(R(1)-R(2)/2)+0.5)
          k0o=K0
          call WAVE4H_original (L,E,k0o,U,Y,JM,AP,R,Z)
          y(:) = y(:)*R(:)
          call parsefunction(k0,Y,nods_actual,minidx,positive_for_big_r)
      else if (relat.eq.0) then
          call integrate_radial_schrodinger(l,Z,E,k0,R,U,Y)
          call parsefunction(k0,Y,nods_actual,minidx,positive_for_big_r)
      else if (relat.eq.1) then
          stop "Scalar relativistic case not implemented yet"
      else if (relat.eq.2 .or. relat.eq.3) then
          if (relat.eq.3) then
              if (l.eq.0) then
                    stop "for l=0 only spin up (relat=2) is allowed"
              endif
              kappa=l
          else
              kappa=-l-1
          endif
          call integrate_radial_dirac_r(kappa,Z,E,k0,R,U,f,g)
          call parsefunction(k0,g,nods_actual,minidx,positive_for_big_r)
          y(:) = sqrt(f(:)**2+g(:)**2)
      else
          stop "wrong value of relat."
      endif

      nods_theory=n-l-1

      if (nods_actual .eq. nods_theory) then
          if (mod(nods_theory,2) .eq. 0) then 
              Eabove=.not. positive_for_big_r
          else
              Eabove=positive_for_big_r
          endif
      else
          Eabove=nods_actual .gt. nods_theory
      endif

C      write (*,*) E,nods_theory,nods_actual,positive_for_big_r,minidx,
C     *    Eabove
C      if ((nods_actual.eq.5).and.(minidx.eq.497)) write (*,*) Y
      end

      subroutine solve_radial_eigenproblem(n,l,Ein,EII,eps,k0,u,r,e,
     *    Z,relat, y)
C Solves the radial Dirac (Schrodinger) equation and returns the eigen value
C (E) and normalized eigenvector (Y) for the given "n" and "l".
C
C    Finds the wavefunction with defined "n" and "l". The potential is "u".
C    rel ... 0 nonrelat (runge-kutta)
C            1 semirelat (runge-kutta)
C            2 relat (runge-kutta) spin up
C            3 relat (runge-kutta) spin down
C           -1 nonrelat (polynomial approximation)
C    Ein,EII gives the window, where to look at the solution by halving the
C      interval. if there is no solution in the interval, the interval is
C      automatically adjusted, so you don't have to worry about it. The
C      calculation is a little faster, if you set Ein and EII near the correct
C      energy, but it's not necessary.
C    eps ... the solver halves the difference between Emin and Emax until
C            |Emax-Emin|<eps
C    r ... the grid on which to solve the equation, must be hyperbolic (the
C        solver counts with it)
C    Z is the coulomb term in the potential V=-Z/r+..., it is used for
C    the asymptotic
C
C    returns a wavefunction (i.e., not y*r, but y)
Cf2py intent(out) Y
Cf2py intent(out) E
Cf2py real*8 Ein
Cf2py real*8 EII
Cf2py real*8 eps
Cf2py real*8 u
Cf2py real*8 r
Cf2py real*8 e
Cf2py real*8 y
        implicit none
      integer n,l,k0,relat,Z
      real u(k0),y(k0),r(k0),e,eps,Ein,EII
C      integer i
      if (.not.(n.gt.0)) stop "n>0 not satisfied"
      if (.not.((0.le. l).and.(l.lt. n))) then
         stop "0<=l<n not satisfied"
      endif
      call finde(n,l,Ein,EII,eps,k0,u,r,e,Z,relat,y)
Cfinde returns y = r*R, and normalize will just norm int y^2 = int r^2 R^2=1,
Cwhich is correct
      call normalize(k0,y,r)
      !at the end, calculate R=y/r and return R.
      y(:) = y(:)/r(:)
      end

      subroutine finde(n,l,Ein,EII,eps,k0,u,r,e,Z,relat, y)
C Helper function, which returns unnormalized solution to the eigen
C problem. It's called from solve_radial_eigenproblem.
        implicit none
      integer n,l,k0,relat,Z
      real u(k0),y(k0),r(k0),e,eps,Ein,EII,EIIc

      real Emin,Emax
      integer counter,minidx,j
      logical isbig

C      write (*,*) "counter MAX:",Ein,EII,eps

C     rewrite this fucking loop

      if (EII.le.1e-10) EII=1e-10
      Emax=Ein+EII
      call checke(n,l,Emax,k0,u,r,Z,relat, y,isbig,minidx)
      counter=0
      EIIc=EII
      do while (.not. isbig)
          counter=counter+1
          Emin=Emax
          EIIc=EIIc*2
          Emax=Emax+EIIc
C          write (*,*) "counter MAX:",counter,Emin,Emax,EIIc
          if (Emax.ge.1e6) then
              stop "finde: Emax>1e6"
          endif
          call checke(n,l,Emax,k0,u,r,Z,relat, y,isbig,minidx)
      enddo

      Emin=Ein-EII
      call checke(n,l,Emin,k0,u,r,Z,relat, y,isbig,minidx)
      counter=0
      do while (isbig)
          counter=counter+1
          Emax=Emin
          Emin=Emin-EII*(2**counter)
C          write (*,*) "counter MIN:",counter
          if (Emin.le.-1e6) then
              stop "finde: Emin<-1e6"
          endif
          call checke(n,l,Emin,k0,u,r,Z,relat, y,isbig,minidx)
      enddo

C      Emax=0.0
C      Emin=-3000.0


C      write (*,*) "======"
      counter=0
      do while (counter .lt. 1000)
        counter=counter+1
C        write (*,*) Emin,Emax,isbig,minidx
        E=(Emin+Emax)/2.0
        if (Emax-Emin<eps) goto 10
        call checke(n,l,E,k0,u,r,Z,relat, y,isbig,minidx)
        if (isbig) then
            Emax=E
        else 
            Emin=E
        endif
      enddo
10    continue
      if (minidx.lt.1) then
          stop "finde: minidx<1, something went wrong"
      endif
      do j=minidx,k0-1
        y(j+1)=0.0
      enddo
C      k0=minidx
C      if (minidx.eq.k0-1) then
C          write (*,*) "minidx",minidx,n,l
C          write (*,*) y
C      endif
      end

      subroutine interp(t,X,n,Y,val)
        implicit none
      real t
      integer n
      real X(n),Y(n)
      real val,f,denum
      integer i,j
      val=0.
      do j=1,n
        f=1.
        denum=1
        do 30 i=1,n
          if (i.eq.j) goto 30
          f=f*(t-X(i))
          denum=denum*(X(j)-X(i))
30      continue
        val=val+Y(j)*f/denum
      enddo
      end

        subroutine spline(n,x,ld,f,cf)
        implicit none
        ! arguments
        integer, intent(in) :: n
        real(8), intent(in) :: x(n)
        integer, intent(in) :: ld
        real(8), intent(in) :: f(ld,n)
        real(8), intent(out) :: cf(3,n)
        ! local variables
        integer i
        real(8) t1,t2
        ! automatic arrays
        real(8) w(n)
        do i=1,n-1
          cf(3,i)=1.d0/(x(i+1)-x(i))
          cf(1,i)=cf(3,i)*(f(1,i+1)-f(1,i))
        end do
        t1=0.5d0/(x(3)-x(1))
        cf(2,1)=t1*(cf(1,2)-cf(1,1))
        do i=2,n-1
          t2=x(i)-x(i-1)
          w(i)=t2*t1
          t1=1.d0/(2.d0*(x(i+1)-x(i-1))-t2*w(i))
          cf(2,i)=t1*(cf(1,i)-cf(1,i-1)-t2*cf(2,i-1))
        end do
        ! back-substitution
        do i=n-2,1,-1
          cf(2,i)=cf(2,i)-w(i+1)*cf(2,i+1)
        end do
        ! determine coefficients
        cf(2,n)=0.d0
        do i=1,n-1
          t2=cf(2,i+1)-cf(2,i)
          cf(3,i)=t2*cf(3,i)
          cf(2,i)=3.d0*cf(2,i)
          cf(1,i)=cf(1,i)-(cf(2,i)+t2)*(x(i+1)-x(i))
        end do
        ! end-point coefficients for extrapolation
        t1=x(n)-x(n-1)
        cf(1,n)=(3.d0*cf(3,n-1)*t1+2.d0*cf(2,n-1))*t1+cf(1,n-1)
        cf(2,n)=6.d0*cf(3,n-1)*t1+2.d0*cf(2,n-1)
        cf(3,n)=cf(3,n-1)
        return
        end subroutine

        subroutine fderiv(m,n,x,f,g,cf)
        implicit none
        ! arguments
        integer, intent(in) :: m
        integer, intent(in) :: n
        real(8), intent(in) :: x(n)
        real(8), intent(in) :: f(n)
        real(8), intent(out) :: g(n)
        real(8), intent(out) :: cf(3,n)
        ! local variables
        integer i
        real(8) dx
        if (n.le.0) then
          write(*,*)
          write(*,'("Error(fderiv): invalid number of points : ",I8)') n
          write(*,*)
          stop
        end if
        if (m.eq.0) then
          g(:)=f(:)
          return
        end if
        if (m.ge.4) then
          g(:)=0.d0
          return
        end if
        ! high accuracy (anti-)derivatives from a clamped spline fit to the data
        call spline(n,x,1,f,cf)
        select case(m)
        case(:-1)
          g(1)=0.d0
          do i=1,n-1
            dx=x(i+1)-x(i)
            g(i+1)=g(i)+(((0.25d0*cf(3,i)*dx+
     *          0.3333333333333333333d0*cf(2,i))*dx 
     *          +0.5d0*cf(1,i))*dx+f(i))*dx
          end do
        case(1)
          g(:)=cf(1,:)
        case(2)
          g(:)=2.d0*cf(2,:)
        case(3)
          g(:)=6.d0*cf(3,:)
        end select
        return
        end subroutine

        subroutine integrate_radial_poisson(density,R,nr,V0,V0d,V)
        !solves V''(r) + 2/r*V'(r) + V(r) = density
        !with initial conditions V(R(1)) = V0 and V'(R(1)) = V0d
        !using 4th order runge-kutta method.
        implicit none
        integer, intent(in) :: nr
        real(8), intent(in) :: density(nr)
        real(8), intent(in) :: R(nr)
        real(8), intent(in) :: V0
        real(8), intent(in) :: V0d
        real(8), intent(out) :: V(nr)

        real(8) y(2)
        integer i
        external derivspoisson

        integer nmax
        parameter (nmax=20000)
        common /dataderivspoisson/ c_density,c_R,c_nr,c_i
        real(8) c_density(nmax),c_R(nmax)
        integer c_nr,c_i
        c_nr=nr
        do i=1,nr
            c_R(i)=R(i)
            c_density(i)=density(i)
        enddo

        y(1)=V0
        y(2)=V0d

        V(1)=y(1)
        do i=2,nr
            c_i=i
            call rk4step(y,R(i-1),R(i)-R(i-1),derivspoisson)
            V(i)=y(1)
        enddo
C        V(nr)=y(1)
C       do j=1,nr-1
C           i=nr-j
C           c_i=i
C           call rk4step(y,R(i+1),R(i)-R(i+1),derivspoisson)
C           V(i)=y(1)
C       enddo
        return
        end subroutine

        subroutine derivspoisson(x,y,dydx)
        implicit none
        real(8) x,y(2),dydx(2)
        integer nmax
        parameter (nmax=20000)
        common /dataderivspoisson/ c_density,c_R,c_nr,c_i
        real(8) c_density(nmax),c_R(nmax)
        integer c_nr,c_i
        real(8) dns
        external get_val
        call get_val(c_density,x,c_R,c_nr,dns,c_i)
        dydx(1)=y(2)
        dydx(2)=-dns-2.0*y(2)/x
C        write (*,*) dns,2.0*y(2)/x
        end

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

        subroutine integrate_radial_schrodinger(l,Z,E,nr,R,V,g)
        !integrates the Schrodinger eq., returns r*R, where R is the radial
        !solution.
Cf2py intent(out) g
        implicit none
        integer, intent(in) :: l
        integer, intent(in) :: Z
        real(8), intent(in) :: E
        integer, intent(in) :: nr
        real(8), intent(in) :: R(nr)
        real(8), intent(in) :: V(nr)
        real(8), intent(out) :: g(nr)

        real(8), parameter :: c=137.03599911

        real(8) y(2)
        integer i
        external derivs_schrodinger

        integer nmax
        parameter (nmax=20000)
        common /dataderivs/ c_kappa,c_E,c_V,c_R,c_c,c_nr,c_i
        real(8) c_kappa,c_E,c_c,c_V(nmax),c_R(nmax)
        integer c_nr,c_i
        !copy some parameters to common, so that they can be used from derivs
        c_kappa=l
        c_E=E
        c_c=c
        c_nr=nr
        do i=1,nr
            c_R(i)=R(i)
            c_V(i)=V(i)
        enddo

        if (l.eq.0) then
            y(1)=1-Z*R(1)
            y(2)=-Z
        else
            y(1)=R(1)**l
            y(2)=l*R(1)**(l-1)
        endif

        g(1)=y(1)*R(1)
        do i=2,nr
            c_i=i
            call rk4step(y,R(i-1),R(i)-R(i-1),derivs_schrodinger)
            g(i)=y(1)*R(i)
            if (abs(g(i)).ge.1e3) return
        enddo
        end subroutine

        subroutine derivs_schrodinger(x,y,dydx)
Cf2py intent(inout) dydx
        implicit none
        real(8) x,y(2),dydx(2)
        integer nmax
        parameter (nmax=20000)
        common /dataderivs/ c_kappa,c_E,c_V,c_R,c_c,c_nr,c_i
        real(8) c_kappa,c_E,c_c,c_V(nmax),c_R(nmax)
        integer c_nr,c_i
        real(8) V
        external get_val
        call get_val(c_V,x,c_R,c_nr,V,c_i)
        !y=[g,f]
        dydx(1)=y(2)
        dydx(2)=-2.0/x*y(2)+(2*(V-c_E)+c_kappa*(c_kappa+1)/x**2)*y(1)
        end
