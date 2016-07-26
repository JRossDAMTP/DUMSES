module stratified
  use amr_parameters
  implicit none

  logical :: floor,smooth,addmass,massDiff,stratif
  real(dp) :: zfloor,dfloor,d0,mass0
  DOUBLE PRECISION :: rho_cool, c2_cool, t_cool

end module stratified
!----------------------------------------------------------------
!setup the initial conditions for a stratified disk atmosphere
!----------------------------------------------------------------
subroutine condinit(mype)
  use variables
  use hydro_parameters
  use stratified
  use const
  implicit none
#if WITHMPI==1
#include "mpif.h"
#endif

  integer :: mype

  integer  :: i,j,k,l,iseed
  real(dp) :: H,beta,B0,p0,amp,rhofloor,polytrop, pi
  real(dp) :: rho,rhovx,rhovy,rhovz,Bxc,Byc,Bzc,E,csq,P,P_mid
  real :: rvalue
  !
  ! type :
  ! 'r': random, 'h': henrik's r-modes
  character(len=10) :: type
  !character(len=10) :: nchar,Rm,Re !mode number (character)
  !character(len=256) :: filename,dirname
  !integer :: nz_eigen
  !real(dp) :: interpolated, P_mid, k_eigen, omega_eigen
  !real(dp), dimension(:), allocatable :: z_eigen,ux_eigen,uy_eigen,uz_eigen,rho_eigen,P_eigen
  DOUBLE PRECISION, DIMENSION(1:3) :: y0
  !
  namelist /init_params/d0,beta,gamma,polytrop,amp,type,floor,zfloor,dfloor,smooth,addmass,massDiff,stratif,rho_cool, c2_cool, t_cool,P_mid
  !
  ! Calculate the initial state for uin...
  ! uin   (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)
  ! uin = (\rho, \rho vx, \rho vy, \rho vz, Etot, bx, by, bz)

  pi=2.d0*asin(1.d0)

  gamma = 1.4
  d0=1.d0
  beta=100.d0
  amp=1.d-1
  type='r'
  floor=.false. ; smooth=.false. ; addmass=.false. ; massDiff=.false. ; stratif = .true.
  zfloor=5.
  ! parameters of the cooling function :
  rho_cool=1.d-2  ! density cutoff
  c2_cool=0.25     ! equilibrium temperature
  t_cool = 1.d0  ! cooling time, if t_cool negative then there is no cooling
  P_mid=16/9
  open(unit=1,file='input' ,status='old')
  read(1,init_params)
  close(1)

 
  H=ciso/Omega0
  B0=0.d0

 ! WRITE(*,*) 'H :', H, '   P_mid :', P_mid

 
  uin=0.d0
  iseed=mype

  ! Compute c2_cool such c2 = 1. at midplane correspond to a thermal equilibrium :


  !Set vertical profile
  do l=1,ngrid
     do k=ku1,ku2
        ! Compute stationary state :

        ! store the sound speed squared c2 in the energy variable (it will updated later taking into account the velocity and magnetic field


           !Momentum (Random condition)
           do j=1,ny
              do i=1,nx 
		 uin(l,i,j,k,1)=1.d0

                 uin(l,i,j,k,2) = 0.d0

                 uin(l,i,j,k,4) = 0.d0
 
                 uin(l,i,j,k,3) = 0.d0
              end do
           end do
 
        
     end do
  end do

  ! magnetic field :
  uin(1,iu1:iu2,ju1:ju2,ku1:ku2,6:8) = 0.d0
  uin(1,iu1:iu2,ju1:ju2,ku1:ku2,nvar+1:nvar+3) = 0.d0
 

  !Set total energy density
#if ISO==0
  do k=ku1,ku2-1
     do j=ju1,ju2-1
        do i=iu1,iu2-1
           uin(1,i,j,k,5)=P_mid/(0.4d0)
        end do
     end do
  end do
#endif


  !Calculate initial mass in computational box
  mass0=0.d0
  do k=1,nz
     do j=1,ny
        do i=1,nx
           mass0 = mass0 + uin(1,i,j,k,1)
        end do
     end do
  end do
#if WITHMPI==1
  call sumBcast(mass0,1)
#endif
  mass0=mass0*dx*dy*dz

  return
end subroutine condinit
!====================================================================
!  numerical recipes random number generator ran2
!    requires input seed value=iseed
!    returns real random number=rvalue
!    Also updates iseed for next call 
!
      subroutine ran2(iseed,rvalue)
     
      integer iseed
      real rvalue
      integer idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real AM,EPS,RNMX
      parameter (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&
               & IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,&
               & IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      integer idum2,jj,kk,iv(NTAB),iy
      data idum2/123456789/, iv/NTAB*0/, iy/0/
!
      idum=iseed
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 jj=NTAB+8,1,-1
          kk=idum/IQ1
          idum=IA1*(idum-kk*IQ1)-kk*IR1
          if (idum.lt.0) idum=idum+IM1
          if (jj.le.NTAB) iv(jj)=idum
11      continue
        iy=iv(1)
      endif
      kk=idum/IQ1
      idum=IA1*(idum-kk*IQ1)-kk*IR1
      if (idum.lt.0) idum=idum+IM1
      kk=idum2/IQ2
      idum2=IA2*(idum2-kk*IQ2)-kk*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      jj=1+iy/NDIV
      iy=iv(jj)-idum2
      iv(jj)=idum
      if(iy.lt.1)iy=iy+IMM1
      rvalue=min(AM*iy,RNMX)
      iseed=idum
      return
      end
!====================================================================
!====================================================================
!====================================================================
subroutine GetInterpolated(zloc,floc,z,func,n)
  use amr_parameters
  implicit none
  integer :: n
  real(dp) :: zloc,floc
  real(dp), dimension(n) :: z,func

  integer :: k
  real(dp) :: z0,f0,df,dz

  IF (verbose) WRITE(*,*) 'entering GetInterpolated...'
  
  k=1
  do while (z(k) .GT. zloc)
     k=k+1
  end do
  f0=func(k) ; z0=z(k)
  df=func(k)-func(k-1)
  dz=z(k)-z(k-1)
  floc=f0+(zloc-z0)*df/dz

  return
end subroutine GetInterpolated
 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE DERIVE0(z0,y0,dzy0)
  ! This subroutine gives the differential system governing the vertical structure
  ! It returns in dzy0 the derivative of y0 as a function z_loc
  ! y0(1) : pressure P
  ! y0(2) : P/rho
  ! y0(3) : vertical energy flux (from thermal diffusion)
  USE VARIABLES  
  use hydro_parameters
  use stratified
  IMPLICIT NONE
  ! ------------- input and output variables ----------------------
  DOUBLE PRECISION, INTENT(IN)   :: z0
  DOUBLE PRECISION, DIMENSION(1:3), INTENT(IN)   :: y0
  DOUBLE PRECISION, DIMENSION(1:3), INTENT(OUT) :: dzy0
  !-------------- local variables--------------------------------
  DOUBLE PRECISION :: L_cool, rho, c2, kappa_loc
  !----------------------------------------------------------	


  !IF (VERBOSE) WRITE(*,*) 'y0 :', y0
  rho = y0(1)/y0(2)
  c2 = gamma*y0(2)
  CALL GET_KAPPA(kappa_loc,rho,y0(1))
 
  dzy0(1) = - rho*Omega0**2*z0
  dzy0(2) = - (gamma-1.d0)/(gamma*kappa_loc)*y0(3)/rho
  !IF (VERBOSE) WRITE(*,*) 'dzy0(2) =',dzy0(2), - (gamma-1.d0)/(gamma*kappa)*y0(3)*y0(2)/y0(1),  - (gamma-1.d0)/(gamma*kappa), y0(3)*y0(2)/y0(1), y0(3), y0(2), y0(1)
  CALL GET_COOL(L_cool,rho,c2)
  dzy0(3) = rho*nu*(1.5d0*Omega0)**2 + L_cool

  !IF (VERBOSE) WRITE(*,*) 'In derive, z :', z0, '  y0 :', y0, '  dzy0 :', dzy0, '  gamma :', gamma, '  kappa :', kappa, '  nu :', nu, '  Lcool :', L_cool

  RETURN 
END SUBROUTINE DERIVE0
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE INTEGRE0(y0)
  ! Integrate the equations giving var from the midplane to far far away.. and give the boundary condition there (c2 = c2_cool)
  USE VARIABLES  
  use hydro_parameters
  use stratified
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(1:3), INTENT(out) :: y0
  DOUBLE PRECISION :: z0,dz0, zfar
  INTEGER          :: ncount, ncountmax, naccur, i
  !DOUBLE PRECISION :: x0, v0, c0, M2, M, dvx, dphidxsh
 
  naccur = 10000
  ncount = 0
  ncountmax = 0
  zfar = 4.d0


  ! ++++++++++++++++++++++++++++++++ Initialisation at the midplane +++++++++++++++++++++++++++++++
  z0 = 0.d0
  y0(1) = 1./gamma 
  y0(2) = 1./gamma
  y0(3) = 0.d0
  !WRITE(*,*) 'y:', y
  
  ! ++++++++++++++++++++++++++++++++ Integration to far far away ++++++++++++++++++++++++++++++
  dz0 = zfar/naccur
  DO WHILE( (z0 .LE. zfar).AND.(gamma*y0(2) .GT. c2_cool).AND.(y0(3).GE.0.d0) )
     CALL RKUTTA0(z0,y0,dz0,ncount)
     z0=z0+dz0
     IF (ncount .GT. ncountmax) THEN
        ncountmax = ncount
     ENDIF
     !IF (VERBOSE) WRITE(*,*) 'z :', z0, '  y :', y0
     !IF (VERBOSE) WRITE(*,*) 'rho :', y0(1)/y0(2), '  c2 :', gamma*y0(2)
  END DO
  
  ! Boundary condition  :
  WRITE(*,*) 'z =', z0,' y0 :', y0
  IF (VERBOSE) WRITE(*,*) 'Runge Kutta a-t-il convergé? ncountmax :', ncountmax


  RETURN
END SUBROUTINE INTEGRE0
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE INTEGRE_z(zfar,y0)
  ! Integrate the equations giving var from the midplane to far far away.. and give the boundary condition there (c2 = c2_cool)
  USE VARIABLES  
  use hydro_parameters
  use stratified
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(1:3), INTENT(OUT) :: y0
  DOUBLE PRECISION, INTENT(IN) :: zfar
  DOUBLE PRECISION :: z0,dz0
  INTEGER          :: ncount, ncountmax, naccur, i
  !DOUBLE PRECISION :: x0, v0, c0, M2, M, dvx, dphidxsh
 
  naccur = 10000
  ncount = 0
  ncountmax = 0

  ! ++++++++++++++++++++++++++++++++ Initialisation at the midplane +++++++++++++++++++++++++++++++
  z0 = 0.d0
  y0(1) = 1./gamma 
  y0(2) = 1./gamma
  y0(3) = 0.d0
  !WRITE(*,*) 'y:', y
  
  ! ++++++++++++++++++++++++++++++++ Integration to far far away ++++++++++++++++++++++++++++++
  dz0 = zfar/naccur
  DO i = 1, naccur
     CALL RKUTTA0(z0,y0,dz0,ncount)
     z0=z0+dz0
     IF (ncount .GT. ncountmax) THEN
        ncountmax = ncount
     ENDIF
  END DO
  
 
  IF (VERBOSE)  WRITE(*,*) 'y0 in z =', z0,' :', y0
  IF (VERBOSE)  WRITE(*,*) 'Runge Kutta a-t-il convergé? ncountmax :', ncountmax!, 'y :', y


  RETURN
END SUBROUTINE INTEGRE_Z
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE c2DICHOTOMY
  USE VARIABLES  
  use hydro_parameters
  use stratified
  IMPLICIT NONE
  ! Look for the right value of c2_cool with a dichotomy :
  ! local variables
  DOUBLE PRECISION :: c2_coolmin, c2_coolmax
  !DOUBLE PRECISION :: maxdc2_cool, prec, conv
  DOUBLE PRECISION, DIMENSION(1:3) :: y0
  INTEGER :: niteration

 
  c2_coolmax = 1.d0
  c2_coolmin = -1.d0
  niteration = 0
  ! Divide c2_cool by 2 until it is below what we want :
  DO WHILE((c2_coolmin .LT. 0.d0).AND.(niteration .LT. 40))
     niteration = niteration + 1
     c2_cool = c2_coolmax/2.d0
     CALL INTEGRE0(y0)
     WRITE(*,*) 'C2_COOL :', c2_cool, '  c2 :', gamma*y0(2)
     WRITE(*,*) 'y0 : ', y0
     IF (gamma*y0(2).LT.c2_cool) THEN
        c2_coolmax = c2_cool
     ELSE
        c2_coolmin = c2_cool
     END IF
  END DO 
  IF (niteration .GE. 40) WRITE(*,*) 'Dichotomy has not converged !!'

  niteration = 0
  ! Find c2_cool by dichotomy :
  DO WHILE(((c2_coolmax-c2_coolmin)/c2_coolmin .GT. 1.d-12).AND.(niteration .LT. 60))
     niteration = niteration + 1
     c2_cool = (c2_coolmax+c2_coolmin)/2.d0
     CALL INTEGRE0(y0)
     WRITE(*,*) 'C2_COOL :', c2_cool, '  c2 :', gamma*y0(2)
     WRITE(*,*) 'y0 : ', y0
     IF (gamma*y0(2).LT.c2_cool) THEN
        c2_coolmax = c2_cool
     ELSE
        c2_coolmin = c2_cool
     END IF
  END DO
  IF (niteration .GE. 40) WRITE(*,*) 'Dichotomy has not converged !!'
  
  ! Final integration for the fun :
  c2_cool = (c2_coolmax+c2_coolmin)/2.d0
  CALL INTEGRE0(y0)
  WRITE(*,*) 'C2_COOL :', c2_cool
  WRITE(*,*) 'y0 at the end of the dichotomy : ', y0


  RETURN
END SUBROUTINE c2DICHOTOMY
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE rkutta0(x0,y0,dx0,ncount)
  USE hydro_parameters
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(1:3) :: y0
  REAL(DP)   ::x0,dx0
  INTEGER*4	ncount,nvar0,nst,maxcount
  REAL(DP)   ::a(2,2),b(2),c(2)
  !--------------------variables locales--------------	
  INTEGER*4 	i,n,icount,j
  DOUBLE PRECISION, DIMENSION(1:3) :: yst,dxy
  DOUBLE PRECISION, DIMENSION(1:2,1:3) :: k,knew,dk	
  DOUBLE PRECISION, DIMENSION(1:2) :: xst
  DOUBLE PRECISION :: sup
  !-----------------------------------------------
  
  CALL initrk(a,b,c)
  nst = 2
  nvar0 = 3
  maxcount = 100
 

  DO i = 1,nst
     xst(i) = x0+c(i)*dx0
     DO n = 1,nvar0
        k(i,n) = 0.d0
     ENDDO
  ENDDO
  
  icount = 0
3 sup = 0.d0
  DO i = 1,nst
     DO n = 1,nvar0
        yst(n) = y0(n)
     ENDDO
     DO n = 1,nvar0
        DO j = 1,nst
           yst(n) = yst(n)+a(i,j)*k(j,n)
        ENDDO
     ENDDO
     CALL derive0(xst(i),yst,dxy)
     
     DO n = 1,nvar0
        knew(i,n) = dx0*dxy(n)
        dk(i,n)   = ABS(knew(i,n)-k(i,n))
        k(i,n)    = knew(i,n)
        IF (dk(i,n).GT.sup) sup = dk(i,n)
     ENDDO
  ENDDO
  
  icount = icount+1
  IF ((sup.GE.1.d-10).AND.(icount.LT.maxcount)) GOTO 3
  
  ncount = icount
  DO n = 1,nvar0
     DO i = 1,nst
        y0(n) = y0(n)+b(i)*k(i,n)
     ENDDO
  ENDDO
  

  RETURN
END SUBROUTINE rkutta0
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE initRK(a,b,c)
  USE hydro_parameters
  IMPLICIT NONE
  REAL(DP) :: a(2,2),b(2),c(2)
  !------variables locales------------------	
  REAL(DP) :: p,q,r,t			
  !-------------------------------------	
  
  
  ! coeffs. de la methode implicite d'ordre 4,nst = 2 
  a(1,1) = 0.25d0
  a(1,2) = 0.25d0-dsqrt(3.d0)/6.d0
  a(2,1) = 0.25d0+dsqrt(3.d0)/6.d0
  a(2,2) = 0.25d0
  
  b(1) = 0.5d0
  b(2) = 0.5d0
  
  c(1) = 0.5d0-dsqrt(3.d0)/6.d0
  c(2) = 0.5d0+dsqrt(3.d0)/6.d0
  
  
  RETURN
END SUBROUTINE initRK
!====================================================================
!====================================================================
!====================================================================
