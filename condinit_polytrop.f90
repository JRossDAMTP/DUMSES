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
  real(dp) :: H,beta,B0,p0,amp,rhofloor,polytrop
  real(dp) :: rho,rhovx,rhovy,rhovz,Bxc,Byc,Bzc,E,csq,P
  real :: rvalue
  !
  ! type :
  ! 'r': random, 'h': henrik's r-modes
  character(len=10) :: type
  !character(len=10) :: nchar,Rm,Re !mode number (character)
  !character(len=256) :: filename,dirname
  integer :: nz_eigen
  real(dp) :: interpolated, P_mid, H_henrik, k_eigen, omega_eigen, pi
  real(dp), dimension(:), allocatable :: z_eigen,ux_eigen,uy_eigen,uz_eigen,rho_eigen,P_eigen

  !
  namelist /init_params/d0,beta,gamma,polytrop,amp,type,floor,zfloor,dfloor,smooth,addmass,massDiff,stratif,rho_cool, c2_cool, t_cool
  !
  ! Calculate the initial state for uin...
  ! uin   (1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)
  ! uin = (\rho, \rho vx, \rho vy, \rho vz, Etot, bx, by, bz)

  pi=2.d0*asin(1.d0)

  gamma = 5./3.
  d0=1.d0
  beta=100.d0
  amp=1.d-1
  type='r'
  floor=.false. ; smooth=.false. ; addmass=.false. ; massDiff=.false. ; stratif = .true.
  zfloor=5.
  ! parameters of the cooling function :
  rho_cool=1.d-2  ! density cutoff
  c2_cool=0.1     ! equilibrium temperature
  t_cool = -1.d0  ! cooling time, if t_cool negative then there is no cooling

  open(unit=1,file='input' ,status='old')
  read(1,init_params)
  close(1)

  IF (type.EQ.'h') THEN
     ! read parameters of the linear modes :
     open(unit=1,file='AdiabWaveParam.dat',status='old')
     read(1,*) polytrop,gamma,k_eigen,omega_eigen,nz_eigen
     close(1)
     WRITE(*,*) 'polytrop :', polytrop, 'gamma :', gamma, 'k_eigen :', k_eigen, 'omega_eigen :', omega_eigen, 'nz_eigen :', nz_eigen
     ! read the grid :
     allocate(z_eigen(nz_eigen))
     open(unit=1,file='AdiabWaveGrid.dat',status='old')
     read(1,*) z_eigen
     close(1)
     IF (verbose) WRITE(*,*) 'z_eigen :', z_eigen
     ! read the eigenfunctions :
     allocate(ux_eigen(nz_eigen),uy_eigen(nz_eigen),uz_eigen(nz_eigen),rho_eigen(nz_eigen),P_eigen(nz_eigen))
     open(unit=1,file='AdiabWave.dat',status='old')
     read(1,*) ux_eigen,uy_eigen,uz_eigen,rho_eigen,P_eigen
     close(1)
     !IF (verbose) WRITE(*,*) 'ux_eigen :', ux_eigen
     !IF (verbose) WRITE(*,*) 'P_eigen :', P_eigen
  END IF
!###################
  polytrop=1
!###################
  H=ciso/Omega0
  H_henrik = H*sqrt(2.d0*(polytrop + 1.d0)/gamma)
  P_mid = 1.d0/gamma
  B0=0.d0

  WRITE(*,*) 'H_henrik :', H_henrik, '   P_mid :', P_mid

 
  uin=0.d0
  iseed=mype


  !Set vertical profile
  do l=1,ngrid
     do k=ku1,ku2
        ! background density (polytropic disc)
        uin(l,:,:,k,1)=d0*(1.D0 - gamma/(2.D0*(polytrop+1.d0))*z(k)**2/H**2 )**(polytrop)
        ! perturbations :
        IF (type .EQ. 'h') THEN
           do j=1,ny
              do i=1,nx
                 CALL GetInterpolated(ABS(z(k)/H_henrik),interpolated,z_eigen,rho_eigen,nz_eigen)
                 uin(l,i,j,k,1) = uin(l,i,j,k,1) + amp*sign(1.d0,z(k))*interpolated*cos(k_eigen*x(i)/H_henrik)               
              END do
           END do
        END IF
        
     end do
  end do

  do l=1,ngrid
     do k=1,nz+1
        IF (type .EQ. 'r') THEN
           !Momentum (Random condition)
           do j=1,ny
              do i=1,nx
                 call ran2(iseed,rvalue)
                 uin(l,i,j,k,2) = uin(l,i,j,k,1)*amp*(rvalue-0.5)*sqrt(ciso**2)
                 call ran2(iseed,rvalue)
                 uin(l,i,j,k,4) = uin(l,i,j,k,1)*amp*(rvalue-0.5)*sqrt(ciso**2)
                 call ran2(iseed,rvalue)
                 uin(l,i,j,k,3) = uin(l,i,j,k,1)*amp*(rvalue-0.5)*sqrt(ciso**2)
              end do
           end do
        ELSEIF(type .EQ. 'h') THEN
           do j=1,ny
              do i=1,nx
                 CALL GetInterpolated(ABS(z(k)/H_henrik),interpolated,z_eigen,ux_eigen,nz_eigen)
                 uin(l,i,j,k,2) = uin(l,i,j,k,1)*amp*sign(1.d0,z(k))*interpolated*cos(k_eigen*x(i)/H_henrik)*H_henrik*Omega0
                 CALL GetInterpolated(ABS(z(k)/H_henrik),interpolated,z_eigen,uy_eigen,nz_eigen)
                 uin(l,i,j,k,3) = uin(l,i,j,k,1)*amp*sign(1.d0,z(k))*interpolated*sin(k_eigen*x(i)/H_henrik)*H_henrik*Omega0
                 CALL GetInterpolated(ABS(z(k)/H_henrik),interpolated,z_eigen,uz_eigen,nz_eigen)
                 uin(l,i,j,k,4) = uin(l,i,j,k,1)*amp*interpolated*sin(k_eigen*x(i)/H_henrik)*H_henrik*Omega0
              END do
           END do

        END IF
        


     enddo
  enddo

  ! magnetic field :
  uin(1,iu1:iu2,ju1:ju2,ku1:ku2,6:8) = 0.d0
  uin(1,iu1:iu2,ju1:ju2,ku1:ku2,nvar+1:nvar+3) = 0.d0
 

  !Set total energy density
#if ISO==0
  do k=ku1,ku2-1
     do j=ju1,ju2-1
        do i=iu1,iu2-1
           
           rho=uin(1,i,j,k,1)
           rhovx=uin(1,i,j,k,2) ; Bxc=half*(uin(1,i,j,k,6)+uin(1,i+1,j  ,k  ,6))
           rhovy=uin(1,i,j,k,3) ; Byc=half*(uin(1,i,j,k,7)+uin(1,i  ,j+1,k  ,7))
           rhovz=uin(1,i,j,k,4) ; Bzc=half*(uin(1,i,j,k,8)+uin(1,i  ,j  ,k+1,8))
           ! sound speed of background state :
           csq = ciso**2*(1.d0 - gamma/(2.d0*(polytrop+1.d0))*z(k)**2/H**2)
           ! Add pressure perturbations
           IF (type .EQ. 'h') THEN
              CALL GetInterpolated(ABS(z(k)/H_henrik),interpolated,z_eigen,P_eigen,nz_eigen)
              P = d0*(1.D0 - gamma/(2.D0*(polytrop+1.d0))*z(k)**2/H**2 )**(polytrop)*csq/gamma + amp*sign(1.d0,z(k))*interpolated*P_mid*cos(k_eigen*x(i)/H_henrik)
              csq = gamma*P/rho
           END IF
           ! initialise energy :
           call get_E(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,csq)
           uin(1,i,j,k,5)=E
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

