subroutine user_init
  use hydro_parameters
  use const
  use variables , only : gravin,z,dx,dy,dz,ngrid
  use stratified
  implicit none

  integer :: i,j,k,l,ipos,idim
  integer :: ilo,ihi,jlo,jhi,klo,khi
  real(dp) :: xc,yc,phic,dist
  real(dp), dimension(2,ndim) :: phi

  !Compute gravin array
  ilo = min(1,iu1+1) ; ihi = max(1,iu2-1)
  jlo = min(1,ju1+1) ; jhi = max(1,ju2-1)
  klo = min(1,ku1+1) ; khi = max(1,ku2-1)

  gravin=0.d0 ; phi=0.d0
  do k=klo,khi
     do j=jlo,jhi
        do i=ilo,ihi

           IF (stratif) THEN
              phi(1,3)=half*Omega0**2*z(k-1)**2
              phi(2,3)=half*Omega0**2*z(k+1)**2
              if (smooth) then
                 if (abs(z(k-1))>zfloor) phi(1,3)=half*Omega0**2*zfloor**2
                 if (abs(z(k+1))>zfloor) phi(2,3)=half*Omega0**2*zfloor**2
              endif
           ENDIF

           do l=1,ngrid
              !First coordinate gravity force
              gravin(l,i,j,k,1)=-half*(phi(2,1)-phi(1,1))/dx
              !Second coordinate gravity force
              if (ndim>1) gravin(l,i,j,k,2)=-half*(phi(2,2)-phi(1,2))/dy
              !Third coordinate gravity force
              if (ndim>2) gravin(l,i,j,k,3)=-half*(phi(2,3)-phi(1,3))/dz
           end do

        end do
     end do
  end do

  return
end subroutine user_init
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine get_eta(etaval,r)
  use amr_parameters
  use variables , only : eta
  use stratified
  implicit none
  real(dp) :: etaval,r

  if (resist==.true.) then
      if ((r>2.d0) .or. (r<-2.d0)) then
	 etaval=eta
      else 
         etaval=0.d0
      end if
  else 
  	etaval=eta
  endif
  return
end subroutine get_eta
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine get_kappa(kappaval,rho,P)
  use amr_parameters
  use variables , only : kappa
  implicit none
  real(dp) :: kappaval,rho,P

  ! Constant kappa :
  kappaval=kappa

  ! kappa depending on rho and P : to be implemented


  return
end subroutine get_kappa
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine get_P(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,P)
  use amr_parameters
  use hydro_parameters , only : smallc, smallr,gamma
  use const
  implicit none
  real(dp) , intent(in ) :: rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc
  real(dp) , intent(out) :: P

  real(dp) :: Ec,Em,smallp

  smallp = smallr*smallc**2/gamma

  Ec=half*(rhovx*rhovx+rhovy*rhovy+rhovz*rhovz)/rho
  Em=half*(Bxc*Bxc+Byc*Byc+Bzc*Bzc)

  P = (gamma-1.d0)*(E-Ec-Em)
  P = max(P,smallp)
  
!  write (*,*) P,rhovz,rho,E,Ec,Em,Bxc,Byc,Bzc
  
  return
end subroutine get_P
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine get_csq(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,csq)
  use amr_parameters
  use hydro_parameters , only : smallc
  use const
  implicit none
  real(dp) , intent(in ) :: rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma
  real(dp) , intent(out) :: csq

  real(dp) :: Ec,Em

  Ec=half*(rhovx*rhovx+rhovy*rhovy+rhovz*rhovz)/rho
  Em=half*(Bxc*Bxc+Byc*Byc+Bzc*Bzc)

  csq=gamma*(gamma-1.d0)*(E-Ec-Em)/rho
  !csq=(gamma-1.d0)*(E-Ec-Em)/rho
  csq=max(csq,smallc*smallc)
!  write (*,*) csq,rhovz,rho,E,Ec,Em,Bxc,Byc,Bzc
  
  return
end subroutine get_csq
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine get_E(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,csq)
  use amr_parameters
  use const
  implicit none
  real(dp) , intent(in ) :: rho,rhovx,rhovy,rhovz,Bxc,Byc,Bzc,gamma,csq
  real(dp) , intent(out) :: E

  real(dp) :: Ec,Em

  Ec=half*(rhovx*rhovx+rhovy*rhovy+rhovz*rhovz)/rho
  Em=half*(Bxc*Bxc+Byc*Byc+Bzc*Bzc)

  E=rho*csq/gamma/(gamma-1.)+Ec+Em

  return
end subroutine get_E
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine get_cool(cool,rho,csq,rhovx,rhovy,rhovz,E,Bx,By,Bz,gamma)
  use amr_parameters
  use const
  use stratified
  use variables
!  use variables , only : cool_const
  implicit none
  real(dp) , intent(in) :: rho,csq,rhovx,rhovy,rhovz,E,Bx,By,Bz,gamma
  real(dp) , intent(out) :: cool
  real(dp) ::kin,mag,Press !,cool_const

 ! cool_const=0.02!.002!.2
!  real(dp) :: 
  kin=0.5d0*(rhovx**2+rhovy**2+rhovz**2)/rho
  mag=(Bx**2+By**2+Bz**2)/2
  Press=(E-kin-mag)*(1-gamma)
  cool = - rho*(csq - c2_cool)/t_cool*exp(-rho/rho_cool)
   cool = - cool_const*((Press)**2.d0)

  return
end subroutine get_cool

subroutine get_nu(nu,rho,rhovx,rhovy,rhovz,E,Bx,By,Bz)
  use amr_parameters
  use const
  use stratified
  use hydro_parameters
!  use variables , only : alpha
  implicit none
  real(dp) , intent(in) :: rho,rhovx,rhovy,rhovz,E,Bx,By,Bz
  real(dp) , intent(out) :: nu
  real(dp) :: kin,mag,Press !,alpha

!  alpha=0.05
!  real(dp) :: 
  kin=0.5d0*(rhovx**2+rhovy**2+rhovz**2)/rho
  mag=(Bx**2+By**2+Bz**2)/2
  Press=(E-kin-mag)*(gamma-1.d0)
  nu=Press*0.004d0/(Omega0*rho)
 ! nu=max(nu,1.d0)
  nu=1.d-1
  return
end subroutine get_nu

