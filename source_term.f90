! uin = (\rho, \rho vx, \rho vy, \rho vz, Etot, bx, by, bz)
subroutine source_term
  use hydro_parameters
  use const
  use variables , only : dt,uin,old_uin,gravin,ngrid,z,dx,dy,dz
  use stratified
  implicit none

  integer :: i,j,k,l

  real(dp) :: rhon,rhonp1,rho,vxn,vyn,vxnp1,vynp1,vx,vy ! Gravity force parameters
  real(dp) :: Bxn,Byn,Bxnp1,Bynp1,Bzn,Bznp1,Bx,By,Bz,rhovx,rhovy,rhovz,E,csq,cool,P
  real(dp) :: Bxc,Byc,Bzc,rho_old,Eint
  real(dp) :: theta,smallH,pi,mass !Addmass parameters

  if (verbose) write (*,*) 'Entering source_term subroutine...'

  
  do k=1,nz
     do j=1,ny
        do i=1,nx
           do l=1,ngrid
              !Gravity term
              !density
              rhon  =old_uin(l,i,j,k,1)
              rhonp1=    uin(l,i,j,k,1)

              !update momentum
                          uin(l,i,j,k,2)=uin(l,i,j,k,2)+half*(rhon+rhonp1)*gravin(l,i,j,k,1)*dt
              if (ndim>1) uin(l,i,j,k,3)=uin(l,i,j,k,3)+half*(rhon+rhonp1)*gravin(l,i,j,k,2)*dt
              if (ndim>2) uin(l,i,j,k,4)=uin(l,i,j,k,4)+half*(rhon+rhonp1)*gravin(l,i,j,k,3)*dt


#if ISO==0
              !update energy
              ! Gravity source term
              uin(1,i,j,k,5) = uin(1,i,j,k,5) + &
                      half*(old_uin(l,i,j,k,2)+uin(1,i,j,k,2))*gravin(l,i,j,k,1)*dt + &
                      half*(old_uin(l,i,j,k,3)+uin(1,i,j,k,3))*gravin(l,i,j,k,2)*dt + &
                      half*(old_uin(l,i,j,k,4)+uin(1,i,j,k,4))*gravin(l,i,j,k,3)*dt

              ! Energy source term in shearing box
              ! compute vx,vy and rho centered in time
              vxn  = old_uin(1,i,j,k,2)/rhon
              vyn  = old_uin(1,i,j,k,3)/rhon
              vxnp1=     uin(1,i,j,k,2)/rhonp1
              vynp1=     uin(1,i,j,k,3)/rhonp1
              vx  = half*(vxn+vxnp1)
              vy  = half*(vyn+vynp1)
              rho = half*(rhon+rhonp1)
              ! compute Bx and By centered in the cell and in time
              Bxn   = half*(old_uin(1,i,j,k,6) + old_uin(1,i,j,k,nvar+1))
              Bxnp1 = half*(    uin(1,i,j,k,6) +     uin(1,i,j,k,nvar+1))
              Bx  = half*(Bxn+Bxnp1)
              Byn   = half*(old_uin(1,i,j,k,7) + old_uin(1,i,j,k,nvar+2))
              Bynp1 = half*(    uin(1,i,j,k,7) +     uin(1,i,j,k,nvar+2))
              By  = half*(Byn+Bynp1)
              ! update energy
              uin(1,i,j,k,5) = uin(1,i,j,k,5) + 1.5d0*Omega0*dt*(rho*vx*vy - Bx*By)

              ! Cooling source term :
              IF (t_cool .GT. 0.d0) THEN
		
                 ! compute csq centered in time :
                 Bzn   = half*(old_uin(1,i,j,k,8) + old_uin(1,i,j,k,nvar+3))
                 Bznp1 = half*(    uin(1,i,j,k,8) +     uin(1,i,j,k,nvar+3))
                 Bz    = half*(Bzn+Bznp1)
                 E = half*(uin(1,i,j,k,5) + old_uin(1,i,j,k,5))
                 rhovx = half*(uin(1,i,j,k,2) + old_uin(1,i,j,k,2))
                 rhovy = half*(uin(1,i,j,k,3) + old_uin(1,i,j,k,3))
                 rhovz = half*(uin(1,i,j,k,4) + old_uin(1,i,j,k,4))
                 CALL GET_csq(rho,rhovx,rhovy,rhovz,E,Bx,By,Bz,gamma,csq)
                 ! compute cooling function
                  CALL get_cool(cool,rho,csq,rhovx,rhovy,rhovz,E,Bx,By,Bz,gamma)
                 ! update energy
                 uin(1,i,j,k,5) = uin(1,i,j,k,5) + cool*dt
              ENDIF
              
#endif

           end do
        end do
     end do
  end do

  !Set density field to floor
  if (floor) then
     do k=1,nz
        do j=1,ny
           do i=1,nx
              do l=1,ngrid
                 
                 !density
#if ISO==1
                 uin(l,i,j,k,1)=max(uin(l,i,j,k,1),exp(-zfloor**2/2.))
#else
                 uin(l,i,j,k,1)=max(uin(l,i,j,k,1),dfloor)
#endif
                 
              end do
           end do
        end do
     end do
  endif  

  !Add source term in continuity equation
  if (addmass) then

     smallH=0.2 ; pi=2.d0*asin(1.d0)

     mass=0.d0
     do k=1,nz
        do j=1,ny
           do i=1,nx
              mass = mass + uin(1,i,j,k,1)
           end do
        end do
     end do
#if WITHMPI==1
     call sumBcast(mass,1)
#endif
     mass=mass*dx*dy*dz
     theta=(mass0-mass)/sqrt(2.*pi)/smallH/(xmax-xmin)/(ymax-ymin)
     
     !restore initial mass
     do k=1,nz
        do j=1,ny
           do i=1,nx

#if ISO ==0

	      E=uin(l,i,j,k,5)
	      Bxc=half*(    uin(1,i,j,k,6) +     uin(1,i,j,k,nvar+1))
              Byc=half*(    uin(1,i,j,k,7) +     uin(1,i,j,k,nvar+2))
              Bzc=half*(    uin(1,i,j,k,8) +     uin(1,i,j,k,nvar+3))
	      Eint=E-half*(Bxc*Bxc+Byc*Byc+Bzc*Bzc)
	      Eint=Eint*(uin(1,i,j,k,1)+theta*exp(-z(k)**2/2.d0/smallH**2))/uin(1,i,j,k,1)
	      uin(l,i,j,k,5)=Eint+half*(Bxc*Bxc+Byc*Byc+Bzc*Bzc)
     
#endif

              uin(1,i,j,k,2:4)=uin(1,i,j,k,2:4)/uin(1,i,j,k,1)
              uin(1,i,j,k,1)=uin(1,i,j,k,1)*(1+(mass0-mass)/mass)!theta*exp(-z(k)**2/2.d0/smallH**2)
              uin(1,i,j,k,2:4)=uin(1,i,j,k,2:4)*uin(1,i,j,k,1)




           end do
        end do
     end do
              
  endif

  if (verbose) write (*,*) 'End of source_term subroutine...'

end subroutine source_term
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine checkLowBeta(q,dq,dbf,flagCell,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer :: ngrid
  real(dp), dimension(1:nvector,iu1:iu2  ,ju1:ju2  ,ku1:ku2  ,1:nvar       ) :: q 
  real(dp), dimension(1:nvector,iu1:iu2  ,ju1:ju2  ,ku1:ku2  ,1:nvar,1:ndim) :: dq
  real(dp), dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2      ,1:3   ,1:ndim) :: dbf
  logical , dimension(1:nvector,iu1:iu2  ,ju1:ju2  ,ku1:ku2 ) :: flagCell

  real(dp) :: rho,bx,by,bz,pmag,pth,beta
  real(dp) :: minBeta
  integer  :: i,j,k,l,ncell
  integer  :: ilo,ihi,jlo,jhi,klo,khi

  minBeta=1.d-3 ; ncell=0 ; flagCell=.False.

  ilo=min(1,iu1+3); ihi=max(1,iu2-3)
  jlo=min(1,ju1+3); jhi=max(1,ju2-3)
  klo=min(1,ku1+3); khi=max(1,ku2-3)

  do k=klo,khi
     do j=jlo,jhi
        do i=ilo,ihi
           do l=1,ngrid
           
              rho= q(l,i,j,k,1)
              bx = q(l,i,j,k,6)
              by = q(l,i,j,k,7)
              bz = q(l,i,j,k,8)
           
              pth  = rho*ciso**2
              pmag = half*(bx*bx+by*by+bz*bz)

              beta=pth/pmag

              if (beta<minBeta) then
                 !dq(l,i,j,k,:,:)=zero
                 !dbf(l,i:i+1,j:j+1,k:k+1,:,:)=zero
                 !flagCell(l,i:i+1,j:j+1,k:k+1)=.True.
                 !ncell=ncell+1
              endif

           end do
        end do
     end do
  end do

  if (ncell>0) write (*,*) '   Nb of cells affected by low beta:',ncell

  return
end subroutine checkLowBeta
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine gravity_predictor(v,gravin,dt,igrav,jgrav,kgrav,ngrid)
  use hydro_parameters
  use const
  implicit none

  integer :: ngrid
  logical :: igrav,jgrav,kgrav
  real(dp), dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3   ) :: v
  real(dp), dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim) :: gravin
  real(dp) :: dt

  integer :: i,j,k,l
  integer :: ilo,ihi,jlo,jhi,klo,khi

  ilo = min(1,iu1+1) ; ihi = max(1,iu2-1)
  jlo = min(1,ju1+1) ; jhi = max(1,ju2-1)
  klo = min(1,ku1+1) ; khi = max(1,ku2-1)

  !v=v+gravin*dt*half
  do k=klo,khi
     do j=jlo,jhi
        do i=ilo,ihi
           do l=1,ngrid

              if (igrav)                v(l,i,j,k,1) = v(l,i,j,k,1) + half*dt*gravin(l,i,j,k,1)
              if ((jgrav).and.(ndim>1)) v(l,i,j,k,2) = v(l,i,j,k,2) + half*dt*gravin(l,i,j,k,2)
              if ((kgrav).and.(ndim>2)) v(l,i,j,k,3) = v(l,i,j,k,3) + half*dt*gravin(l,i,j,k,3)

           end do
        end do
     end do
  end do

  return
end subroutine gravity_predictor
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!subroutine get_phi_SB(phi_SB,x,z)
!  use hydro_parameters , only : Omega0
!  implicit none
!  DOUBLE PRECISION , intent(in ) :: x,z
!  DOUBLE PRECISION , intent(out) :: phi_SB
!
!
!  phi_SB = -1.5d0*Omega0**2*x**2 + 0.5d0*Omega0**2*z**2
!  
!  
!  return
!end subroutine get_phi_SB
