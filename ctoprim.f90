!> \file
!! Contains subroutine ctoprim()
!<

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Subroutine CTOPRIM
!
!> Converts conservative variables to primitive in mhd.
!<
subroutine ctoprim(uin,q,bf,gravin,dt,ngrid)
  use hydro_parameters
  use const  
  use variables , only : cartesian
  implicit none

  integer :: ngrid

  real(dp), dimension(1:nvector,iu1:iu2  ,ju1:ju2  ,ku1:ku2  ,1:nvar+3) :: uin
  real(dp), dimension(1:nvector,iu1:iu2  ,ju1:ju2  ,ku1:ku2  ,1:ndim  ) :: gravin
  real(dp), dimension(1:nvector,iu1:iu2  ,ju1:ju2  ,ku1:ku2  ,1:nvar  ) :: q  
  real(dp), dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3     ) :: bf
  real(dp), dimension(1:nvector) :: eken,emag

  integer :: i,j,k,l,n,idim

  real(dp) :: eint,smalle,smallp,dt
  real(dp) :: dvx,dvy,lambda,ratio,rc
  
#if ISO==1
  smalle = 1.
  smallp = smallr*ciso**2
#else
  smalle = smallc**2/gamma/(gamma-one)
  smallp = smallr*smallc**2/gamma
#endif

  ! Store face centered magnetic field
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2+1
           do l = 1, ngrid
              if(i<=iu2)then
                 bf(l,i,j,k,1) = uin(l,i,j,k,6)
              else
                 bf(l,i,j,k,1) = uin(l,i-1,j,k,nvar+1)
              endif
           end do
        end do
     end do
  end do
  do k = ku1, ku2
     do j = ju1, ju2+1
        do i = iu1, iu2
           do l = 1, ngrid
              if(j<=ju2)then
                 bf(l,i,j,k,2) = uin(l,i,j,k,7)
              else
                 bf(l,i,j,k,2) = uin(l,i,j-1,k,nvar+2)
              endif
           end do
        end do
     end do
  end do
  do k = ku1, ku2+1
     do j = ju1, ju2
        do i = iu1, iu2
           do l = 1, ngrid
              if(k<=ku2)then
                 bf(l,i,j,k,3) = uin(l,i,j,k,8)
              else
                 bf(l,i,j,k,3) = uin(l,i,j,k-1,nvar+3)
              endif
           end do
        end do
     end do
  end do
  
  ! convert to primitive variable
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2

           ! compute density
           do l = 1, ngrid
              q(l,i,j,k,1) = max(uin(l,i,j,k,1),smallr)
           end do

           ! compute velocities
           do l = 1, ngrid
              q(l,i,j,k,2) = uin(l,i,j,k,2)/uin(l,i,j,k,1)
              q(l,i,j,k,3) = uin(l,i,j,k,3)/uin(l,i,j,k,1)
              q(l,i,j,k,4) = uin(l,i,j,k,4)/uin(l,i,j,k,1)
           end do

           ! compute cell centered magnetic field
           do l = 1, ngrid
              q(l,i,j,k,6) = half*(uin(l,i,j,k,6)+uin(l,i,j,k,nvar+1))
              q(l,i,j,k,7) = half*(uin(l,i,j,k,7)+uin(l,i,j,k,nvar+2))
              q(l,i,j,k,8) = half*(uin(l,i,j,k,8)+uin(l,i,j,k,nvar+3))
           end do

           ! Compute specific kinetic energy and magnetic energy
           do l = 1, ngrid
              eken(l) = half*(q(l,i,j,k,2)**2+q(l,i,j,k,3)**2+q(l,i,j,k,4)**2)
              emag(l) = half*(q(l,i,j,k,6)**2+q(l,i,j,k,7)**2+q(l,i,j,k,8)**2)
           end do

           ! Compute pressure
           do l = 1, ngrid
#if ISO==1
              q(l,i,j,k,5)=q(l,i,j,k,1)*ciso**2
#else
              eint = ( uin(l,i,j,k,5) - emag(l) ) /uin(l,i,j,k,1)-eken(l)
              q(l,i,j,k,5)=max((gamma-one)*q(l,i,j,k,1)*eint,smallp)
#endif
           end do

           ! Coriolis force predictor step
           if ((cartesian).and.(Omega0>0.d0)) then
              do l = 1, ngrid
                 dvx= 2.d0 *Omega0* q(l,i,j,k,3)
                 dvy=-5.d-1*Omega0* q(l,i,j,k,2)
                 q(l,i,j,k,2) = q(l,i,j,k,2) + dvx*dt*half
                 q(l,i,j,k,3) = q(l,i,j,k,3) + dvy*dt*half
              end do
           endif



!           ! gravity predictor step
!           do idim = 1, ndim
!              do l = 1, ngrid
!                 q(l,i,j,k,idim+1) = q(l,i,j,k,idim+1) + gravin(l,i,j,k,idim)*dt*half
!              end do
!           end do

        end do
     end do
  end do

!  ! Passive scalar
!  do n = 9, nvar
!     do k = ku1, ku2
!        do j = ju1, ju2
!           do i = iu1, iu2
!              do l = 1, ngrid
!                 q(l,i,j,k,n) = uin(l,i,j,k,n)/uin(l,i,j,k,1)
!              end do
!           end do
!        end do
!     end do
!  end do

  return

end subroutine ctoprim
