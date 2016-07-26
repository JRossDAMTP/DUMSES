!###########################################################
!###########################################################
!###########################################################
!###########################################################
! uin = (\rho, \rho vx, \rho vy, \rho vz, Etot, bx, by, bz)
! q   = (\rho,      vx,      vy,      vz, P   , bx, by, bz)
subroutine compute_dt(uin,dx,dy,dz,dt,courant,ngrid)
  use hydro_parameters
  use const
  use variables , only : x,y,cartesian,cylindrical,spherical,fargo,nu,eta,kappa
  implicit none
#if WITHMPI==1
#include "mpif.h"
#endif
 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3) ::uin
  real(dp) :: dx,dy,dz,dt,courant, dt_kappa, dt_nu, dt_eta, kappa_loc
  integer  :: ngrid

  real(dp),dimension(:,:,:,:,:),allocatable :: q,bf,gravin
  real(dp) :: vx,vy,vz,bx,by,bz,p,rho,velx,vely,velz
  integer  :: i,j,k,l
#if WITHMPI==1
  real(dp) :: dt_mype
  integer :: ierr
#endif

  if (verbose) write (*,*) 'Entering compute_dt...'

  allocate(q     (1:nvector,iu1:iu2  ,ju1:ju2  ,ku1:ku2  ,1:nvar)) ; q=0.d0
  allocate(gravin(1:nvector,iu1:iu2  ,ju1:ju2  ,ku1:ku2  ,1:ndim)) ; gravin=0.d0
  allocate(bf    (1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3   )) ; bf=0.d0
  
  call ctoprim(uin,q,bf,gravin,0.d0,ngrid)

  dt=courant*dx/smallc
  dt_kappa=dt
  dt_nu=dt
  dt_eta=dt

  do l=1,ngrid

     do k=1,nz
        do j=1,ny
           do i=1,nx

              rho= q(l,i,j,k,1)
              vx = q(l,i,j,k,2)
              vy = q(l,i,j,k,3)
              vz = q(l,i,j,k,4)
              p  = q(l,i,j,k,5)
              bx = q(l,i,j,k,6)
              by = q(l,i,j,k,7)
              bz = q(l,i,j,k,8)

              !compute fastest signal speed (x dir)
              call find_speed_info((/rho,p,vx,bx,vy,by,vz,bz/),velx)

#if NDIM > 1
              !compute fastest signal speed (y dir)
              call find_speed_info((/rho,p,vy,by,vx,bx,vz,bz/),vely)
#endif

#if NDIM == 3
              !compute fastest signal speed (z dir)
              call find_speed_info((/rho,p,vz,bz,vx,bx,vy,by/),velz)
#endif

              ! Compute the thermal diffusivity kappa :
               CALL GET_KAPPA(kappa_loc,rho,P)

#if NDIM == 1
              dt=min(dt,dx/velx)
              if (kappa .GT. 0.d0) dt_kappa=min(dt_kappa,0.5d0*dx**2/kappa_loc)
              if (nu .GT. 0.d0)    dt_nu=min(dt_nu,0.5d0*dx**2/nu)
              if (eta .GT. 0.d0)   dt_eta=min(dt_nu,0.5d0*dx**2/eta)
#endif

#if NDIM == 2
              if (cartesian  ) then
                 dt=min(dt,one/(velx/dx + vely/dy))
                 if (kappa .GT. 0.d0) dt_kappa=min(dt_kappa,0.5d0*dx**2/kappa_loc,0.5d0*dy**2/kappa_loc)
                 if (nu .GT. 0.d0)    dt_nu   =min(dt_nu,   0.5d0*dx**2/nu,       0.5d0*dy**2/nu)
                 if (eta .GT. 0.d0)   dt_eta  =min(dt_nu,   0.5d0*dx**2/eta,      0.5d0*dy**2/eta)
              else
                 dt=min(dt,one/(velx/dx + vely/dy/x(i)))
                 if (kappa .GT. 0.d0) dt_kappa=min(dt_kappa,0.5d0*dx**2/kappa_loc,0.5d0*(dy*x(i))**2/kappa_loc)
                 if (nu .GT. 0.d0)    dt_nu   =min(dt_nu,   0.5d0*dx**2/nu,       0.5d0*(dy*x(i))**2/nu)
                 if (eta .GT. 0.d0)   dt_eta  =min(dt_nu,   0.5d0*dx**2/eta,      0.5d0*(dy*x(i))**2/eta)
              endif
#endif

#if NDIM == 3

              if ((cartesian).and.(Omega0>0.d0).and..not.(fargo)) vely=vely+1.5*Omega0*(xmax-xmin)/2.

              if (cartesian  ) then 
                 dt=min(dt,one/(velx/dx + vely/dy      + velz/dz               ) )
                 if (kappa .GT. 0.d0) dt_kappa=min(dt_kappa,0.5d0*dx**2/kappa_loc,0.5d0*dy**2/kappa_loc,0.5d0*dz**2/kappa_loc)
                 if (nu .GT. 0.d0)    dt_nu   =min(dt_nu,   0.5d0*dx**2/nu,       0.5d0*dy**2/nu,       0.5d0*dz**2/nu)
                 if (eta .GT. 0.d0)   dt_eta  =min(dt_nu,   0.5d0*dx**2/eta,      0.5d0*dy**2/eta,      0.5d0*dz**2/eta)
              elseif (cylindrical) then
                 dt=min(dt,one/(velx/dx + vely/dy/x(i) + velz/dz               ) )
                 if (kappa .GT. 0.d0) dt_kappa=min(dt_kappa,0.5d0*dx**2/kappa_loc,0.5d0*(dy*x(i))**2/kappa_loc,0.5d0*dz**2/kappa_loc)
                 if (nu .GT. 0.d0)    dt_nu   =min(dt_nu,   0.5d0*dx**2/nu,       0.5d0*(dy*x(i))**2/nu,       0.5d0*dz**2/nu)
                 if (eta .GT. 0.d0)   dt_eta  =min(dt_nu,   0.5d0*dx**2/eta,      0.5d0*(dy*x(i))**2/eta,      0.5d0*dz**2/eta)
              elseif (spherical  ) then
                 dt=min(dt,one/(velx/dx + vely/dy/x(i) + velz/dz/x(i)/sin(y(j))) )
                 if (kappa .GT. 0.d0) dt_kappa=min(dt_kappa,0.5d0*dx**2/kappa_loc,0.5d0*(dy*x(i))**2/kappa_loc,0.5d0*(dz*x(i)*sin(y(j)))**2/kappa_loc)
                 if (nu .GT. 0.d0)    dt_nu   =min(dt_nu,   0.5d0*dx**2/nu,       0.5d0*(dy*x(i))**2/nu,       0.5d0*(dz*x(i)*sin(y(j)))**2/nu)
                 if (eta .GT. 0.d0)   dt_eta  =min(dt_nu,   0.5d0*dx**2/eta,      0.5d0*(dy*x(i))**2/eta,      0.5d0*(dz*x(i)*sin(y(j)))**2/eta)
              end if
#endif
              
           end do
        end do
     end do

   end do

   ! factor ndim for the diffusion terms :
   dt_kappa = dt_kappa/ndim
   dt_nu = dt_nu/ndim
   dt_eta = dt_eta/ndim

#if WITHMPI==1
   dt_mype=courant*min(dt,dt_kappa,dt_nu,dt_eta)
   call MPI_Allreduce(dt_mype,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
     &                MPI_COMM_WORLD,ierr)
#else
   dt=courant*min(dt,dt_kappa,dt_nu,dt_eta)
#endif      

  deallocate(q,gravin,bf)

  return
end subroutine compute_dt
