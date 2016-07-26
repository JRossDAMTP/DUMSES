! uin = (\rho, \rho vx, \rho vy, \rho vz, Etot, bx, by, bz)
subroutine thermal_diffusion(uin,flux,dx,dy,dz,dt,nu)
  use stratified
  use hydro_parameters
  use const
  implicit none
  real(dp)::dx,dy,dz,dt,nu
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin 
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
!-----------------------------------------------------------------------
  integer :: i,j,k,i0,j0,k0,idim
  real(dp) :: rho,P,rhoR,rhoL,PR,PL,Bxc,Byc,Bzc,kappa_loc
  real(dp) :: gradT
  
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
!=========================================================
!
  if (verbose) write (*,*) 'Entering thermal_diffusion...'

  do k=kf1,kf2
     do j=jf1,jf2
        do i=if1,if2

           !1st direction energy flux
           rhoR = uin(1,i,j,k,1)
           rhoL = uin(1,i-1,j,k,1)
           rho=half*(rhoL + rhoR)
           Bxc = half*(uin(1,i,j,k,6)+uin(1,i,j,k,nvar+1))
           Byc = half*(uin(1,i,j,k,7)+uin(1,i,j,k,nvar+2))
           Bzc = half*(uin(1,i,j,k,8)+uin(1,i,j,k,nvar+3))
           CALL get_P(uin(1,i,j,k,1),uin(1,i,j,k,2),uin(1,i,j,k,3),uin(1,i,j,k,4),uin(1,i,j,k,5),Bxc,Byc,Bzc,PR)
           Bxc = half*(uin(1,i-1,j,k,6)+uin(1,i-1,j,k,nvar+1))
           Byc = half*(uin(1,i-1,j,k,7)+uin(1,i-1,j,k,nvar+2))
           Bzc = half*(uin(1,i-1,j,k,8)+uin(1,i-1,j,k,nvar+3))
           CALL get_P(uin(1,i-1,j,k,1),uin(1,i-1,j,k,2),uin(1,i-1,j,k,3),uin(1,i-1,j,k,4),uin(1,i-1,j,k,5),Bxc,Byc,Bzc,PL)
           P = half*(PL + PR)
           gradT = (PR/rhoR - PL/rhoL)/dx
           CALL GET_KAPPA(kappa_loc,rho,P)

           flux(1,i,j,k,5,1)=-gamma/(gamma-1.D0)*rho*kappa_loc*gradT
           flux(1,i,j,k,5,1)=flux(1,i,j,k,5,1)*dt/dx


           !2nd direction energy flux
#if NDIM > 1
           rhoL = uin(1,i,j-1,k,1)
           rho=half*(rhoL + rhoR)
           Bxc = half*(uin(1,i,j-1,k,6)+uin(1,i,j-1,k,nvar+1))
           Byc = half*(uin(1,i,j-1,k,7)+uin(1,i,j-1,k,nvar+2))
           Bzc = half*(uin(1,i,j-1,k,8)+uin(1,i,j-1,k,nvar+3))
           CALL get_P(uin(1,i,j-1,k,1),uin(1,i,j-1,k,2),uin(1,i,j-1,k,3),uin(1,i,j-1,k,4),uin(1,i,j-1,k,5),Bxc,Byc,Bzc,PL)
           P = half*(PL + PR)
           gradT = (PR/rhoR - PL/rhoL)/dy
           CALL GET_KAPPA(kappa_loc,rho,P)

           flux(1,i,j,k,5,2)=-gamma/(gamma-1.D0)*rho*kappa_loc*gradT
           flux(1,i,j,k,5,2)=flux(1,i,j,k,5,2)*dt/dy

#endif

           !3rd direction energy flux
#if NDIM==3
           rhoL = uin(1,i,j,k-1,1)
           rho=half*(rhoL + rhoR)
           Bxc = half*(uin(1,i,j,k-1,6)+uin(1,i,j,k-1,nvar+1))
           Byc = half*(uin(1,i,j,k-1,7)+uin(1,i,j,k-1,nvar+2))
           Bzc = half*(uin(1,i,j,k-1,8)+uin(1,i,j,k-1,nvar+3))
           CALL get_P(uin(1,i,j,k-1,1),uin(1,i,j,k-1,2),uin(1,i,j,k-1,3),uin(1,i,j,k-1,4),uin(1,i,j,k-1,5),Bxc,Byc,Bzc,PL)
           P = half*(PL + PR)
           gradT = (PR/rhoR - PL/rhoL)/dz
           CALL GET_KAPPA(kappa_loc,rho,P)

           flux(1,i,j,k,5,3)=-gamma/(gamma-1.D0)*rho*kappa_loc*gradT
           flux(1,i,j,k,5,3)=flux(1,i,j,k,5,3)*dt/dz

#endif

        end do
     end do
  end do
  !
  ! Conservative update for energy
  !
  do idim=1,ndim
     i0=0; j0=0; k0=0 
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1

     do k=1,nz
        do j=1,ny
           do i=1,nx
              uin(1,i,j,k,5)=uin(1,i,j,k,5)+ &
                   & (flux(1,i   ,j   ,k   ,5,idim) &
                   & -flux(1,i+i0,j+j0,k+k0,5,idim))
           end do
        end do
     end do
     
  end do


  
  
!
  return
end subroutine thermal_diffusion
!
