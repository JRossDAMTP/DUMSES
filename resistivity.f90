! uin = (\rho, \rho vx, \rho vy, \rho vz, Etot, bx, by, bz)
subroutine resistivity(uin,flux,emfx,emfy,emfz,dx,dy,dz,dt,ngrid)
  use const
  use hydro_parameters
  use variables , only : x,y,z,cartesian,cylindrical,spherical
  implicit none
!
  real(dp), intent(in ) :: dx,dy,dz,dt,ngrid
  real(dp), intent(out),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin 
  real(dp), intent(out),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp), intent(out), dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::emfx,emfy,emfz
!----------------------------------------------------------------------------
  integer :: i,j,k,l,ibx,iby,ibz,i0,j0,k0,idim
  real(dp) :: eta,jx,jy,jz,xc,xL
  real(dp) :: dbxdy,dbxdz,dbydx,dbydz,dbzdx,dbzdy
#if ISO==0
  real(dp) :: jxp1,jyp1,jzp1,bx,by,bz
#endif
!
!=======================================================
!
  if (verbose) write (*,*) 'Entering resistivity...'

  ibx=6 ; iby=7 ; ibz=8

  dbxdy=0.d0 ; dbzdy=0.d0
  dbxdz=0.d0 ; dbydz=0.d0

  emfx=0.d0 ; emfy=0.d0 ; emfz=0.d0

  !Calculate J=curl(B)
  do k=kf1,kf2
     do j=jf1,jf2
        do i=if1,if2

           if (cartesian) then

              dbydx=(uin(1,i,j,k,iby)-uin(1,i-1,j  ,k  ,iby))/dx
              dbzdx=(uin(1,i,j,k,ibz)-uin(1,i-1,j  ,k  ,ibz))/dx
#if NDIM>1
              dbxdy=(uin(1,i,j,k,ibx)-uin(1,i  ,j-1,k  ,ibx))/dy
              dbzdy=(uin(1,i,j,k,ibz)-uin(1,i  ,j-1,k  ,ibz))/dy
#endif
#if NDIM==3
              dbxdz=(uin(1,i,j,k,ibx)-uin(1,i  ,j  ,k-1,ibx))/dz
              dbydz=(uin(1,i,j,k,iby)-uin(1,i  ,j  ,k-1,iby))/dz
#endif

#if NDIM>1
              jx = dbzdy-dbydz
#endif
              jy = dbxdz-dbzdx ; jz = dbydx-dbxdy

           endif
!
           if (cylindrical) then

              xL=half*(x(i)+x(i-1))
              xc=x(i)

              jz= 2.d0*(x(i)*uin(1,i,j,k,iby)-x(i-1)*uin(1,i-1,j,k,iby))/(x(i)**2-x(i-1)**2) - &
                     & (     uin(1,i,j,k,ibx)-       uin(1,i,j-1,k,ibx))/dy/xL
#if NDIM>2 
              jx = - (uin(1,i,j,k,iby)-uin(1,i,j,k-1,iby))/dz &
                 & + (uin(1,i,j,k,ibz)-uin(1,i,j-1,k,ibz))/dy/xc
              jy = (uin(1,i,j,k,ibx)-uin(1,i,j,k-1,ibx))/dz &
               & - (uin(1,i,j,k,ibz)-uin(1,i-1,j,k,ibz))/dx
#endif
           endif
           !
#if NDIM>1
           call get_eta(eta,z(k))
           emfx(1,i,j,k)=-eta*jx*dt
#endif
           call get_eta(eta,half*(z(k)+z(k-1)))
           emfy(1,i,j,k)=-eta*jy*dt
           emfz(1,i,j,k)=-eta*jz*dt
           !
        end do
     end do
  end do
  !
#if ISO==0
  do k=kf1,kf2
     do j=jf1,jf2
        do i=if1,if2
           
           !1st direction energy flux
#if NDIM==1
           call get_eta(eta,z(k))
           by=half*(uin(1,i,j,k,7)+uin(1,i-1,j,k,7))
#endif
#if NDIM>1
           by=forth*( uin(1,i,j  ,k,7)+uin(1,i-1,j  ,k,7) + &
                    & uin(1,i,j+1,k,7)+uin(1,i-1,j+1,k,7) )
           call get_eta(eta,half*(z(k)+z(k-1)))
#endif
#if NDIM<3
           bz=half*(uin(1,i,j,k,8)+uin(1,i-1,j,k,8))
#endif
#if NDIM==3
           bz=forth*( uin(1,i,j,k  ,8)+uin(1,i-1,j,k  ,8) + &
                    & uin(1,i,j,k+1,8)+uin(1,i-1,j,k+1,8))
#endif
           
#if NDIM<3        
           jy   = - (uin(1,i,j,k  ,8)-uin(1,i-1,j  ,k  ,8))/dx
#endif
#if NDIM==3           
           jy   = (uin(1,i,j,k  ,6)-uin(1,i  ,j  ,k-1,6))/dz &
            &   - (uin(1,i,j,k  ,8)-uin(1,i-1,j  ,k  ,8))/dx
           jyp1 = (uin(1,i,j,k+1,6)-uin(1,i  ,j  ,k  ,6))/dz &
            &   - (uin(1,i,j,k+1,8)-uin(1,i-1,j  ,k+1,8))/dx
           jy   = half*(jy+jyp1)
#endif
           
#if NDIM==1          
           jz   =   (uin(1,i,j  ,k,7)-uin(1,i-1,j  ,k  ,7))/dx
#endif
#if NDIM>1           
           jz   = (uin(1,i,j  ,k,7)-uin(1,i-1,j  ,k  ,7))/dx &
            &   - (uin(1,i,j  ,k,6)-uin(1,i  ,j-1,k  ,6))/dy
           jzp1 = (uin(1,i,j+1,k,7)-uin(1,i-1,j+1,k  ,7))/dx &
            &   - (uin(1,i,j+1,k,6)-uin(1,i  ,j  ,k  ,6))/dy
           jz   = half*(jz+jzp1)
#endif
           
           flux(1,i,j,k,5,1)=flux(1,i,j,k,5,1) - eta*(jy*bz-jz*by)*dt/dx
           
#if NDIM>1
           !2nd direction energy flux
           bx=forth*( uin(1,i  ,j,k,6)+uin(1,i  ,j-1,k,6) + &
                    & uin(1,i+1,j,k,6)+uin(1,i+1,j-1,k,6) )
#if NDIM==2
           bz=half*(uin(1,i,j,k,8)+uin(1,i,j-1,k,8))
#endif
#if NDIM==3
           bz=forth*( uin(1,i,j,k  ,8)+uin(1,i,j-1,k  ,8) + &
                    & uin(1,i,j,k+1,8)+uin(1,i,j-1,k+1,8))
#endif

#if NDIM==2
           jx   = (uin(1,i,j,k  ,8)-uin(1,i  ,j-1,k  ,8))/dy
#endif
#if NDIM==3
           jx   = (uin(1,i,j,k  ,8)-uin(1,i  ,j-1,k  ,8))/dy &
            &   - (uin(1,i,j,k  ,7)-uin(1,i  ,j  ,k-1,7))/dz
           jxp1 = (uin(1,i,j,k+1,8)-uin(1,i  ,j-1,k+1,8))/dy &
            &   - (uin(1,i,j,k+1,7)-uin(1,i  ,j  ,k  ,7))/dz
           jx=half*(jx+jxp1)
#endif
           
           jz   = (uin(1,i  ,j,k,7)-uin(1,i-1,j  ,k  ,7))/dx &
            &   - (uin(1,i  ,j,k,6)-uin(1,i  ,j-1,k  ,6))/dy
           jzp1 = (uin(1,i+1,j,k,7)-uin(1,i  ,j  ,k  ,7))/dx &
            &   - (uin(1,i+1,j,k,6)-uin(1,i+1,j-1,k  ,6))/dy
           jz=half*(jz+jzp1)
           
           flux(1,i,j,k,5,2)=flux(1,i,j,k,5,2) - eta*(jz*bx-jx*bz)*dt/dy
#endif
           
#if NDIM==3
           !3rd direction energy flux
           bx=forth*( uin(1,i  ,j,k,6)+uin(1,i  ,j,k-1,6) + &
                    & uin(1,i+1,j,k,6)+uin(1,i+1,j,k-1,6) )
           by=forth*( uin(1,i,j  ,k,7)+uin(1,i,j  ,k-1,7) + &
                    & uin(1,i,j+1,k,7)+uin(1,i,j+1,k-1,7) )
           
           jx   = (uin(1,i,j  ,k,8)-uin(1,i  ,j-1,k  ,8))/dy &
            &   - (uin(1,i,j  ,k,7)-uin(1,i  ,j  ,k-1,7))/dz
           jxp1 = (uin(1,i,j+1,k,8)-uin(1,i  ,j  ,k  ,8))/dy &
            &   - (uin(1,i,j+1,k,7)-uin(1,i  ,j+1,k-1,7))/dz
           jx   = half*(jx+jxp1)
           jy   = (uin(1,i  ,j,k,6)-uin(1,i  ,j  ,k-1,6))/dz &
            &   - (uin(1,i  ,j,k,8)-uin(1,i-1,j  ,k  ,8))/dx
           jyp1 = (uin(1,i+1,j,k,6)-uin(1,i+1,j  ,k-1,6))/dz &
            &   - (uin(1,i+1,j,k,8)-uin(1,i  ,j  ,k  ,8))/dx
           jy   = half*(jy+jyp1)
           
           flux(1,i,j,k,5,3)=flux(1,i,j,k,5,3) - eta*(jx*by-jy*bx)*dt/dz
#endif

        end do
     end do
  end do
#endif
!
  ! Magnetic field update
  call ct(uin,emfx,emfy,emfz,x,y,z,dx,dy,dz,cartesian,cylindrical,spherical,ngrid)

  ! Energy update :
#if ISO==0
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
#endif
!
  return
end subroutine resistivity
