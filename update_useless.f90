!> \file
!! contains subroutine update()
!<

!###########################################################
!###########################################################
!###########################################################
!###########################################################

!  Subroutine UPDATE
!
!> Updates the conservative variables after the mhd step
!! using the computed mhd fluxes.
!<
subroutine update(uin,flux,flux_pre_tot,emfx,emfy,emfz,bval_type,ngrid,dt)
  use amr_parameters
  use hydro_parameters
  use const
  use variables , only : x,y,z,dx,dy,dz,cartesian,cylindrical,spherical
  implicit none

  integer :: ngrid
  real(dp), dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3       ) :: uin 
  real(dp), dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar  ,1:ndim) :: flux
  real(dp), dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2         ,1:ndim) :: flux_pre_tot
  real(dp), dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2                ) :: emfx,emfy,emfz
  integer,dimension(ndim,2) :: bval_type
  real(dp) :: lambda,ratio,dsx,dsy,dfn2,dfn3,dt
  real(dp) :: phi_m, phi_p, phi_c

  integer :: i,j,k,l,ivar,idim,ibx,iby,ibz
  integer :: i0,j0,k0,ihi,jhi,khi,ivit

  real(dp), dimension(3,2,3) :: radius
  real(dp)                   :: xL,xR,xc,sinxL,sinxR,sinxc,yc,dv

  if (verbose) write (*,*) 'Entering update subroutine...'

  if ((cartesian).and.(Omega0>0)) then
     lambda=Omega0*dt ; lambda=2.5D-1*lambda*lambda
     ratio=(1.d0-lambda)/(1.d0+lambda)
  endif

  ibx=6 ; iby=7 ; ibz=8

  ihi=max(1,if2-1)
  jhi=max(1,jf2-1)
  khi=max(1,kf2-1)

  do k=1,khi
     do j=1,jhi
        do i=1,ihi
           radius = 1.
           if (cylindrical .or. spherical) then
              radius(2,1:2,1:3) = x(i)
              radius(2,1  ,1  ) = half*(x(i)+x(i-1))
              radius(2,2  ,1  ) = half*(x(i)+x(i+1))
           endif
           if (spherical) then
              yc = y(j)
              radius(3,1  ,1  ) = half*(x(i)+x(i-1)) * sin(yc)
              radius(3,2  ,1  ) = half*(x(i)+x(i+1)) * sin(yc)
              radius(3,1  ,2  ) = x(i) * sin(half*(y(j)+y(j-1)))
              radius(3,2  ,2  ) = x(i) * sin(half*(y(j)+y(j+1)))
              radius(3,1  ,3  ) = x(i) * sin(yc)
              radius(3,2  ,3  ) = x(i) * sin(yc)
           endif
           if (cartesian)   dv=dx*dy*dz
           if (cylindrical) dv=dx*x(i)*dy*dz
           if (spherical)   dv=dx*x(i)*dy*x(i)*sin(y(j))*dz

           do l=1,ngrid
              do idim=1,ndim

                 i0=0; j0=0; k0=0 
                 if(idim==1)i0=1
                 if(idim==2)j0=1
                 if(idim==3)k0=1

                 !rho
                 uin(l,i,j,k,1)=uin(l,i,j,k,1) + (flux(l,i   ,j   ,k   ,1,idim) &
                      &                        -  flux(l,i+i0,j+j0,k+k0,1,idim))/dv

                 !rho*vert_v
                 ivit=4
                 uin(l,i,j,k,ivit)=uin(l,i,j,k,ivit) + (flux(l,i   ,j   ,k   ,ivit,idim)*radius(ivit-1,1,idim)  &
                      &                              -  flux(l,i+i0,j+j0,k+k0,ivit,idim)*radius(ivit-1,2,idim)) &
                      &                                                    /dv/radius(ivit-1,2,3)

                 !Etot
                 uin(l,i,j,k,5)=uin(l,i,j,k,5) + (flux(l,i   ,j   ,k   ,5,idim) &
                      &                        -  flux(l,i+i0,j+j0,k+k0,5,idim))/dv

                 !Cell-centered B-field (only in 1 & 2D)
                 if (ndim==1) then
                    uin(l,i,j,k,6:nvar)=uin(l,i,j,k,6:nvar)+ &
                         &          (flux(l,i   ,j   ,k   ,6:nvar,idim) &
                         &          -flux(l,i+i0,j+j0,k+k0,6:nvar,idim))/dv
                 else if (ndim==2) then
                    uin(l,i,j,k,nvar)=uin(l,i,j,k,nvar)+ &
                         &          (flux(l,i   ,j   ,k   ,nvar,idim) &
                         &          -flux(l,i+i0,j+j0,k+k0,nvar,idim))/dv
                 endif

              end do

              !Add presure gradient and geometrical terms
              if (ndim>2) uin(l,i,j,k,4) = uin(l,i,j,k,4) + flux_pre_tot(l,i,j,k,3)/dv

              if ((cartesian).and.(Omega0>0.d0)) then

                 !Source term in case of shearing box
                 ! Crank Nicholson for the horizontal momentum source terms
                 dsx= two *Omega0*dt*uin(l,i,j,k,3)/(1.d0+lambda)
                 dsy=-half*Omega0*dt*uin(l,i,j,k,2)/(1.d0+lambda)
                 
                 uin(l,i,j,k,2)=uin(l,i,j,k,2)*ratio + dsx 
                 uin(l,i,j,k,3)=uin(l,i,j,k,3)*ratio + dsy 

                 do idim=1,ndim
                    i0=0; j0=0; k0=0 
                    if(idim==1)i0=1
                    if(idim==2)j0=1
                    if(idim==3)k0=1
                    dfn2=(flux(l,i,j,k,2,idim)-flux(l,i+i0,j+j0,k+k0,2,idim))/dv
                    if (idim==1) dfn2=dfn2+flux_pre_tot(l,i,j,k,1)/dv
                    dfn3=(flux(l,i,j,k,3,idim)-flux(l,i+i0,j+j0,k+k0,3,idim))/dv
                    if (idim==2) dfn3=dfn3+flux_pre_tot(l,i,j,k,2)/dv
                    uin(l,i,j,k,2)=uin(l,i,j,k,2) + dfn2/(1.d0+lambda) + Omega0*dt/(1+lambda)*dfn3
                    uin(l,i,j,k,3)=uin(l,i,j,k,3) + dfn3/(1.d0+lambda) - 2.5D-1*Omega0*dt/(1+lambda)*dfn2
                 end do
                 
                 ! Source terms for the energy equation in stratified shearing box :
                 ! Horizontal potential (gravity+centrifugal) :
                 CALL GET_PHI_SB(phi_c,x(i),z(k))
                 CALL GET_PHI_SB(phi_m,0.5d0*(x(i)+x(i-1)),z(k))
                 CALL GET_PHI_SB(phi_p,0.5d0*(x(i)+x(i+1)),z(k))
                 uin(l,i,j,k,nvar)=uin(l,i,j,k,nvar) - (flux(l,i,j,k,1,1)*(phi_c-phi_m) +  flux(l,i+1,j,k,1,1)*(phi_p-phi_c))/dx/dv
                 ! Vertical gravity :
                 CALL GET_PHI_SB(phi_c,x(i),z(k))
                 CALL GET_PHI_SB(phi_m,x(i),0.5d0*(z(k)+z(k-1)))
                 CALL GET_PHI_SB(phi_p,x(i),0.5d0*(z(k)+z(k+1)))
                 uin(l,i,j,k,nvar)=uin(l,i,j,k,nvar) - (flux(l,i,j,k,1,3)*(phi_c-phi_m) +  flux(l,i+1,j,k,1,3)*(phi_p-phi_c))/dz/dv

              else

                 do idim=1,ndim
                    i0=0; j0=0; k0=0 
                    if(idim==1)i0=1
                    if(idim==2)j0=1
                    if(idim==3)k0=1
                    do ivit=2,3
                       uin(l,i,j,k,ivit)=uin(l,i,j,k,ivit) + (flux(l,i   ,j   ,k   ,ivit,idim)*radius(ivit-1,1,idim)  &
                            &                              -  flux(l,i+i0,j+j0,k+k0,ivit,idim)*radius(ivit-1,2,idim)) &
                            &                                                    /dv/radius(ivit-1,2,3)
                    end do
                 end do
                 !Add presure gradient and geometrical terms
                             uin(l,i,j,k,2) = uin(l,i,j,k,2) + flux_pre_tot(l,i,j,k,1)/dv
                 if (ndim>1) uin(l,i,j,k,3) = uin(l,i,j,k,3) + flux_pre_tot(l,i,j,k,2)/dv

              endif

           end do
        end do
     end do
  end do

  call ct(uin,emfx,emfy,emfz,x,y,z,dx,dy,dz,cartesian,cylindrical,spherical,ngrid)

  if (verbose) write (*,*) 'End of update subroutine'

  return
end subroutine update
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine conservative_update(uin,flux,bval_type,ngrid)
  use amr_parameters
  use hydro_parameters
  implicit none

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin 
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  integer,dimension(ndim,2) :: bval_type
  integer ::ngrid

  integer i,j,k,l,ivar,idim
  integer i0,j0,k0

  if (verbose) write (*,*) 'Entering conservative_update...'

  do idim=1,ndim
     i0=0; j0=0; k0=0 
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1

     do l=1,ngrid
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 uin(l,i,j,k,1:nvar)=uin(l,i,j,k,1:nvar)+ &
                      & (flux(l,i   ,j   ,k   ,1:nvar,idim) &
                      & -flux(l,i+i0,j+j0,k+k0,1:nvar,idim))
              end do
           end do
        end do
     end do

  end do
!
! Need to explicitely update Bx at the outer edge of the grid
! to conserve divB in case of shearing boundary conditions.
!
#if NDIM==3
  if (bval_type(1,2)==3) then

     do l=1,ngrid
        do k=1,nz
           do j=1,ny
              uin(l,nx+1,j,k,6)=uin(l,nx+1,j,k,6) +&
                      & (flux(l,nx+1,j  ,k  ,6,2)  &
                      & -flux(l,nx+1,j+1,k  ,6,2))+&
                      & (flux(l,nx+1,j  ,k  ,6,3)  &
                      & -flux(l,nx+1,j  ,k+1,6,3))
           end do
        end do
     end do
  endif 
  if (bval_type(3,2)/=1) then
     do l=1,ngrid
        do j=1,ny
           do i=1,nx
              do idim=1,2
                 i0=0; j0=0
                 if(idim==1)i0=1
                 if(idim==2)j0=1
                 uin(l,i,j,nz+1,8)=uin(l,i,j,nz+1,8) + &
                      & (flux(l,i   ,j   ,nz+1,8,idim) &
                      & -flux(l,i+i0,j+j0,nz+1,8,idim))
              end do
           end do
        end do
     enddo
  endif
#endif
  
  !Right faces B-field
  uin(1,iu1:iu2-1,ju1:ju2  ,ku1:ku2  ,nvar+1)=uin(1,iu1+1:iu2,ju1  :ju2,ku1  :ku2,6)
  uin(1,iu1:iu2  ,ju1:ju2-1,ku1:ku2  ,nvar+2)=uin(1,iu1  :iu2,ju1+1:ju2,ku1  :ku2,7)
  uin(1,iu1:iu2  ,ju1:ju2  ,ku1:ku2-1,nvar+3)=uin(1,iu1  :iu2,ju1  :ju2,ku1+1:ku2,8)

  return
end subroutine conservative_update
