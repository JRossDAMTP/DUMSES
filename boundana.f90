subroutine xinner_ana
  use hydro_parameters
  use variables
  implicit none
!
  integer :: j,k
!
  ! Periodic BC
  do k=ku1,ku2
     do j=ju1,ju2      
        uin(1,iu1  ,j,k,:) = uin(1,nx-2,j,k,:)
        uin(1,iu1+1,j,k,:) = uin(1,nx-1,j,k,:)
        uin(1,iu1+2,j,k,:) = uin(1,nx  ,j,k,:)
     end do
  end do
!
  return
end subroutine xinner_ana
!###########################################################
!###########################################################
!###########################################################
subroutine xouter_ana
  use hydro_parameters
  use variables
  implicit none

  integer ::k,j

  ! Periodic BC
  do k=ku1,ku2
     do j=ju1,ju2
        uin(1,iu2-2,j,k,:) = uin(1,1 ,j,k,:)
        uin(1,iu2-1,j,k,:) = uin(1,2 ,j,k,:)
        uin(1,iu2  ,j,k,:) = uin(1,3 ,j,k,:)
     end do
  end do
!
  return
end subroutine xouter_ana
!###########################################################
!###########################################################
!###########################################################
subroutine yinner_ana
#if NDIM>1
  use hydro_parameters
  use variables
  implicit none

  integer ::k,i

  ! Periodic BC
  do k=ku1,ku2
     do i=iu1,iu2
        uin(1,i,ju1  ,k,:) = uin(1,i,ny-2,k,:)
        uin(1,i,ju1+1,k,:) = uin(1,i,ny-1,k,:)
        uin(1,i,ju1+2,k,:) = uin(1,i,ny  ,k,:)
     end do
  end do
!
#endif
  return
end subroutine yinner_ana
!###########################################################
!###########################################################
!###########################################################
subroutine youter_ana
#if NDIM>1
  use hydro_parameters
  use variables
  implicit none

  integer ::k,i

  do k=ku1,ku2
     do i=iu1,iu2
        uin(1,i,ju2-2,k,:) = uin(1,i,1 ,k,:)
        uin(1,i,ju2-1,k,:) = uin(1,i,2 ,k,:)
        uin(1,i,ju2  ,k,:) = uin(1,i,3 ,k,:)
     end do
  end do
!
#endif
  return
end subroutine youter_ana
!###########################################################
!###########################################################
!###########################################################
subroutine zinner_ana
#if NDIM>1
  use hydro_parameters
  use const
  use variables
  use stratified
  implicit none

  integer ::i,j,k
  real(dp)::H,factor,ratio_nyp1,ratio_nyp2,ratio_nyp3,dbxdx,dbydy
  real(dp)::rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,csq

  H=ciso/Omega0
  factor=-dz/2.d0/H/H
  if (floor) then
     ratio_nyp1=1.d0
     ratio_nyp2=1.d0
     ratio_nyp3=1.d0
  else
     ratio_nyp1=exp(factor*(-2*(zmin+half*dz)+     dz))
     ratio_nyp2=exp(factor*(-2*(zmin+half*dz)+3.d0*dz))
     ratio_nyp3=exp(factor*(-2*(zmin+half*dz)+5.d0*dz))
  endif
  do j=ju1,ju2-1
     do i=iu1,iu2-1
#if ISO==0
        rho=uin(1,i,j,ku1+3,1) ; E=uin(1,i,j,ku1+3,5)
        rhovx=uin(1,i,j,ku1+3,2) ; rhovy=uin(1,i,j,ku1+3,3) ; rhovz=uin(1,i,j,ku1+3,4)
        Bxc=half*(uin(1,i,j,ku1+3,6)+uin(1,i+1,j  ,ku1+3,6))
        Byc=half*(uin(1,i,j,ku1+3,7)+uin(1,i  ,j+1,ku1+3,7))
        Bzc=half*(uin(1,i,j,ku1+3,8)+uin(1,i  ,j  ,ku1+4,8))
          !write (*,*) i,j,Bxc,Byc,Bzc
        call get_csq(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,csq)
        H=sqrt(csq)/Omega0/sqrt(gamma)
        factor=-dz/2.d0/H/H
        ratio_nyp1=exp(factor*(-2*(zmin+half*dz)+     dz))
        ratio_nyp2=exp(factor*(-2*(zmin+half*dz)+3.d0*dz))
        ratio_nyp3=exp(factor*(-2*(zmin+half*dz)+5.d0*dz))
#endif
        !extrapolating density
        factor=(log(uin(1,i,j,ku1+4,1))-log(uin(1,i,j,ku1+3,1)))/log(uin(1,i,j,ku1+3,1))
        uin(1,i,j,ku1+2,1) = uin(1,i,j,ku1+3,1)**(1.d0-factor)
        factor=(log(uin(1,i,j,ku1+3,1))-log(uin(1,i,j,ku1+2,1)))/log(uin(1,i,j,ku1+2,1))
        uin(1,i,j,ku1+1,1) = uin(1,i,j,ku1+2,1)**(1.d0-factor)
        factor=(log(uin(1,i,j,ku1+2,1))-log(uin(1,i,j,ku1+1,1)))/log(uin(1,i,j,ku1+1,1))
        uin(1,i,j,ku1  ,1) = uin(1,i,j,ku1+1,1)**(1.d0-factor)
        !uin(1,i,j,ku1+2,1) = max(uin(1,i,j,ku1+3,1)*ratio_nyp1,smallr)
        !uin(1,i,j,ku1+1,1) = max(uin(1,i,j,ku1+2,1)*ratio_nyp2,smallr)
        !uin(1,i,j,ku1  ,1) = max(uin(1,i,j,ku1+1,1)*ratio_nyp3,smallr)
     end do
  end do
  do j=ju1,ju2
     do i=iu1,iu2
        !zero gradient BC for velocities
        uin(1,i,j,ku1  ,2:4) = uin(1,i,j,ku1+3,2:4)/uin(1,i,j,ku1+3,1)
        uin(1,i,j,ku1+1,2:4) = uin(1,i,j,ku1+3,2:4)/uin(1,i,j,ku1+3,1)
        uin(1,i,j,ku1+2,2:4) = uin(1,i,j,ku1+3,2:4)/uin(1,i,j,ku1+3,1)
        uin(1,i,j,ku1  ,2:3) = uin(1,i,j,ku1  ,2:3)*uin(1,i,j,ku1  ,1)
        uin(1,i,j,ku1+1,2:3) = uin(1,i,j,ku1+1,2:3)*uin(1,i,j,ku1+1,1)
        uin(1,i,j,ku1+2,2:3) = uin(1,i,j,ku1+2,2:3)*uin(1,i,j,ku1+2,1)
        uin(1,i,j,ku1  ,4) = min(uin(1,i,j,ku1  ,4)*uin(1,i,j,ku1  ,1),0.d0)
        uin(1,i,j,ku1+1,4) = min(uin(1,i,j,ku1+1,4)*uin(1,i,j,ku1+1,1),0.d0)
        uin(1,i,j,ku1+2,4) = min(uin(1,i,j,ku1+2,4)*uin(1,i,j,ku1+2,1),0.d0)
        !uin(1,i,j,ku1  ,4) = min(uin(1,i,j,ku1+3,4),0.d0)
        !uin(1,i,j,ku1+1,4) = min(uin(1,i,j,ku1+3,4),0.d0)
        !uin(1,i,j,ku1+2,4) = min(uin(1,i,j,ku1+3,4),0.d0)
        !uin(1,i,j,ku1  ,4) = uin(1,i,j,ku1+3,4)
        !uin(1,i,j,ku1+1,4) = uin(1,i,j,ku1+3,4)
        !uin(1,i,j,ku1+2,4) = uin(1,i,j,ku1+3,4)
        !Tangential magnetic field
        uin(1,i,j,ku1  ,6:7) = 0. !Vanishing tangential field
        uin(1,i,j,ku1+1,6:7) = 0. !Vanishing tangential field
        uin(1,i,j,ku1+2,6:7) = 0. !Vanishing tangential field
        !Tangential magnetic field
        !uin(1,i,j,ku1  ,6:7) = uin(1,i,j,ku1+5,6:7)
        !uin(1,i,j,ku1+1,6:7) = uin(1,i,j,ku1+4,6:7)
        !uin(1,i,j,ku1+2,6:7) = uin(1,i,j,ku1+3,6:7)
!        uin(1,i,j,ku1  ,6:7) = uin(1,i,j,ku1+3,6:7)*(uin(1,i,j,ku1  ,1)/uin(1,i,j,ku1+3,1))**half
!        uin(1,i,j,ku1+1,6:7) = uin(1,i,j,ku1+3,6:7)*(uin(1,i,j,ku1+1,1)/uin(1,i,j,ku1+3,1))**half
!        uin(1,i,j,ku1+2,6:7) = uin(1,i,j,ku1+3,6:7)*(uin(1,i,j,ku1+2,1)/uin(1,i,j,ku1+3,1))**half
     end do
  end do
  do j=ju1,ju2-1
     do i=iu1,iu2-1
        !Normal magnetic field
        dbxdx=(uin(1,i+1,j  ,ku1+2,6)-uin(1,i,j,ku1+2,6))/dx
        dbydy=(uin(1,i  ,j+1,ku1+2,7)-uin(1,i,j,ku1+2,7))/dy
        uin(1,i,j,ku1+2,8) = uin(1,i,j,ku1+3,8) + dz*(dbxdx+dbydy)
        dbxdx=(uin(1,i+1,j  ,ku1+1,6)-uin(1,i,j,ku1+1,6))/dx
        dbydy=(uin(1,i  ,j+1,ku1+1,7)-uin(1,i,j,ku1+1,7))/dy
        uin(1,i,j,ku1+1,8) = uin(1,i,j,ku1+2,8) + dz*(dbxdx+dbydy)
        dbxdx=(uin(1,i+1,j  ,ku1  ,6)-uin(1,i,j,ku1  ,6))/dx
        dbydy=(uin(1,i  ,j+1,ku1  ,7)-uin(1,i,j,ku1  ,7))/dy
        uin(1,i,j,ku1  ,8) = uin(1,i,j,ku1+1,8) + dz*(dbxdx+dbydy)
     end do
  end do
!
#if ISO==0
  !Zero gradient Temperature
  do j=ju1,ju2-1
     do i=iu1,iu2-1
        rho=uin(1,i,j,ku1+3,1) ; E=uin(1,i,j,ku1+3,5)
        rhovx=uin(1,i,j,ku1+3,2) ; rhovy=uin(1,i,j,ku1+3,3) ; rhovz=uin(1,i,j,ku1+3,4)
        Bxc=half*(uin(1,i,j,ku1+3,6)+uin(1,i+1,j  ,ku1+3,6))
        Byc=half*(uin(1,i,j,ku1+3,7)+uin(1,i  ,j+1,ku1+3,7))
        Bzc=half*(uin(1,i,j,ku1+3,8)+uin(1,i  ,j  ,ku1+4,8))
        call get_csq(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,csq)
	csq=1.d0!! overide set to constant value
        do k=ku1+2,ku1,-1
           rho=uin(1,i,j,k,1)
           rhovx=uin(1,i,j,k,2) ; rhovy=uin(1,i,j,k,3) ; rhovz=uin(1,i,j,k,4)
           Bxc=half*(uin(1,i,j,k,6)+uin(1,i+1,j  ,k  ,6))
           Byc=half*(uin(1,i,j,k,7)+uin(1,i  ,j+1,k  ,7))
           Bzc=half*(uin(1,i,j,k,8)+uin(1,i  ,j  ,k+1,8))
           call get_E(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,csq)
           uin(1,i,j,k,5)=E!uin(1,i,j,ku1+3,5)!E
        end do
     end do
  end do
#endif
!
#endif
  return
end subroutine zinner_ana
!###########################################################
!###########################################################
!###########################################################
subroutine zouter_ana
#if NDIM>1
  use hydro_parameters
  use const
  use variables
  use stratified
  implicit none

  integer ::i,j,k
  real(dp)::H,factor,ratio_nyp1,ratio_nyp2,ratio_nyp3,dbxdx,dbydy
  real(dp)::rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,csq

  H=ciso/Omega0
  factor=-dz/2.d0/H/H
  if (floor) then
     ratio_nyp1=1.d0
     ratio_nyp2=1.d0
     ratio_nyp3=1.d0
  else
     ratio_nyp1=exp(factor*(2*(zmax-half*dz)+     dz))
     ratio_nyp2=exp(factor*(2*(zmax-half*dz)+3.d0*dz))
     ratio_nyp3=exp(factor*(2*(zmax-half*dz)+5.d0*dz))
  endif
  do j=ju1,ju2-1
     do i=iu1,iu2-1
#if ISO==0
        rho=uin(1,i,j,ku2-3,1) ; E=uin(1,i,j,ku2-3,5)
        rhovx=uin(1,i,j,ku2-3,2) ; rhovy=uin(1,i,j,ku2-3,3) ; rhovz=uin(1,i,j,ku2-3,4)
        Bxc=half*(uin(1,i,j,ku2-3,6)+uin(1,i+1,j  ,ku2-3,6))
        Byc=half*(uin(1,i,j,ku2-3,7)+uin(1,i  ,j+1,ku2-3,7))
        Bzc=half*(uin(1,i,j,ku2-3,8)+uin(1,i  ,j  ,ku2-2,8))
        call get_csq(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,csq)
        H=sqrt(csq)/Omega0/sqrt(gamma)
        factor=-dz/2.d0/H/H
        ratio_nyp1=exp(factor*(2*(zmax-half*dz)+     dz))
        ratio_nyp2=exp(factor*(2*(zmax-half*dz)+3.d0*dz))
        ratio_nyp3=exp(factor*(2*(zmax-half*dz)+5.d0*dz))
#endif
        !extrapolating density
        factor=(log(uin(1,i,j,ku2-3,1))-log(uin(1,i,j,ku2-4,1)))/log(uin(1,i,j,ku2-3,1))
        uin(1,i,j,ku2-2,1) = uin(1,i,j,ku2-3,1)**(1.d0+factor)
        factor=(log(uin(1,i,j,ku2-2,1))-log(uin(1,i,j,ku2-3,1)))/log(uin(1,i,j,ku2-2,1))
        uin(1,i,j,ku2-1,1) = uin(1,i,j,ku2-2,1)**(1.d0+factor)
        factor=(log(uin(1,i,j,ku2-1,1))-log(uin(1,i,j,ku2-2,1)))/log(uin(1,i,j,ku2-1,1))
        uin(1,i,j,ku2  ,1) = uin(1,i,j,ku2-1,1)**(1.d0+factor)
!        uin(1,i,j,ku2-2,1) = max(uin(1,i,j,ku2-3,1)*ratio_nyp1,smallr)
!        uin(1,i,j,ku2-1,1) = max(uin(1,i,j,ku2-2,1)*ratio_nyp2,smallr)
!        uin(1,i,j,ku2  ,1) = max(uin(1,i,j,ku2-1,1)*ratio_nyp3,smallr)
     end do
  end do
  do j=ju1,ju2
     do i=iu1,iu2
        !zero gradient BC for velocities
        uin(1,i,j,ku2-2,2:4) = uin(1,i,j,ku2-3,2:4)/uin(1,i,j,ku2-3,1)
        uin(1,i,j,ku2-1,2:4) = uin(1,i,j,ku2-3,2:4)/uin(1,i,j,ku2-3,1)
        uin(1,i,j,ku2  ,2:4) = uin(1,i,j,ku2-3,2:4)/uin(1,i,j,ku2-3,1)
        uin(1,i,j,ku2-2,2:3) = uin(1,i,j,ku2-2,2:3)*uin(1,i,j,ku2-2,1)
        uin(1,i,j,ku2-1,2:3) = uin(1,i,j,ku2-1,2:3)*uin(1,i,j,ku2-1,1)
        uin(1,i,j,ku2  ,2:3) = uin(1,i,j,ku2  ,2:3)*uin(1,i,j,ku2  ,1)
        uin(1,i,j,ku2-2,4) = max(uin(1,i,j,ku2-2,4)*uin(1,i,j,ku2-2,1),0.d0)
        uin(1,i,j,ku2-1,4) = max(uin(1,i,j,ku2-1,4)*uin(1,i,j,ku2-1,1),0.d0)
        uin(1,i,j,ku2  ,4) = max(uin(1,i,j,ku2  ,4)*uin(1,i,j,ku2  ,1),0.d0)
        !uin(1,i,j,ku2-2,4) = max(uin(1,i,j,ku2-3,4),0.d0)
        !uin(1,i,j,ku2-1,4) = max(uin(1,i,j,ku2-3,4),0.d0)
        !uin(1,i,j,ku2  ,4) = max(uin(1,i,j,ku2-3,4),0.d0)
        !uin(1,i,j,ku2-2,4) = uin(1,i,j,ku2-3,4)
        !uin(1,i,j,ku2-1,4) = uin(1,i,j,ku2-3,4)
        !uin(1,i,j,ku2  ,4) = uin(1,i,j,ku2-3,4)
        !Tangential magnetic field
        !uin(1,i,j,ku2-2,6:7) = uin(1,i,j,ku2-3,6:7)
        !uin(1,i,j,ku2-1,6:7) = uin(1,i,j,ku2-3,6:7)
        !uin(1,i,j,ku2  ,6:7) = uin(1,i,j,ku2-3,6:7)
        !Tangential magnetic field
        uin(1,i,j,ku2-2,6:7) = 0.!- uin(1,i,j,ku2-3,6:7)
        uin(1,i,j,ku2-1,6:7) = 0.!- uin(1,i,j,ku2-4,6:7)
        uin(1,i,j,ku2  ,6:7) = 0.!- uin(1,i,j,ku2-5,6:7)
        !Tangential magnetic field
        !uin(1,i,j,ku2-2,6:7) = uin(1,i,j,ku2-3,6:7)*(uin(1,i,j,ku2-2,1)/uin(1,i,j,ku2-3,1))**0.5
        !uin(1,i,j,ku2-1,6:7) = uin(1,i,j,ku2-3,6:7)*(uin(1,i,j,ku2-1,1)/uin(1,i,j,ku2-3,1))**0.5
        !uin(1,i,j,ku2  ,6:7) = uin(1,i,j,ku2-3,6:7)*(uin(1,i,j,ku2  ,1)/uin(1,i,j,ku2-3,1))**0.5
     end do
  end do
  do j=ju1,ju2-1
     do i=iu1,iu2-1
        !Normal magnetic field
        dbxdx=(uin(1,i+1,j,ku2-2,6)-uin(1,i,j,ku2-2,6))/dx
        dbydy=(uin(1,i  ,j,ku2-2,7)-uin(1,i,j,ku2-2,7))/dy
        uin(1,i,j,ku2-1,8) = uin(1,i,j,ku2-2,8)-dz*(dbxdx+dbydy)
        dbxdx=(uin(1,i+1,j,ku2-1,6)-uin(1,i,j,ku2-1,6))/dx
        dbydy=(uin(1,i  ,j,ku2-1,7)-uin(1,i,j,ku2-1,7))/dy
        uin(1,i,j,ku2  ,8) = uin(1,i,j,ku2-1,8)-dz*(dbxdx+dbydy)
     end do
  end do
!
#if ISO==0
  !Zero gradient Temperature
  do j=ju1,ju2-1
     do i=iu1,iu2-1
        rho=uin(1,i,j,ku2-3,1) ; E=uin(1,i,j,ku2-3,5)
        rhovx=uin(1,i,j,ku2-3,2) ; rhovy=uin(1,i,j,ku2-3,3) ; rhovz=uin(1,i,j,ku2-3,4)
        Bxc=half*(uin(1,i,j,ku2-3,6)+uin(1,i+1,j  ,ku2-3,6))
        Byc=half*(uin(1,i,j,ku2-3,7)+uin(1,i  ,j+1,ku2-3,7))
        Bzc=half*(uin(1,i,j,ku2-3,8)+uin(1,i  ,j  ,ku2-2,8))
        call get_csq(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,csq)
	csq=1.d0!! overide set to constant value
        do k=ku2-2,ku2-1
           rho=uin(1,i,j,k,1)
           rhovx=uin(1,i,j,k,2) ; rhovy=uin(1,i,j,k,3) ; rhovz=uin(1,i,j,k,4)
           Bxc=half*(uin(1,i,j,k,6)+uin(1,i+1,j  ,k  ,6))
           Byc=half*(uin(1,i,j,k,7)+uin(1,i  ,j+1,k  ,7))
           Bzc=half*(uin(1,i,j,k,8)+uin(1,i  ,j  ,k+1,8))
           call get_E(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,csq)
           uin(1,i,j,k,5)=E!uin(1,i,j,ku2-3,5)!E
        end do
     end do
  end do
#endif
!
#endif
  return
end subroutine zouter_ana
