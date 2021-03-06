!###########################################################
!###########################################################
!###########################################################
subroutine history(uin,x,z,time,dt,dx,dy,dz,mype,npes)
  use hydro_parameters
  implicit none
#if WITHMPI==1
#include "mpif.h"
#endif

  integer :: mype,npes
  real(dp), dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3) :: uin
  real(dp),dimension(iu1:iu2) ::x
  real(dp),dimension(ku1:ku2) ::z
  real(dp) ::time,dt,dx,dy,dz

  real(dp)::mass,maxwell,reynolds,dtau,magp,divB,gasp, Ec, Em, P_mid,num_mid,dM, dR, Hs_num, Hs_den,Hs
  real(dp),dimension(ndim)::mean_B,rms_B,rms_v
  real(dp),dimension(:,:,:,:),allocatable :: loc_var
  real(dp),dimension(:,:,:  ),allocatable :: loc_var_z_mean
  character(LEN=255) :: filename='history.txt'
  integer ::i,j,k,nvar_mean,ipos

  if (verbose) write (*,*) 'Entering history...'

#if NDIM==1
  ipos=1
  do while (x(ipos)<0.5)
     ipos=ipos+1
  end do
  if (mype==0) then
     if (time==0.d0) then
        open(unit=2,file=filename,status='unknown')
     else
        open(unit=2,file=filename,status='old',position='append')
     endif
     write(2,110) time,uin(1,ipos,1,1,1),uin(1,ipos,1,1,2),uin(1,ipos,1,1,5)
     close(2)
  endif
110  format  (1x,1pe12.5,1x,3e14.5)  
#endif

#if NDIM==3

  mass = 0.d0 ; magp = 0.d0 ; maxwell = 0.d0
  mean_B = 0.d0 ; rms_B=0.d0 ; rms_v=0.d0
  gasp=0.d0; P_mid=0.d0; num_mid=0.d0
  Hs=0.d0; Hs_num=0.d0; Hs_den=0.d0; dM=0.d0; dR=0
  do k=1,nz
      do j=1,ny
         do i=1,nx
            mass = mass + uin(1,i,j,k,1)
            magp = magp + &
                 & 2.5d-1*(uin(1,i,j,k,6)+uin(1,i+1,j  ,k  ,6))**2 + &
                 & 2.5d-1*(uin(1,i,j,k,7)+uin(1,i  ,j+1,k  ,7))**2 + &
                 & 2.5d-1*(uin(1,i,j,k,8)+uin(1,i  ,j  ,k+1,8))**2
            maxwell = maxwell - &
                 & 5.d-1*(uin(1,i,j,k,6)+uin(1,i+1,j,k  ,6)) * &
                 & 5.d-1*(uin(1,i,j,k,7)+uin(1,i  ,j,k+1,7))
            mean_B(1) = mean_B(1) + uin(1,i,j,k,6)
            mean_B(2) = mean_B(2) + uin(1,i,j,k,7)
            mean_B(3) = mean_B(3) + uin(1,i,j,k,8)
            rms_B(1)  = rms_B(1)  + uin(1,i,j,k,6)**2
            rms_B(2)  = rms_B(2)  + uin(1,i,j,k,7)**2
            rms_B(3)  = rms_B(3)  + uin(1,i,j,k,8)**2
            rms_v(1)  = rms_v(1)  + (uin(1,i,j,k,2)/uin(1,i,j,k,1))**2
            rms_v(2)  = rms_v(2)  + (uin(1,i,j,k,3)/uin(1,i,j,k,1))**2
            rms_v(3)  = rms_v(3)  + (uin(1,i,j,k,4)/uin(1,i,j,k,1))**2

 	    Ec=0.5*(uin(1,i,j,k,2)**2+uin(1,i,j,k,3)**2+uin(1,i,j,k,4)**2)/uin(1,i,j,k,1)
  	    Em=0.5*(uin(1,i,j,k,6)**2+uin(1,i,j,k,7)**2+uin(1,i,j,k,8)**2)

            gasp = gasp + (gamma-1.d0)*(uin(1,i,j,k,5)-Ec-Em)
	 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   if (abs(z(k))<1.1*dz) then
		P_mid=P_mid+(gamma-1.d0)*(uin(1,i,j,k,5)-Ec-Em)
		num_mid=num_mid+1.0
	   endif


	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         end do
      end do
  end do
  dtau    = dx*dy*dz/(xmax-xmin)/(ymax-ymin)/(zmax-zmin)
  magp    = magp*dtau
  mass    = mass*dtau
  maxwell = maxwell*dtau
  mean_B  = mean_B*dtau
  rms_B   = rms_B*dtau
  rms_v   = rms_v*dtau
  gasp    = gasp*dtau


  reynolds = 0.d0

  nvar_mean=3
  allocate(loc_var       (iu1:iu2,ju1:ju2,ku1:ku2,nvar_mean))
  allocate(loc_var_z_mean(iu1:iu2,ju1:ju2        ,nvar_mean))
  loc_var(:,:,:,1)=uin(1,:,:,:,1)
  loc_var(:,:,:,2)=uin(1,:,:,:,2)/uin(1,:,:,:,1)
  loc_var(:,:,:,3)=uin(1,:,:,:,3)/uin(1,:,:,:,1)
  call compute_z_mean(loc_var,loc_var_z_mean,nvar_mean,mype,npes)
  deallocate(loc_var)

  do k=1,nz
      do j=1,ny
         do i=1,nx
!            reynolds=reynolds+uin(1,i,j,k,1)*&
!                    & (uin(1,i,j,k,2)/uin(1,i,j,k,1)-loc_var_z_mean(i,j,2))*&
!                    & (uin(1,i,j,k,4)/uin(1,i,j,k,1)-loc_var_z_mean(i,j,3))
	    dR=uin(1,i,j,k,1)*(uin(1,i,j,k,2)/uin(1,i,j,k,1))*&
                    &                        (uin(1,i,j,k,3)/uin(1,i,j,k,1))

            reynolds=reynolds+dR
	    dM=-  5.d-1*(uin(1,i,j,k,6)+uin(1,i+1,j,k  ,6)) * &
                 & 5.d-1*(uin(1,i,j,k,7)+uin(1,i  ,j,k+1,7))

	    Hs_num=Hs_num+(dR+dM)*abs(z(k))
	    Hs_den=Hs_den+(dR+dM)



         end do
      end do
  end do
  reynolds=reynolds*dtau
  deallocate(loc_var_z_mean)

  divB=0.d0
  do k=1,nz
      do j=1,ny
         do i=1,nx
            divB=divB + &
              (uin(1,i+1,j  ,k  ,6)-uin(1,i,j,k,6))/dx + &
              (uin(1,i  ,j+1,k  ,7)-uin(1,i,j,k,7))/dy + &
              (uin(1,i  ,j  ,k+1,8)-uin(1,i,j,k,8))/dz
         end do
      end do
  end do


#if WITHMPI==1
  call vol_average(mass,1)
  call vol_average(reynolds,1)
  call vol_average(maxwell,1)
  call vol_average(magp,1)
  call vol_average(divB,1)
  call vol_average(mean_B,3)
  call vol_average(rms_B ,3)
  call vol_average(rms_v ,3)
  call vol_average(gasp,1)
  call vol_average(Hs_den,1)
  call vol_average(Hs_num,1)
#endif   

Hs=Hs_num/Hs_den

#if WITHMPI==1
     call sumBcast(P_mid,1)
     call sumBcast(num_mid,1)
     P_mid=P_mid/num_mid
#endif   
      
  if (mype==0) then
     if (time==0.d0) then
        open(unit=2,file=filename,status='unknown')
     else
        open(unit=2,file=filename,status='old',position='append')
     endif
     write(2,120) time,dt,mass,maxwell,reynolds,maxwell+reynolds,magp,&
                & mean_B(1),mean_B(2),mean_B(3), gasp,P_mid,rms_B(2),rms_B(3),&
                Hs,rms_v(2),rms_v(3)!,divB
     close(2)
  endif
120  format  (1x,1pe12.5,1x,1pe12.5,15e14.5)
#endif

end subroutine history
