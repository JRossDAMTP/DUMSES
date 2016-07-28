!###########################################################
!###########################################################
!###########################################################
! uin = (\rho, \rho vx, \rho vy, \rho vz, Etot, bx, by, bz)
#if HDF5==1
subroutine special(uin,x,y,z,time,dt,dx,dy,dz,nspec,mype,npes)
  use hdf5
  use hydro_parameters
  use rdwrt_h5
#if WITHMPI==1
  use mpi_var
#endif
  implicit none

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3) ::uin
  real(dp),dimension(iu1:iu2) ::x
  real(dp),dimension(ju1:ju2) ::y
  real(dp),dimension(ku1:ku2) ::z
  real(dp) ::time,dt,dx,dy,dz
  integer :: nspec,mype,npes

  character(LEN=80) :: filename

  ! HDF5 vars
  integer        :: ierr
  integer(hid_t) :: file_id       ! File ID

  integer, parameter :: nb_real=1
  real(dp), dimension(nb_real) :: para_real

  integer :: nzglob
  real(dp), dimension(:), allocatable :: zglob
  real(dp), dimension(:,:), allocatable :: zuin
  real(dp), dimension(:), allocatable :: vxsq,vysq,vzsq,rhovxvy,rhovx,rhovy,rhovz,maxwell

  if (verbose) write (*,*) 'Entering special...'

  call get_specname(nspec,filename,'h5',mype)

  !Compute dimensions and allocate arrays
#if WITHMPI==1
  nzglob=nz*nzslice
#else
  nzglob=nz
#endif

  allocate(zglob(nzglob),zuin(1:nzglob,1:nvar))
  call get_zprof(uin,z,zuin,zglob,nzglob,npes)

  allocate(vxsq(nzglob),vysq(nzglob),vzsq(nzglob))
  allocate(rhovxvy(nzglob),rhovx(nzglob),rhovy(nzglob),rhovz(nzglob))
  allocate(maxwell(nzglob))

  call get_arrayglob((uin(:,:,:,:,2)/uin(:,:,:,:,1))**2,vxsq,nzglob,npes)
  call get_arrayglob((uin(:,:,:,:,3)/uin(:,:,:,:,1))**2,vysq,nzglob,npes)
  call get_arrayglob((uin(:,:,:,:,4)/uin(:,:,:,:,1))**2,vzsq,nzglob,npes)
  call get_arrayglob(uin(:,:,:,:,2)*uin(:,:,:,:,3)/uin(:,:,:,:,1),rhovxvy,nzglob,npes)
  call get_arrayglob(uin(:,:,:,:,2),rhovx,nzglob,npes)
  call get_arrayglob(uin(:,:,:,:,3),rhovy,nzglob,npes)
  call get_arrayglob(uin(:,:,:,:,4),rhovz,nzglob,npes)
  call get_arrayglob(-uin(:,:,:,:,6)*uin(:,:,:,:,7),maxwell,nzglob,npes)

  if (mype==0) then

     call H5open_f(ierr)
     call H5Fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,ierr)

     !Dataset "para_real"; Type real; dim=1; size=1
     para_real=(/ time /)
     call dump_1d_array_h5(file_id,"para_real",para_real,nb_real)

     !Dataset "yglob"; Type real; dim=1, size=nzglob
     call dump_1d_array_h5(file_id,"zglob",zglob,nzglob)

     !Dataset "yuin"; Type real; dim=2; size=(nzglob,nvar)
     call dump_2d_array_h5(file_id,"zuin",zuin,nzglob,nvar)

     !Dataset "vxsq","vysq","vzsq"; Type real; dim=1; size=nzglob
     call dump_1d_array_h5(file_id,"vxsq",vxsq,nzglob)
     call dump_1d_array_h5(file_id,"vysq",vysq,nzglob)
     call dump_1d_array_h5(file_id,"vzsq",vzsq,nzglob)
     !Dataset "rhovxvz","rhovx","rhovy","rhovz"; Type real; dim=1; size=nzglob
     call dump_1d_array_h5(file_id,"rhovxvy",rhovxvy,nzglob)
     call dump_1d_array_h5(file_id,"rhovx",rhovx,nzglob)
     call dump_1d_array_h5(file_id,"rhovy",rhovy,nzglob)
     call dump_1d_array_h5(file_id,"rhovz",rhovz,nzglob)
     !Dataset "maxwell"; Type real; dim=1; size=nzglob
     call dump_1d_array_h5(file_id,"maxwell",maxwell,nzglob)

     ! Terminate access to HDF5
     call H5Fclose_f(file_id,ierr)
     call H5close_f(ierr)

  endif

  deallocate(zglob,zuin)

  return
end subroutine special
#else
!###########################################################
!###########################################################
!###########################################################
subroutine special(uin,x,y,z,time,dt,dx,dy,dz,nspec,mype,npes)
#if WITHMPI==1
  use mpi_var
#endif
  use hydro_parameters
  implicit none

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar) ::uin
  real(dp),dimension(iu1:iu2) ::x
  real(dp),dimension(ju1:ju2) ::y
  real(dp),dimension(ku1:ku2) ::z
  real(dp) ::time,dt,dx,dy,dz
  integer :: nspec,mype,npes

  return
end subroutine special
#endif
!###########################################################
!###########################################################
!###########################################################
subroutine local_xymean(array,meanarray)
  use amr_parameters
  implicit none
  real(dp), dimension(nx,ny,nz) :: array
  real(dp), dimension(nz) :: meanarray

  integer :: i,j

  meanarray=0.d0
  do j=1,ny
     do i=1,nx
        meanarray(1:nz)=meanarray(1:nz)+array(i,j,1:nz)
     end do
  end do
  meanarray=meanarray/nx/ny

end subroutine local_xymean
!###########################################################
!###########################################################
!###########################################################
subroutine get_arrayglob(array,arrayglob,nzglob,npes)
  use hydro_parameters
#if WITHMPI==1
  use mpi_var
#endif
  implicit none
#if WITHMPI==1
#include "mpif.h"
#endif

  integer :: nzglob,npes
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2) :: array
  real(dp),dimension(1:nzglob) :: arrayglob

  real(dp),dimension(1:nz       ) :: array_loc
  real(dp),dimension(1:nz,1:npes) :: array_gather

  integer :: iproc,k,kglob

#if WITHMPI==1
  integer :: ierr
#endif

  arrayglob=0.d0 ; array_loc=0.d0

  call local_xymean(array(1,1:nx,1:ny,1:nz),array_loc)

#if WITHMPI==1
  
  call MPI_allgather(array_loc   ,nz,MPI_DOUBLE_PRECISION, &
     &               array_gather,nz,MPI_DOUBLE_PRECISION, &
     &               MPI_COMM_WORLD,ierr)

  do iproc=1,npes
     do k=1,nz
        kglob=k+nz*mod((iproc-1)/nxslice, nzslice)
        arrayglob(kglob) = arrayglob(kglob) + array_gather(k,iproc)
     end do
  end do
  arrayglob=arrayglob/nxslice/nyslice
#else
  arrayglob=array_loc
#endif

  return
end subroutine get_arrayglob
!###########################################################
!###########################################################
!###########################################################
subroutine get_zprof(uin,z,zuin,zglob,nzglob,npes)
  use hydro_parameters
#if WITHMPI==1
  use mpi_var
#endif
  implicit none
#if WITHMPI==1
#include "mpif.h"
#endif

  integer :: nzglob,npes
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3) :: uin
  real(dp),dimension(ku1:ku2) :: z
  real(dp),dimension(1:nzglob,1:nvar) :: zuin
  real(dp),dimension(1:nzglob) :: zglob

  real(dp),dimension(1:nz,1:nvar       ) :: zuin_loc
  real(dp),dimension(1:nz,1:nvar,1:npes) :: zuin_gather
  real(dp),dimension(1:nz       ) :: z_loc
  real(dp),dimension(1:nz,1:npes) :: z_gather

  integer :: iproc,i,j,k,kglob

#if WITHMPI==1
  integer :: ierr
#endif

  zuin=0.d0 ; zglob=0.d0

  zuin_loc=0.d0 ; z_loc=0.d0
  z_loc=z(1:nz)
  do j=1,ny
     do i=1,nx
        zuin_loc(1:nz,1:nvar)=zuin_loc(1:nz,1:nvar) + uin(1,i,j,1:nz,1:nvar)
     end do
  end do
  zuin_loc=zuin_loc/nx/ny

#if WITHMPI==1
  
  call MPI_allgather(zuin_loc   ,nz*nvar,MPI_DOUBLE_PRECISION, &
     &               zuin_gather,nz*nvar,MPI_DOUBLE_PRECISION, &
     &               MPI_COMM_WORLD,ierr)
  call MPI_allgather(z_loc   ,nz,MPI_DOUBLE_PRECISION, &
     &               z_gather,nz,MPI_DOUBLE_PRECISION, &
     &               MPI_COMM_WORLD,ierr)

  do iproc=1,npes
     do k=1,nz
        kglob=k+nz*mod((iproc-1)/nxslice, nzslice)
        zuin (kglob,1:nvar) = zuin (kglob,1:nvar) + zuin_gather(k,1:nvar,iproc)
        zglob(kglob)        = zglob(kglob)        + z_gather   (k       ,iproc)
     end do
  end do
  zuin =zuin /nxslice/nyslice
  zglob=zglob/nxslice/nyslice
#else
  zuin =zuin_loc
  zglob=z_loc
#endif

  return
end subroutine get_zprof
!###########################################################
!###########################################################
!###########################################################
subroutine get_specname(nfile,filename,suffix,mype)
  use hydro_parameters
  implicit none

  integer :: nfile,mype
  character(LEN=80) :: filename
  character(LEN=2) :: suffix

  character(LEN=6) :: snumfile
  character(LEN=80) :: filedir,filecmd
  integer :: info

  call convtoasc(nfile,snumfile)
  filename='zprof/zprof_'//trim(snumfile)//'.'//trim(suffix)

  if (nfile==1) then
     filedir='zprof/'
     filecmd='mkdir -p '//trim(filedir)
#ifdef NOSYSTEM
     if (mype==0) call PXFMKDIR(trim(filedir),len(trim(filedir)),O'755',info)
#else
     if (mype==0) call system(filecmd)
#endif
  endif

  return
end subroutine get_specname
