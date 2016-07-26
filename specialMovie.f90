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

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar) ::uin
  real(dp),dimension(iu1:iu2) ::x
  real(dp),dimension(ju1:ju2) ::y
  real(dp),dimension(ku1:ku2) ::z
  real(dp) ::time,dt,dx,dy,dz
  integer :: nspec,mype,npes

  character(LEN=80) :: filename

  ! HDF5 vars
  integer        :: ierr
  integer(hid_t) :: file_id       ! File ID

  integer, parameter :: nb_real=5,nb_int=7,nb_mpi=9

  integer , dimension(nb_int ) :: para_int
  integer , dimension(nb_mpi ) :: para_mpi
  real(dp), dimension(nb_real) :: para_real

  integer :: n1,n2,n3

  if (verbose) write (*,*) 'Entering output (hdf5)...'

  call get_specname(nspec,mype,filename)

  call H5open_f(ierr)
  call H5Fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,ierr)

  !Dataset "para_real"; Type real; dim=1; size=nb_real
  para_real=(/ time,dt,dx,dy,dz /)
  call dump_1d_array_h5(file_id,"para_real",para_real,nb_real)

  !Dataset "para_int"; Type integer; dim=1; size=nb_int
#if WITHMPI==1
  para_int=(/ nspec,nx,ny,nz,nxslice,nyslice,nzslice/)
#else
  para_int=(/ nspec,nx,ny,nz,1,1,1/)
#endif
  call dump_1d_array_int_h5(file_id,"para_int",para_int,nb_int)

  !Dataset "x"; Type real; dim=1, size=n1
  n1=iu2-iu1+1
  call dump_1d_array_h5(file_id,"x",x,n1)

  !Dataset "y"; Type real; dim=1, size=n2
  n2=ju2-ju1+1
  call dump_1d_array_h5(file_id,"y",y,n2)

  !Dataset "z"; Type real; dim=1, size=n3
  n3=ku2-ku1+1
  call dump_1d_array_h5(file_id,"z",z,n3)

  !Dataset rho; Type real; dim=3; size=(n1,n2,n3)
  call dump_3d_array_h5(file_id,"rho",uin(1,:,:,:,1),n1,n2,n3)

#if WITHMPI==1
  !Dataset "para_mpi"; Type integer; dim=1; size=nb_mpi
  para_mpi=(/ xleft,xright,yleft,yright,zleft,zright,&
            & xposition,yposition,zposition/)
  call dump_1d_array_int_h5(file_id,"para_mpi",para_mpi,nb_mpi)
#endif
  
  ! Terminate access to HDF5
  call H5Fclose_f(file_id,ierr)
  call H5close_f(ierr)

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
subroutine get_specname(nspec,mype,filename)
  use hydro_parameters
  implicit none
#if WITHMPI==1
#include "mpif.h"
#endif

  integer :: nspec,mype

  integer :: info
  character(LEN=6) :: snumfile,n_mype
  character(LEN=80) :: filename,filedir,filecmd
#if WITHMPI==1
  integer :: ierr
#endif

  call convtoasc(nspec,snumfile)
  call convtoasc(mype,n_mype  )
  filedir='movie_'//trim(snumfile)//'/'
  filecmd='mkdir -p '//trim(filedir)
#ifdef NOSYSTEM
  if (mype==0) call PXFMKDIR(trim(filedir),len(trim(filedir)),O'755',info)
#else
  if (mype==0) call system(filecmd)
#endif
#if WITHMPI==1
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
  filename=trim(filedir)//'rho.'//trim(n_mype)
  
  return
end subroutine get_specname
