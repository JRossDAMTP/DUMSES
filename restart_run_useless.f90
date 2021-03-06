#if HDF5==1
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine restart_run(uin,x,y,z,time,dt,dx,dy,dz,ndump,nhist,nspec,mype)
  use hdf5
#if WITHMPI==1
  use mpi_var
#endif
  use hydro_parameters
  use rdwrt_h5
  use stratified , only : mass0
  implicit none
  real(dp)::time,dt,dx,dy,dz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin
  real(dp),dimension(iu1:iu2)::x
  real(dp),dimension(ju1:ju2)::y
  real(dp),dimension(ku1:ku2)::z
  integer::ndump,nhist,nspec,mype,i,j,k

  character(LEN=80) :: filename

  ! HDF5 vars
  integer        :: ierr
  integer(hid_t) :: file_id       ! File ID

  integer, parameter :: nb_real=5,nb_int=9,nb_mpi=9

  integer , dimension(nb_int ) :: para_int
  integer , dimension(nb_mpi ) :: para_mpi
  real(dp), dimension(nb_real) :: para_real

  integer :: n1,n2,n3

  real(dp) :: d0,beta,beta_init,B0,B0_init
  namelist /restart_params/d0,beta,beta_init

  if (verbose) write (*,*) 'Entering restart_run...'

  call get_dumpname(ndump,mype,filename)

  call H5open_f(ierr)
  call H5Fopen_f(trim(filename),H5F_ACC_RDONLY_F,file_id,ierr)

  !Dataset "para_real"; Type real; dim=1; size=nb_real
  call get_1d_array_h5(file_id,"para_real",para_real,nb_real)
  time=para_real(1)
  dt  =para_real(2)
  dx  =para_real(3)
  dy  =para_real(4)
  dz  =para_real(5)

  !Dataset "para_int"; Type integer; dim=1; size=nb_int
  call get_1d_array_int_h5(file_id,"para_int",para_int,nb_int)
  ndump  =para_int(1)
  nhist  =para_int(2)
  nspec  =para_int(3)
  nx     =para_int(4)
  ny     =para_int(5)
  nz     =para_int(6)
#if WITHMPI==1
  nxslice=para_int(7)
  nyslice=para_int(8)
  nzslice=para_int(9)
#endif
  
  !Dataset "x"; Type real; dim=1, size=n1
  n1=iu2-iu1+1
  call get_1d_array_h5(file_id,"x",x,n1)

  !Dataset "y"; Type real; dim=1, size=n2
  n2=ju2-ju1+1
  call get_1d_array_h5(file_id,"y",y,n2)

  !Dataset "z"; Type real; dim=1, size=n3
  n3=ku2-ku1+1
  call get_1d_array_h5(file_id,"z",z,n3)

  !Dataset "uin"; Type real; dim=4; size=(n1,n2,n3,nvar)
  call get_4d_array_h5(file_id,"uin",uin,n1,n2,n3,nvar+3)

#if WITHMPI==1
  !Dataset "para_mpi"; Type integer; dim=1; size=nb_mpi
  call get_1d_array_int_h5(file_id,"para_mpi",para_mpi,nb_mpi)
  xleft    =para_mpi(1)
  xright   =para_mpi(2)
  yleft    =para_mpi(3)
  yright   =para_mpi(4)
  zleft    =para_mpi(5)
  zright   =para_mpi(6)
  xposition=para_mpi(7)
  yposition=para_mpi(8)
  zposition=para_mpi(9)
#endif
  
  ! Terminate access to HDF5
  call H5Fclose_f(file_id,ierr)
  call H5close_f(ierr)

  !Reset vertical field
  d0=1.d0 ; beta=1.d3 ; beta_init=1.d3
  open(unit=1,file='input' ,status='old')
  read(1,restart_params)
  close(1)

  if (beta_init/=0.d0) B0_init=2.0*sqrt(d0*ciso**2/beta_init)
  if (beta     /=0.d0) B0     =2.0*sqrt(d0*ciso**2/beta     )
  uin(:,:,:,:,8     )=uin(:,:,:,:,8     )-B0_init+B0
  uin(:,:,:,:,nvar+3)=uin(:,:,:,:,nvar+3)-B0_init+B0
  
  !Calculate initial mass in computational box
  mass0=0.d0
  do k=1,nz
     do j=1,ny
        do i=1,nx
           mass0 = mass0 + uin(1,i,j,k,1)
        end do
     end do
  end do
#if WITHMPI==1
  call sumBcast(mass0,1)
#endif
  mass0=mass0*dx*dy*dz

  return
end subroutine restart_run
#else
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine restart_run(uin,x,y,z,time,dt,dx,dy,dz,ndump,nhist,nspec,mype,cool_const,kappa)
#if WITHMPI==1
  use mpi_var
#endif
  use hydro_parameters
  use stratified
  use const
  implicit none
  real(dp)::time,dt,dx,dy,dz,cool_const, kappa, rho,rhovx,rhovy,rhovz,Bxc,Byc,Bzc,E,P,csq,Ec, Em
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin
  real(dp),dimension(iu1:iu2)::x
  real(dp),dimension(ju1:ju2)::y
  real(dp),dimension(ku1:ku2)::z
  integer::ndump,nhist,nspec,mype, i,j, k

  DOUBLE PRECISION :: amp
  character(len=10) :: type 

  logical :: IsoRestart
  integer :: a,b,c
  character(LEN=6) :: snumfile,n_mype
  character(LEN=80) :: filename,filedir

!  real(dp) :: d0,beta,beta_init,B0,B0_init
  real(dp) ::B0,beta
  namelist /init_params/d0,beta,t_cool,type,stratif,gamma,amp,cool_const,floor,kappa,c2_cool,addmass,massDiff,resist,IsoRestart

  IsoRestart=.false.

  if (verbose) write (*,*) 'Entering newstart...'

  call get_dumpname(ndump,mype,filename)

  open(unit=2,file=trim(filename),status='unknown',form='unformatted')
  read (2) time,dt,dx,dy,dz
  read (2) ndump,nhist,nspec
  read (2) nx,ny,nz
#if WITHMPI==1
  read(2) nxslice,nyslice,nzslice
#else
  read(2) a,b,c
#endif
  read (2) x,y,z
  read (2) uin

 

#if WITHMPI==1
  read(2) xleft,xright,yleft,yright,zleft,zright
  read(2) xposition,yposition,zposition
#endif

  close(2)

  !Reset vertical field
  d0=1.d0 ; beta=1.d3 
  open(unit=1,file='input' ,status='old')
  read(1,init_params)
  close(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Restart from Isothermal dump, 

#if ISO==0
  if (IsoRestart==.true.) then 
	  do k=ku1,ku2-1
	     do j=ju1,ju2-1
		do i=iu1,iu2-1

		   rho=uin(1,i,j,k,1)
		   rhovx=uin(1,i,j,k,2) ; Bxc=half*(uin(1,i,j,k,6)+uin(1,i+1,j  ,k  ,6))
		   rhovy=uin(1,i,j,k,3) ; Byc=half*(uin(1,i,j,k,7)+uin(1,i  ,j+1,k  ,7))
		   rhovz=uin(1,i,j,k,4) ; Bzc=half*(uin(1,i,j,k,8)+uin(1,i  ,j  ,k+1,8))

		   call get_E(rho,rhovx,rhovy,rhovz,E,Bxc,Byc,Bzc,gamma,1.d0)
		   uin(1,i,j,k,5)=E

		end do
	     end do
	  end do
  endif
#endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


!  B0=sqrt(2.d0*1.d0*ciso**2/beta)
!  uin(:,:,:,:,8:nvar+3:3)= B0
 ! uin(:,:,:,:,7)= 0.d0 
 ! uin(:,:,:,:,6)= 0.d0
!  if (beta_init/=0.d0) B0_init=2.0*sqrt(d0*ciso**2/beta_init)
!  if (beta     /=0.d0) B0     =2.0*sqrt(d0*ciso**2/beta     )
!  uin(:,:,:,:,7)=uin(:,:,:,:,7)-B0_init+B0

  return
end subroutine restart_run
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine get_dumpname(ndump,mype,filename)
  implicit none

  integer::ndump,mype
  character(LEN=*) :: filename

  character(LEN=6) :: snumfile,n_mype
  character(LEN=80) :: filedir

  call convtoasc(ndump,snumfile)
  call convtoasc(mype,n_mype  )
  filedir='output_'//trim(snumfile)//'/'
  filename=trim(filedir)//'slices.'//trim(n_mype)

end subroutine get_dumpname
