! ---------------------------------------------------------------
!  This is the version of ramses for the dummy....
! ---------------------------------------------------------------
program dumses
  use variables
  implicit none

  integer :: ndump,nhist,nspec
  real(dp) :: thist,tdump,tspec
  real(dp) ::cpu1,cpu2
  integer :: mype,npes
! ---------------------------------------------------------------

#if WITHMPI==1
  call init_mpi(mype,npes)
#else
  mype=0 ; npes=1
#endif

  if (mype==0) write (*,*)'Welcome in DUMSES the RAMSES FOR THE DUMMY!'

  call init(ndump,nhist,nspec,mype,npes)

  if (time==0.d0) then
     call output (uin,x,y,z,time,dt,dx,dy,dz,nhist,ndump,nspec,mype,npes)
     call history(uin,x,z,time,dt,dx,dy,dz,mype,npes)
     call special(uin,x,y,z,time,dt,dx,dy,dz,nspec,mype,npes)
  endif
  tdump=time ; thist=time ; tspec=time
  ndump=ndump+1 ; nspec=nspec+1

  if (mype==0) write (*,*) 'Start main loop...'

  cpu1=0.d0 ; cpu2=0.d0

  do

      call compute_dt(uin,dx,dy,dz,dt,courant,ngrid)

      call umuscl(uin,gravin,flux,flux_pre_tot,emfx,emfy,emfz,dx,dy,dz,dt,ngrid)
      if (bval_type(1,1)==3) call bval_shear_flux(flux,emfy,time+dt/2.,dy,ngrid,mype)

      old_uin=uin
      call update(uin,flux,flux_pre_tot,emfx,emfy,emfz,bval_type,ngrid,dt)
      call source_term(mype,npes)

      if (fargo) call fargo_update(uin,emfx,emfz,dx,dy,dz,dt,ngrid,mype)

      call bval(uin,ngrid,time,dt,dx,dy,dz,bval_type,mype)

!      nu=1
     if ((nu>0.d0).or.(eta>0.d0).or.(kappa>0.d0))      call dissip(uin,emfx,emfy,emfz,ngrid,mype)

      if (mype==0) write (*,*) time,dt      
      time = time + dt

      if (( (time-dt .le. (thist+dthist)).and.(time.gt.(thist+dthist)).and.(dthist>0.d0)  ).or. &
         (dthist==0)) then
         call history(uin,x,z,time,dt,dx,dy,dz,mype,npes)
         thist = thist + dthist
         nhist = nhist + 1
      endif

      if ( (time-dt .le. tdump+dtdump) .and. (time .gt. (tdump+dtdump)).and.(dtdump>0.d0)  ) then
         call output(uin,x,y,z,time,dt,dx,dy,dz,nhist,ndump,nspec,mype,npes)
         tdump = tdump + dtdump
         ndump = ndump + 1
      endif
     
      if ( (time-dt .le. tspec+dtspec) .and. (time .gt. (tspec+dtspec)).and.(dtspec>0.d0)  ) then
         call special(uin,x,y,z,time,dt,dx,dy,dz,nspec,mype,npes)
         tspec = tspec + dtspec
         nspec = nspec + 1
      endif
     
      if (time > tlim) exit

      call cpu_time(cpu1)
      if (mype==0) write(*,*)"CPU time required (per timestep): ",cpu1-cpu2
      cpu2=cpu1

  end do

  if (mype==0) write (*,*) 'Execution terminated...'

#if WITHMPI==1
  call finalize_mpi
#endif

end program dumses
