! uin = (\rho, \rho vx, \rho vy, \rho vz, Etot, bx, by, bz)
subroutine dissip(uin,emfx,emfy,emfz,ngrid,mype)
  use hydro_parameters
  use variables , only : flux,nu,dx,dy,dz,dt,bval_type,cartesian,cylindrical,time,eta,kappa
  implicit none
  real(dp), intent(out), dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)          :: emfx,emfy,emfz
  real(dp), intent(out), dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3) :: uin
  integer, intent(in) :: ngrid,mype
!-----------------------------------------------------------------------
  if (verbose) write (*,*) 'Entering dissip...'
  ! Ohmic resistivity
  flux = 0.d0
  if (eta>0.d0) call resistivity(uin,flux,emfx,emfy,emfz,dx,dy,dz,dt,ngrid) 

  ! Navier-Stokes viscosity
  if (nu>0.d0) then
     if (cartesian) then
        flux=0.d0
        call viscosity(uin,flux,dx,dy,dz,dt,nu)
     endif
     if (cylindrical) call viscosity_II(uin,dx,dy,dz,dt,nu)
  endif

  ! Thermal diffusion 
#if ISO==0
  IF (kappa>0.d0) THEN
      flux=0.d0
      call thermal_diffusion(uin,flux,dx,dy,dz,dt)
  ENDIF
#endif

  ! Apply boundary conditions
  call bval(uin,ngrid,time,dt,dx,dy,dz,bval_type,mype)

  return
end subroutine dissip
