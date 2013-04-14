program ideal
	use init                ! init_*  routines
! 	use boundary_conditions ! bc_*    routines, constant in an ideal simulation? 
	use data_structures     ! *_type  types
	use microphysics        ! mp*     routines
	use wind                ! update_winds
	use advection           ! advect
	use output              ! write_* routines
	
	implicit none
	type(options_type) :: options
	type(domain_type)  :: domain
	real :: dt
	integer::i
	
	call init_model("ideal_options.namelist",options,domain,boundary)
	call mp_init(options%physics%microphysics)
	
	! linear theory winds if specified else just balance dudx+dvdy+dwdz = 0
	call update_winds(domain,options)
	
	! courant condition for 3D advection... should make 3 x 1D to maximize dt?
	dt=options%dx/maxval(domain%u)/3.0
	write(*,*) "dt = ",dt
	
	do i=0,options%ntimesteps
		call advect(domain,options,dt,options%dx)
		call mp(domain,options,dt)
! 		call lsm(domain)
! 		call pbl(domain)
! 		call radiation(domain)
		call write_domain(domain,options,i)
	end do
	
end program ideal

! # simple python code to visualize results.  
! def vis(q,n=100,vmax=None,vmin=None):
!   try:
!     for i,cq in enumerate(q[:n]):
!       plt.clf()
!       plt.imshow(cq[:,:,1].T,vmax=vmax,vmin=vmin)
!       plt.colorbar()
!       plt.title(str(i))
!       plt.draw()
!       time.sleep(0.1)
!     
!   # allow us to break-in without spitting traceback error to the screen. 
!   except KeyboardInterrupt:
!     pass
! 
! w=io.read_files("ideal*","w")
! u=io.read_files("ideal*","u")
! th=io.read_files("ideal*","th")
! qg=io.read_files("ideal*","qg")
! qv=io.read_files("ideal*","qv")
! qc=io.read_files("ideal*","qc")
! qi=io.read_files("ideal*","qi")
! qr=io.read_files("ideal*","qr")
! qs=io.read_files("ideal*","qs")
! rain=io.read_files("ideal*","rain")
! snow=io.read_files("ideal*","snow")
! graupel=io.read_files("ideal*","graupel")
! raintotal=np.concatenate([r[np.newaxis,:,:] for r in rain])
! snowtotal=np.concatenate([s[np.newaxis,:,:] for s in snow])
! graupeltotal=np.concatenate([g[np.newaxis,:,:] for g in graupel])
! plt.clf();plt.plot(raintotal.sum(axis=0)[:,1])
! plt.plot(snowtotal.sum(axis=0)[:,1])
