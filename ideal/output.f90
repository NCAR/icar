module output
	use netcdf
	use io_routines
	use data_structures
	implicit none
contains
! 	simple routine to write all domain data from this current time step to the output file. 
!   note these are instantaneous fields, precip etc are accumulated fluxes. 
! 	u/v are destaggered first
! 	We could accumulated multiple time periods per file at some point, but this routine would
!   still serve as a good restart file
	subroutine write_domain(domain,options,timestep)
		implicit none
	    ! This is the name of the data file and variable we will read. 
		type(domain_type),intent(in)::domain
		type(options_type),intent(in)::options
		integer,intent(in)::timestep
		character(len=255) :: filename
		
		! We are writing 3D data, a ny x nz x nx grid. 
		integer :: nx,ny,nz,i
		integer, parameter :: ndims = 3
		integer, parameter :: nvars=17
		! This will be the netCDF ID for the file and data variable.
		integer :: ncid, varid(nvars),temp_id,dimids(ndims)

		nx=size(domain%qv,1)
		nz=size(domain%qv,2)
		ny=size(domain%qv,3)
		
		! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
		! the file.
		write(filename,"(A,I5.5)") trim(options%output_file),timestep
		if (options%debug) then
			write(*,*) trim(filename)
		endif
! 		create the file (clobbering any existing files!)
		call check( nf90_create(filename, NF90_CLOBBER, ncid) )
		
! 		define the dimensions
		call check( nf90_def_dim(ncid, "x", nx, temp_id) )
		dimids(1)=temp_id
		call check( nf90_def_dim(ncid, "z", nz, temp_id) )
		dimids(2)=temp_id
		call check( nf90_def_dim(ncid, "y", ny, temp_id) )
		dimids(3)=temp_id
		
		! Create the variable returns varid of the data variable
		call check( nf90_def_var(ncid, "qv", NF90_REAL, dimids, temp_id) )
		varid(1)=temp_id
		call check( nf90_def_var(ncid, "qc", NF90_REAL, dimids, temp_id) )
		varid(2)=temp_id
		call check( nf90_def_var(ncid, "qi", NF90_REAL, dimids, temp_id) )
		varid(3)=temp_id
		call check( nf90_def_var(ncid, "qr", NF90_REAL, dimids, temp_id) )
		varid(4)=temp_id
		call check( nf90_def_var(ncid, "qs", NF90_REAL, dimids, temp_id) )
		varid(5)=temp_id
		call check( nf90_def_var(ncid, "qg", NF90_REAL, dimids, temp_id) )
		varid(6)=temp_id
		call check( nf90_def_var(ncid, "nr", NF90_REAL, dimids, temp_id) )
		varid(7)=temp_id
		call check( nf90_def_var(ncid, "ni", NF90_REAL, dimids, temp_id) )
		varid(8)=temp_id
		call check( nf90_def_var(ncid, "u",  NF90_REAL, dimids, temp_id) )
		varid(9)=temp_id
		call check( nf90_def_var(ncid, "v",  NF90_REAL, dimids, temp_id) )
		varid(10)=temp_id
		call check( nf90_def_var(ncid, "w",  NF90_REAL, dimids, temp_id) )
		varid(11)=temp_id
		call check( nf90_def_var(ncid, "p",  NF90_REAL, dimids, temp_id) )
		varid(12)=temp_id
		call check( nf90_def_var(ncid, "th", NF90_REAL, dimids, temp_id) )
		varid(13)=temp_id
		
! 		surface precip fluxes
		call check( nf90_def_var(ncid, "rain", NF90_REAL, dimids(1:3:2), temp_id) )
		varid(14)=temp_id
		call check( nf90_def_var(ncid, "snow", NF90_REAL, dimids(1:3:2), temp_id) )
		varid(15)=temp_id
		call check( nf90_def_var(ncid, "graupel", NF90_REAL, dimids(1:3:2), temp_id) )
		varid(16)=temp_id
		
		! End define mode. This tells netCDF we are done defining metadata.
		call check( nf90_enddef(ncid) )
		
! 		write the actual data
		call check( nf90_put_var(ncid, varid(1),  domain%qv) )
		call check( nf90_put_var(ncid, varid(2),  domain%cloud) )
		call check( nf90_put_var(ncid, varid(3),  domain%ice) )
		call check( nf90_put_var(ncid, varid(4),  domain%qrain) )
		call check( nf90_put_var(ncid, varid(5),  domain%qsnow) )
		call check( nf90_put_var(ncid, varid(6),  domain%qgrau) )
		call check( nf90_put_var(ncid, varid(7),  domain%nrain) )
		call check( nf90_put_var(ncid, varid(8),  domain%nice) )
		call check( nf90_put_var(ncid, varid(9), (domain%u(1:nx,:,:)+domain%u(2:nx+1,:,:))/2) )
		call check( nf90_put_var(ncid, varid(10),(domain%v(:,:,1:ny)+domain%v(:,:,2:ny+1))/2) )
		call check( nf90_put_var(ncid, varid(11), domain%w) )
		call check( nf90_put_var(ncid, varid(12), domain%p) )
		call check( nf90_put_var(ncid, varid(13), domain%th) )
		call check( nf90_put_var(ncid, varid(14), domain%rain) )
		call check( nf90_put_var(ncid, varid(15), domain%snow) )
		call check( nf90_put_var(ncid, varid(16), domain%graupel) )
		
		! Close the file, freeing all resources.
		call check( nf90_close(ncid) )
	end subroutine write_domain
end module output