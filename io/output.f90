module output
	use netcdf
	use io_routines
	use data_structures
	implicit none
	private
	public :: write_domain
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
		real,allocatable,dimension(:,:,:)::pii,rho
		character(len=255) :: filename
		character(len=19) :: todays_date_time
		integer,dimension(8) :: date_time
		character(len=49) :: date_format
		character(len=5) :: UTCoffset
		! We are writing 3D data, a ny x nz x nx grid. 
		integer :: nx,ny,nz,i
		integer, parameter :: ndims = 3
		integer, parameter :: nvars=30
		! This will be the netCDF ID for the file and data variable.
		integer :: ncid, varid(nvars),temp_id,x_id,y_id,xu_id,yv_id,lat_id,lon_id,dimids(ndims)

		nx=size(domain%qv,1)
		nz=size(domain%qv,2)
		ny=size(domain%qv,3)
		
		! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
		! the file.
		if (timestep.eq.(-1)) then
			write(filename,"(A,A)") trim(options%output_file),".restart"
		else
			write(filename,"(A,I5.5)") trim(options%output_file),timestep
		endif
		if (options%debug) then
			write(*,*) trim(filename)
		endif
		
		call date_and_time(values=date_time,zone=UTCoffset)
		date_format='(I4,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2,":",I2.2)'
		write(todays_date_time,date_format),date_time(1:3),date_time(5:7)
		
		! create the file (clobbering any existing files!)
		call check( nf90_create(filename, NF90_CLOBBER, ncid) )
		
		call check( nf90_put_att(ncid,NF90_GLOBAL,"Conventions","CF-1.6"))
		call check( nf90_put_att(ncid,NF90_GLOBAL,"title","Intermediate Complexity Atmospheric Research Model output"))
		call check( nf90_put_att(ncid,NF90_GLOBAL,"institution","National Center for Atmospheric Research"))
		call check( nf90_put_att(ncid,NF90_GLOBAL,"source","Intermediate Complexity Atmospheric Model version:"//trim(options%version)))
		call check( nf90_put_att(ncid,NF90_GLOBAL,"history","Created:"//todays_date_time//UTCoffset))
		call check( nf90_put_att(ncid,NF90_GLOBAL,"references", &
					"Gutmann et al. 2014: The Intermediate Complexity Atmospheric Model. JHM (in prep)"))
		call check( nf90_put_att(ncid,NF90_GLOBAL,"comment",trim(options%comment)))
		call check( nf90_put_att(ncid,NF90_GLOBAL,"contact","Ethan Gutmann : gutmann@ucar.edu"))
		call check( nf90_put_att(ncid,NF90_GLOBAL,"git",VERSION))
		
		call check( nf90_put_att(ncid,NF90_GLOBAL,"dx",options%dx))
		call check( nf90_put_att(ncid,NF90_GLOBAL,"wind_smoothing",options%smooth_wind_distance))
		call check( nf90_put_att(ncid,NF90_GLOBAL,"vert_smooth",options%vert_smooth))
		call check( nf90_put_att(ncid,NF90_GLOBAL,"microphysics",options%physics%microphysics))
		call check( nf90_put_att(ncid,NF90_GLOBAL,"advection",options%physics%advection))
		call check( nf90_put_att(ncid,NF90_GLOBAL,"boundarylayer",options%physics%boundarylayer))
		call check( nf90_put_att(ncid,NF90_GLOBAL,"landsurface",options%physics%landsurface))
		call check( nf90_put_att(ncid,NF90_GLOBAL,"radiation",options%physics%radiation))
		call check( nf90_put_att(ncid,NF90_GLOBAL,"convection",options%physics%convection))
		call check( nf90_put_att(ncid,NF90_GLOBAL,"windtype",options%physics%windtype))
		
		if (options%ideal) then
			call check( nf90_put_att(ncid,NF90_GLOBAL,"ideal","True"))
		endif
		if (options%readz) then
			call check( nf90_put_att(ncid,NF90_GLOBAL,"readz","True"))
		endif
		if (options%readdz) then
			call check( nf90_put_att(ncid,NF90_GLOBAL,"readdz","True"))
		endif
		if (options%debug) then
			call check( nf90_put_att(ncid,NF90_GLOBAL,"debug","True"))
		endif
		if (options%external_winds) then
			call check( nf90_put_att(ncid,NF90_GLOBAL,"external_winds","True"))
		endif
		if (options%remove_lowres_linear) then
			call check( nf90_put_att(ncid,NF90_GLOBAL,"remove_lowres_linear","True"))
		endif
		if (options%mean_winds) then
			call check( nf90_put_att(ncid,NF90_GLOBAL,"mean_winds","True"))
		endif
		if (options%mean_fields) then
			call check( nf90_put_att(ncid,NF90_GLOBAL,"mean_fields","True"))
		endif
		if (options%restart) then
			call check( nf90_put_att(ncid,NF90_GLOBAL,"restart","True"))
		endif
		if (options%add_low_topo) then
			call check( nf90_put_att(ncid,NF90_GLOBAL,"add_low_topo","True"))
		endif
		if (options%advect_density) then
			call check( nf90_put_att(ncid,NF90_GLOBAL,"advect_density","True"))
		endif
		
		
		
		! define the dimensions
		call check( nf90_def_dim(ncid, "lon", nx, x_id) )
		dimids(1)=x_id
		call check( nf90_def_dim(ncid, "lev", nz, temp_id) )
		dimids(2)=temp_id
		call check( nf90_def_dim(ncid, "lat", ny, y_id) )
		dimids(3)=y_id
		call check( nf90_def_dim(ncid, "lon_u", nx+1, xu_id) )
		call check( nf90_def_dim(ncid, "lat_v", ny+1, yv_id) )
		
		! Create the variable returns varid of the data variable
		call check( nf90_def_var(ncid, "lat", NF90_REAL, dimids(1:3:2), lat_id) )
		call check( nf90_put_att(ncid,lat_id,"standard_name","latitude"))
		call check( nf90_put_att(ncid,lat_id,"long_name","latitude"))
		call check( nf90_put_att(ncid,lat_id,"units","degree_north"))
		
		call check( nf90_def_var(ncid, "lon", NF90_REAL, dimids(1:3:2), lon_id) )
		call check( nf90_put_att(ncid,lon_id,"standard_name","longitude"))
		call check( nf90_put_att(ncid,lon_id,"long_name","longitude"))
		call check( nf90_put_att(ncid,lon_id,"units","degree_east"))
		
		call check( nf90_def_var(ncid, "qv", NF90_REAL, dimids, temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","specific_humidity"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Specific humidity"))
		call check( nf90_put_att(ncid,temp_id,"units","kg/kg"))
		varid(1)=temp_id
		
		call check( nf90_def_var(ncid, "qc", NF90_REAL, dimids, temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","cloud_liquid_water_mixing_ratio"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Cloud liquid water content"))
		call check( nf90_put_att(ncid,temp_id,"units","kg/kg"))
		varid(2)=temp_id
		
		call check( nf90_def_var(ncid, "qi", NF90_REAL, dimids, temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","cloud_ice_mixing_ratio"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Cloud ice content"))
		call check( nf90_put_att(ncid,temp_id,"units","kg/kg"))
		varid(3)=temp_id
		
		call check( nf90_def_var(ncid, "qr", NF90_REAL, dimids, temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","mass_fraction_of_rain_with_air"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Rain water content"))
		call check( nf90_put_att(ncid,temp_id,"WARNING","Could be mixing ratio, not mass fraction w/thompson scheme"))
		call check( nf90_put_att(ncid,temp_id,"units","kg/kg"))
		varid(4)=temp_id
		
		call check( nf90_def_var(ncid, "qs", NF90_REAL, dimids, temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","mass_fraction_of_snow_with_air"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Snow ice content"))
		call check( nf90_put_att(ncid,temp_id,"WARNING","Could be mixing ratio, not mass fraction w/thompson scheme"))
		call check( nf90_put_att(ncid,temp_id,"units","kg/kg"))
		varid(5)=temp_id
		
		call check( nf90_def_var(ncid, "qg", NF90_REAL, dimids, temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","mass_fraction_of_graupel_in_air"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Graupel ice content"))
		call check( nf90_put_att(ncid,temp_id,"WARNING","Could be mixing ratio, not mass fraction w/thompson scheme"))
		call check( nf90_put_att(ncid,temp_id,"units","kg/kg"))
		varid(6)=temp_id
		
		call check( nf90_def_var(ncid, "nr", NF90_REAL, dimids, temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","number_concentration_of_rain_particles_in_air"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Rain number concentration"))
		call check( nf90_put_att(ncid,temp_id,"units","cm^-3"))
		varid(7)=temp_id
		
		call check( nf90_def_var(ncid, "ni", NF90_REAL, dimids, temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","number_concentration_of_ice_crystals_in_air"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Cloud ice number concentration"))
		call check( nf90_put_att(ncid,temp_id,"units","cm^-3"))
		varid(8)=temp_id
		
		dimids(1)=xu_id
		call check( nf90_def_var(ncid, "u",  NF90_REAL, dimids, temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","grid_eastward_wind"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Grid relative eastward wind"))
		call check( nf90_put_att(ncid,temp_id,"units","m/s"))
		varid(9)=temp_id
		
		dimids(1)=x_id
		dimids(3)=yv_id
		call check( nf90_def_var(ncid, "v",  NF90_REAL, dimids, temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","grid_northward_wind"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Grid relative northward wind"))
		call check( nf90_put_att(ncid,temp_id,"units","m/s"))
		varid(10)=temp_id
		
		dimids(3)=y_id
		call check( nf90_def_var(ncid, "w",  NF90_REAL, dimids, temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","upward_air_velocity"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Vertical wind"))
		call check( nf90_put_att(ncid,temp_id,"WARNING","Grid relative (i.e. add u*dz/dx) and scaled by dx/dz"))
		call check( nf90_put_att(ncid,temp_id,"units","m/s"))
		varid(11)=temp_id
		
		call check( nf90_def_var(ncid, "p",  NF90_REAL, dimids, temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","air_pressure"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Pressure"))
		call check( nf90_put_att(ncid,temp_id,"units","Pa"))
		varid(12)=temp_id
		
		call check( nf90_def_var(ncid, "th", NF90_REAL, dimids, temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","air_potential_temperature"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Potential temperature"))
		call check( nf90_put_att(ncid,temp_id,"units","K"))
		varid(13)=temp_id
		
		! surface precip fluxes
		call check( nf90_def_var(ncid, "rain", NF90_REAL, dimids(1:3:2), temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","precipitation_amount"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Combined large scale and convective rain, snow and graupel (accumulated)"))
		call check( nf90_put_att(ncid,temp_id,"units","kg/m^2"))
		varid(14)=temp_id
		
		call check( nf90_def_var(ncid, "snow", NF90_REAL, dimids(1:3:2), temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","snowfall_amount"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Combined large scale and convective snow (accumulated)"))
		call check( nf90_put_att(ncid,temp_id,"units","kg/m^2"))
		varid(15)=temp_id
		
		call check( nf90_def_var(ncid, "graupel", NF90_REAL, dimids(1:3:2), temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","graupel_amount"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Combined large scale and convective graupel (accumulated)"))
		call check( nf90_put_att(ncid,temp_id,"units","kg/m^2"))
		varid(16)=temp_id
		
		call check( nf90_def_var(ncid, "crain", NF90_REAL, dimids(1:3:2), temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","convective_rainfall_amount"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Convective rain (accumulated)"))
		call check( nf90_put_att(ncid,temp_id,"units","kg/m^2"))
		varid(17)=temp_id
		
		! surface fluxes
		call check( nf90_def_var(ncid, "rsds", NF90_REAL, dimids(1:3:2), temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","surface_downwelling_shortwave_flux_in_air"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Shortwave downward radiation energy flux at the surface"))
		call check( nf90_put_att(ncid,temp_id,"units","W/m^2"))
		varid(18)=temp_id

		call check( nf90_def_var(ncid, "rlds", NF90_REAL, dimids(1:3:2), temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","surface_downwelling_longwave_flux_in_air"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Longwave downward radiation energy flux at the surface"))
		call check( nf90_put_att(ncid,temp_id,"units","W/m^2"))
		varid(19)=temp_id

		call check( nf90_def_var(ncid, "z",  NF90_REAL, dimids, temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","height"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Model level height (AGL)"))
		call check( nf90_put_att(ncid,temp_id,"units","m"))
		varid(20)=temp_id
		
		call check( nf90_def_var(ncid, "rho", NF90_REAL, dimids, temp_id) )
		call check( nf90_put_att(ncid,temp_id,"standard_name","air_density"))
		call check( nf90_put_att(ncid,temp_id,"long_name","Density of dry air"))
		call check( nf90_put_att(ncid,temp_id,"units","kg/m^3"))
		varid(21)=temp_id

		
		! End define mode. This tells netCDF we are done defining metadata.
		call check( nf90_enddef(ncid) )
		
		! write the actual data
		call check( nf90_put_var(ncid, lat_id,    domain%lat) )
		call check( nf90_put_var(ncid, lon_id,    domain%lon) )
		call check( nf90_put_var(ncid, varid(1),  domain%qv) )
		call check( nf90_put_var(ncid, varid(2),  domain%cloud) )
		call check( nf90_put_var(ncid, varid(3),  domain%ice) )
		call check( nf90_put_var(ncid, varid(4),  domain%qrain) )
		call check( nf90_put_var(ncid, varid(5),  domain%qsnow) )
		call check( nf90_put_var(ncid, varid(6),  domain%qgrau) )
		call check( nf90_put_var(ncid, varid(7),  domain%nrain) )
		call check( nf90_put_var(ncid, varid(8),  domain%nice) )
		call check( nf90_put_var(ncid, varid(9),  domain%u) )
		call check( nf90_put_var(ncid, varid(10), domain%v) )
		call check( nf90_put_var(ncid, varid(11), domain%w) )
		call check( nf90_put_var(ncid, varid(12), domain%p) )
		call check( nf90_put_var(ncid, varid(13), domain%th) )
		call check( nf90_put_var(ncid, varid(14), domain%rain) )
		call check( nf90_put_var(ncid, varid(15), domain%snow) )
		call check( nf90_put_var(ncid, varid(16), domain%graupel) )
		call check( nf90_put_var(ncid, varid(17), domain%crain) )
		call check( nf90_put_var(ncid, varid(18), domain%swdown) )
		call check( nf90_put_var(ncid, varid(19), domain%lwdown) )
		call check( nf90_put_var(ncid, varid(20), domain%z) )
		call check( nf90_put_var(ncid, varid(21), domain%rho) )
	
		
		! Close the file, freeing all resources.
		call check( nf90_close(ncid) )
	end subroutine write_domain
end module output