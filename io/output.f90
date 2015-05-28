!>------------------------------------------------------------
!!
!!  Model Output
!!
!!  Writes all model data to a (mostly?) CF compliant netcdf file
!!  Currently does not contain a time variable. 
!!  Ideally this should be added and the output file should just be updated
!!  so ICAR doesn't generate quite so many small output files
!!
!!  Author: Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module output
    use netcdf        ! nc_* routines
    use time          ! calendar_date, calendar, YEAR_ZERO, GREGORIAN, NOLEAP, THREESIXTY
    use string        ! str
    use io_routines   ! io_* routines, file_exists, check
    use data_structures
    implicit none
    private
    public :: write_domain
    
    integer, parameter :: ndims = 4
    integer, parameter :: nvars=33 ! current max = 33
    ! This will be the netCDF ID for the file and data variable.
    integer :: ncid, temp_id
    ! dimension IDs
    integer :: t_id, x_id, y_id, xu_id, yv_id, soil_id
    integer :: dimids(ndims)
    integer :: dimtwo(2)
    integer :: dimtwo_time(3)
    ! variable IDs
    integer :: lat_id,lon_id,time_id
    integer :: varid(nvars)
    
    ! the number of time steps to save in each file
    integer, parameter :: time_steps_per_file = 24
    integer, parameter :: EVERY_STEP = 0
    integer, parameter :: MONTHLY_FREQUENCY = 2
    integer, parameter :: DAILY_FREQUENCY = 1
    integer :: output_frequency = DAILY_FREQUENCY
    integer :: current_step
    character(len=255) :: filename = "default_output.nc"
    
    integer :: start_three_D(4) = [1,1,1,1]
    integer :: start_two_D(3)  = [1,1,1]
    integer :: start_scalar(1) = [1]

    real, allocatable, dimension(:,:) :: last_rain
    logical :: surface_io_only
    ! We are writing 3D data, a (ny x nz x nx) grid or (ny x nsoil x nx) grid
    integer :: nx,ny,nz,i,nsoil
    
contains
    subroutine create_file(filename,options)
        character(len=255), intent(in) :: filename
        type(options_type), intent(in) :: options
        
        character(len=19) :: todays_date_time
        integer,dimension(8) :: date_time
        character(len=49) :: date_format
        character(len=5) :: UTCoffset
        
        
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
                    "Gutmann et al. 2015: The Intermediate Complexity Atmospheric Model. JHM (in prep)"))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"comment",trim(options%comment)))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"contact","Ethan Gutmann : gutmann@ucar.edu"))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"git",VERSION))
    
        call check( nf90_put_att(ncid,NF90_GLOBAL,"dx",options%dx))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"wind_smoothing",options%smooth_wind_distance))
        
        if (options%physics%windtype>0) then
            call check( nf90_put_att(ncid,NF90_GLOBAL,"vert_smooth",options%vert_smooth))
            call check( nf90_put_att(ncid,NF90_GLOBAL,"linear_contribution",options%lt_options%linear_contribution))
            if (options%lt_options%variable_N) then
                call check( nf90_put_att(ncid,NF90_GLOBAL,"variable_N","time varying N"))
            else
                call check( nf90_put_att(ncid,NF90_GLOBAL,"fixed_N",options%lt_options%N_squared))
            endif
            if (options%lt_options%remove_lowres_linear) then
                if (.not.options%lt_options%variable_N) then
                    call check( nf90_put_att(ncid,NF90_GLOBAL,"rm_lin_frac_fixed_N",options%lt_options%rm_N_squared))
                endif
                call check( nf90_put_att(ncid,NF90_GLOBAL,"remove_lin_fraction",options%lt_options%rm_linear_contribution))
            endif
            if (options%lt_options%spatial_linear_fields) then
                call check( nf90_put_att(ncid,NF90_GLOBAL,"variable_N","spatially varying N"))
            endif
        endif
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
        if (options%lt_options%remove_lowres_linear) then
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
        if (options%advect_density) then
            call check( nf90_put_att(ncid,NF90_GLOBAL,"advect_density","True"))
        endif
        if (options%lt_options%variable_N) then
            call check( nf90_put_att(ncid,NF90_GLOBAL,"variable_N","True"))
        endif
        if (options%use_agl_height) then
            call check( nf90_put_att(ncid,NF90_GLOBAL,"use_agl_height","True"))
        endif
    
    
    
        ! define the dimensions
        call check( nf90_def_dim(ncid, "time", NF90_UNLIMITED, t_id) )
        dimids(4)=t_id
        dimtwo_time(3)=t_id
        call check( nf90_def_dim(ncid, "lon", nx, x_id) )
        dimids(1)=x_id
        dimtwo_time(1)=x_id
        dimtwo(1)=x_id
        call check( nf90_def_dim(ncid, "lev", nz, temp_id) )
        dimids(2)=temp_id
        call check( nf90_def_dim(ncid, "lat", ny, y_id) )
        dimids(3)=y_id
        dimtwo_time(2)=y_id
        dimtwo(2)=y_id
        call check( nf90_def_dim(ncid, "lon_u", nx+1, xu_id) )
        call check( nf90_def_dim(ncid, "lat_v", ny+1, yv_id) )
        
        if (options%physics%landsurface>=2) then
            call check( nf90_def_dim(ncid, "depth", nsoil, soil_id) )
        endif
    
        ! Create the variable returns varid of the data variable
        call check( nf90_def_var(ncid, "lat", NF90_REAL, dimtwo, lat_id) )
        call check( nf90_put_att(ncid,lat_id,"standard_name","latitude"))
        call check( nf90_put_att(ncid,lat_id,"long_name","latitude"))
        call check( nf90_put_att(ncid,lat_id,"units","degree_north"))
    
        call check( nf90_def_var(ncid, "lon", NF90_REAL, dimtwo, lon_id) )
        call check( nf90_put_att(ncid,lon_id,"standard_name","longitude"))
        call check( nf90_put_att(ncid,lon_id,"long_name","longitude"))
        call check( nf90_put_att(ncid,lon_id,"units","degree_east"))
    
        call check( nf90_def_var(ncid, "time", NF90_DOUBLE, t_id, time_id) )
        call check( nf90_put_att(ncid,time_id,"standard_name","time"))
        call check( nf90_put_att(ncid,time_id,"UTCoffset","0"))
        if (calendar==GREGORIAN) then
            call check( nf90_put_att(ncid,time_id,"long_name","modified Julian Day"))
            call check( nf90_put_att(ncid,time_id,"units","days since 1858-11-17 00:00:00"))
            call check( nf90_put_att(ncid,time_id,"calendar","gregorian"))
        elseif (calendar==NOLEAP) then
            call check( nf90_put_att(ncid,time_id,"long_name","Time"))
            call check( nf90_put_att(ncid,time_id,"units","days since "//trim(str(YEAR_ZERO))//"-01-01 00:00:00"))
            call check( nf90_put_att(ncid,time_id,"calendar","noleap"))
        elseif (calendar==THREESIXTY) then
            call check( nf90_put_att(ncid,time_id,"long_name","Time"))
            call check( nf90_put_att(ncid,time_id,"units","days since "//trim(str(YEAR_ZERO))//"-01-01 00:00:00"))
            call check( nf90_put_att(ncid,time_id,"calendar","360-day"))
        endif
        
        if (.not.surface_io_only) then
            call check( nf90_def_var(ncid, "qv", NF90_REAL, dimids, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","water_vapor_mixing_ratio"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Water Wapor Mixing Ratio"))
            call check( nf90_put_att(ncid,temp_id,"units","kg kg-1"))
            varid(1)=temp_id
    
            call check( nf90_def_var(ncid, "qc", NF90_REAL, dimids, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","cloud_liquid_water_mixing_ratio"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Cloud liquid water content"))
            call check( nf90_put_att(ncid,temp_id,"units","kg kg-1"))
            varid(2)=temp_id
    
            call check( nf90_def_var(ncid, "qi", NF90_REAL, dimids, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","cloud_ice_mixing_ratio"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Cloud ice content"))
            call check( nf90_put_att(ncid,temp_id,"units","kg kg-1"))
            varid(3)=temp_id
    
            call check( nf90_def_var(ncid, "qr", NF90_REAL, dimids, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","mass_fraction_of_rain_with_air"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Rain water content"))
            call check( nf90_put_att(ncid,temp_id,"WARNING","Could be mixing ratio, not mass fraction w/thompson scheme"))
            call check( nf90_put_att(ncid,temp_id,"units","kg kg-1"))
            varid(4)=temp_id
    
            call check( nf90_def_var(ncid, "qs", NF90_REAL, dimids, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","mass_fraction_of_snow_with_air"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Snow ice content"))
            call check( nf90_put_att(ncid,temp_id,"WARNING","Could be mixing ratio, not mass fraction w/thompson scheme"))
            call check( nf90_put_att(ncid,temp_id,"units","kg kg-1"))
            varid(5)=temp_id
    
            ! these should only be output for thompson microphysics
            if (options%physics%microphysics==1) then
                call check( nf90_def_var(ncid, "qg", NF90_REAL, dimids, temp_id) )
                call check( nf90_put_att(ncid,temp_id,"standard_name","mass_fraction_of_graupel_in_air"))
                call check( nf90_put_att(ncid,temp_id,"long_name","Graupel ice content"))
                call check( nf90_put_att(ncid,temp_id,"WARNING","Could be mixing ratio, not mass fraction w/thompson scheme"))
                call check( nf90_put_att(ncid,temp_id,"units","kg kg-1"))
                varid(6)=temp_id
    
                call check( nf90_def_var(ncid, "nr", NF90_REAL, dimids, temp_id) )
                call check( nf90_put_att(ncid,temp_id,"standard_name","number_concentration_of_rain_particles_in_air"))
                call check( nf90_put_att(ncid,temp_id,"long_name","Rain number concentration"))
                call check( nf90_put_att(ncid,temp_id,"units","cm-3"))
                varid(7)=temp_id
    
                call check( nf90_def_var(ncid, "ni", NF90_REAL, dimids, temp_id) )
                call check( nf90_put_att(ncid,temp_id,"standard_name","number_concentration_of_ice_crystals_in_air"))
                call check( nf90_put_att(ncid,temp_id,"long_name","Cloud ice number concentration"))
                call check( nf90_put_att(ncid,temp_id,"units","cm-3"))
                varid(8)=temp_id
            endif
            dimids(1)=xu_id
            call check( nf90_def_var(ncid, "u",  NF90_REAL, dimids, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","grid_eastward_wind"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Grid relative eastward wind"))
            call check( nf90_put_att(ncid,temp_id,"units","m s-1"))
            varid(9)=temp_id
            dimids(1)=x_id
            
            dimids(3)=yv_id
            call check( nf90_def_var(ncid, "v",  NF90_REAL, dimids, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","grid_northward_wind"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Grid relative northward wind"))
            call check( nf90_put_att(ncid,temp_id,"units","m s-1"))
            varid(10)=temp_id
            dimids(3)=y_id
            
            call check( nf90_def_var(ncid, "w",  NF90_REAL, dimids, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","upward_air_velocity"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Vertical wind"))
            call check( nf90_put_att(ncid,temp_id,"WARNING","Grid relative (i.e. add u*dz/dx) and scaled by dx/dz"))
            call check( nf90_put_att(ncid,temp_id,"units","m s-1"))
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
            
            call check( nf90_def_var(ncid, "z",  NF90_REAL, dimids(1:3), temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","height"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Model level height (AGL)"))
            call check( nf90_put_att(ncid,temp_id,"units","m"))
            varid(20)=temp_id
    
            call check( nf90_def_var(ncid, "rho", NF90_REAL, dimids, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","air_density"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Density of dry air"))
            call check( nf90_put_att(ncid,temp_id,"units","kg m-3"))
            varid(21)=temp_id
        
            if (options%physics%windtype==1) then
                call check( nf90_def_var(ncid, "nsq", NF90_REAL, dimids, temp_id) )
                call check( nf90_put_att(ncid,temp_id,"standard_name","square_of_brunt_vaisala_frequency_in_air"))
                call check( nf90_put_att(ncid,temp_id,"long_name","Brunt Vaisala Frequency (squared)"))
                call check( nf90_put_att(ncid,temp_id,"description", "Frequency is the number of oscillations of a wave per unit time."))
                call check( nf90_put_att(ncid,temp_id,"units","s-2"))
                varid(33)=temp_id
            endif

        else
            call check( nf90_def_var(ncid, "qv", NF90_REAL, dimtwo_time, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","water_vapor_mixing_ratio"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Water Wapor Mixing Ratio"))
            call check( nf90_put_att(ncid,temp_id,"units","kg kg-1"))
            varid(1)=temp_id
            
            dimtwo(1)=xu_id
            call check( nf90_def_var(ncid, "u",  NF90_REAL, dimtwo_time, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","grid_eastward_wind"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Grid relative eastward wind"))
            call check( nf90_put_att(ncid,temp_id,"units","m s-1"))
            varid(9)=temp_id
            dimtwo(1)=x_id
            
            dimtwo(2)=yv_id
            call check( nf90_def_var(ncid, "v",  NF90_REAL, dimtwo_time, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","grid_northward_wind"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Grid relative northward wind"))
            call check( nf90_put_att(ncid,temp_id,"units","m s-1"))
            varid(10)=temp_id
            dimtwo(2)=y_id
            
            call check( nf90_def_var(ncid, "p",  NF90_REAL, dimtwo_time, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","air_pressure"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Pressure"))
            call check( nf90_put_att(ncid,temp_id,"units","Pa"))
            varid(12)=temp_id
            
            call check( nf90_def_var(ncid, "th", NF90_REAL, dimtwo_time, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","air_potential_temperature"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Potential temperature"))
            call check( nf90_put_att(ncid,temp_id,"units","K"))
            varid(13)=temp_id
            
            call check( nf90_def_var(ncid, "z",  NF90_REAL, dimtwo, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","height"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Model level height (AGL)"))
            call check( nf90_put_att(ncid,temp_id,"units","m"))
            varid(20)=temp_id
        endif
    
        ! surface precip fluxes
        call check( nf90_def_var(ncid, "rain", NF90_REAL, dimtwo_time, temp_id) )
        call check( nf90_put_att(ncid,temp_id,"standard_name","precipitation_amount"))
        call check( nf90_put_att(ncid,temp_id,"long_name","Combined large scale and convective rain, snow and graupel (accumulated)"))
        call check( nf90_put_att(ncid,temp_id,"units","kg m-2"))
        varid(14)=temp_id

        call check( nf90_def_var(ncid, "rain_rate", NF90_REAL, dimtwo_time, temp_id) )
        call check( nf90_put_att(ncid,temp_id,"standard_name","precipitation_amount"))
        call check( nf90_put_att(ncid,temp_id,"long_name","Combined large scale and convective rain, snow and graupel"))
        call check( nf90_put_att(ncid,temp_id,"units","kg m-2"))
        varid(32)=temp_id
    
        call check( nf90_def_var(ncid, "snow", NF90_REAL, dimtwo_time, temp_id) )
        call check( nf90_put_att(ncid,temp_id,"standard_name","snowfall_amount"))
        call check( nf90_put_att(ncid,temp_id,"long_name","Combined large scale and convective snow (accumulated)"))
        call check( nf90_put_att(ncid,temp_id,"units","kg m-2"))
        varid(15)=temp_id
        
        if (options%physics%microphysics==1) then
            call check( nf90_def_var(ncid, "graupel", NF90_REAL, dimtwo_time, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","graupel_amount"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Combined large scale and convective graupel (accumulated)"))
            call check( nf90_put_att(ncid,temp_id,"units","kg m-2"))
            varid(16)=temp_id
        endif
        
        if (options%physics%convection>0) then
            call check( nf90_def_var(ncid, "crain", NF90_REAL, dimtwo_time, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","convective_rainfall_amount"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Convective rain (accumulated)"))
            call check( nf90_put_att(ncid,temp_id,"units","kg m-2"))
            varid(17)=temp_id
        endif
        
        ! surface fluxes
        ! these should only be output for radiation packages that compute them
        if (options%physics%radiation>=2) then
            call check( nf90_def_var(ncid, "clt", NF90_REAL, dimtwo_time, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","cloud_area_fraction"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Fractional cloud cover"))
            call check( nf90_put_att(ncid,temp_id,"units","[0-1]"))
            varid(22)=temp_id
    
            call check( nf90_def_var(ncid, "rsds", NF90_REAL, dimtwo_time, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","surface_downwelling_shortwave_flux_in_air"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Shortwave downward radiation energy flux at the surface"))
            call check( nf90_put_att(ncid,temp_id,"positive","down"))
            call check( nf90_put_att(ncid,temp_id,"units","W m-2"))
            varid(18)=temp_id

            call check( nf90_def_var(ncid, "rlds", NF90_REAL, dimtwo_time, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","surface_downwelling_longwave_flux_in_air"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Longwave downward radiation energy flux at the surface"))
            call check( nf90_put_att(ncid,temp_id,"positive","down"))
            call check( nf90_put_att(ncid,temp_id,"units","W m-2"))
            varid(19)=temp_id
        endif
        
        ! these should only be output for lsm packages that compute them
        if (options%physics%landsurface>=2) then
            call check( nf90_def_var(ncid, "rlus", NF90_REAL, dimtwo_time, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","surface_upwelling_longwave_flux_in_air"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Longwave upward radiative energy flux from the surface"))
            call check( nf90_put_att(ncid,temp_id,"positive","up"))
            call check( nf90_put_att(ncid,temp_id,"units","W m-2"))
            varid(29)=temp_id
    
            call check( nf90_def_var(ncid, "shs", NF90_REAL, dimtwo_time, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","sensible_heat_flux"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Sensible heat flux"))
            call check( nf90_put_att(ncid,temp_id,"positive","up"))
            call check( nf90_put_att(ncid,temp_id,"units","W m-2"))
            varid(23)=temp_id
    
            call check( nf90_def_var(ncid, "hfls", NF90_REAL, dimtwo_time, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","surface_upward_latent_heat_flux"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Latent Heat Flux"))
            call check( nf90_put_att(ncid,temp_id,"positive","up"))
            call check( nf90_put_att(ncid,temp_id,"units","W m-2"))
            varid(24)=temp_id
    
            call check( nf90_def_var(ncid, "hfgs", NF90_REAL, dimtwo_time, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","upward_heat_flux_at_ground_level_in_soil"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Ground Heat Flux"))
            call check( nf90_put_att(ncid,temp_id,"positive","up"))
            call check( nf90_put_att(ncid,temp_id,"units","W m-2"))
            varid(25)=temp_id
    
            dimids(2)=soil_id
            call check( nf90_def_var(ncid, "soil_w", NF90_REAL, dimids, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","soil_moisture_content"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Column Soil Moisture"))
            call check( nf90_put_att(ncid,temp_id,"units","kg m-2"))
            varid(26)=temp_id
    
            call check( nf90_def_var(ncid, "soil_t", NF90_REAL, dimids, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","soil_temperature"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Column Soil Temperature"))
            call check( nf90_put_att(ncid,temp_id,"units","K"))
            varid(27)=temp_id
    
            call check( nf90_def_var(ncid, "ts", NF90_REAL, dimtwo_time, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","surface_temperature"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Land surface skin temperature"))
            call check( nf90_put_att(ncid,temp_id,"units","K"))
            varid(28)=temp_id
    
            call check( nf90_def_var(ncid, "snw", NF90_REAL, dimtwo_time, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","surface_snow_amount"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Snow water equivalent"))
            call check( nf90_put_att(ncid,temp_id,"units","kg m-2"))
            varid(30)=temp_id

            call check( nf90_def_var(ncid, "canwat", NF90_REAL, dimtwo_time, temp_id) )
            call check( nf90_put_att(ncid,temp_id,"standard_name","canopy_water_amount"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Canopy water content"))
            call check( nf90_put_att(ncid,temp_id,"units","kg m-2"))
            varid(31)=temp_id
        endif

    
        ! End define mode. This tells netCDF we are done defining metadata.
        call check( nf90_enddef(ncid) )
        
    end subroutine create_file
    
    subroutine output_init(options)
        type(options_type), intent(in) :: options
        
        if (trim(options%output_file_frequency)=="monthly") then
            write(*,*) "Outputing a file per month"
            output_frequency=MONTHLY_FREQUENCY
        else if (trim(options%output_file_frequency)=="daily") then
            write(*,*) "Outputing a file per day"
            output_frequency=DAILY_FREQUENCY
        else
            write(*,*) "Outputing a file per time step"
            output_frequency=EVERY_STEP
        endif
        
    end subroutine output_init
    
    subroutine setup_varids(ncid,options)
        integer, intent(in) :: ncid
        type(options_type),intent(in)::options
        integer :: temp_id
        
        call check( nf90_inq_varid(ncid, "lat", lat_id) )
        call check( nf90_inq_varid(ncid, "lon", lon_id) )
        call check( nf90_inq_varid(ncid, "time", time_id) )
        call check( nf90_inq_varid(ncid, "qv", temp_id) )
        varid(1)=temp_id
        if (.not.surface_io_only) then
            call check( nf90_inq_varid(ncid, "qc", temp_id) )
            varid(2)=temp_id
            call check( nf90_inq_varid(ncid, "qi", temp_id) )
            varid(3)=temp_id
            call check( nf90_inq_varid(ncid, "qr", temp_id) )
            varid(4)=temp_id
            call check( nf90_inq_varid(ncid, "qs", temp_id) )
            varid(5)=temp_id
            ! these should only be output for thompson microphysics
            if (options%physics%microphysics==1) then
                call check( nf90_inq_varid(ncid, "qg", temp_id) )
                varid(6)=temp_id
                call check( nf90_inq_varid(ncid, "nr", temp_id) )
                varid(7)=temp_id
                call check( nf90_inq_varid(ncid, "ni", temp_id) )
                varid(8)=temp_id
            endif
            call check( nf90_inq_varid(ncid, "w",  temp_id) )
            varid(11)=temp_id
            call check( nf90_inq_varid(ncid, "rho", temp_id) )
            varid(21)=temp_id
            if (options%physics%windtype==1) then
                call check( nf90_inq_varid(ncid, "nsq", temp_id) )
                varid(33)=temp_id
            endif
        endif
        call check( nf90_inq_varid(ncid, "u",  temp_id) )
        varid(9)=temp_id
        call check( nf90_inq_varid(ncid, "v",  temp_id) )
        varid(10)=temp_id
        call check( nf90_inq_varid(ncid, "p",  temp_id) )
        varid(12)=temp_id
        call check( nf90_inq_varid(ncid, "th", temp_id) )
        varid(13)=temp_id
        call check( nf90_inq_varid(ncid, "z",  temp_id) )
        varid(20)=temp_id
        ! surface precip fluxes
        call check( nf90_inq_varid(ncid, "rain",temp_id) )
        varid(14)=temp_id
        call check( nf90_inq_varid(ncid, "rain_rate",temp_id) )
        varid(32)=temp_id
        call check( nf90_inq_varid(ncid, "snow",temp_id) )
        varid(15)=temp_id
        if (options%physics%microphysics==1) then
            call check( nf90_inq_varid(ncid, "graupel", temp_id) )
            varid(16)=temp_id
        endif
        if (options%physics%convection>0) then
            call check( nf90_inq_varid(ncid, "crain", temp_id) )
            varid(17)=temp_id
        endif
    
        ! surface fluxes
        ! these should only be output for radiation packages that compute them
        if (options%physics%radiation>=2) then
            call check( nf90_inq_varid(ncid, "clt", temp_id) )
            varid(22)=temp_id
            call check( nf90_inq_varid(ncid, "rsds", temp_id) )
            varid(18)=temp_id
            call check( nf90_inq_varid(ncid, "rlds", temp_id) )
            varid(19)=temp_id
        endif
        ! these should only be output for lsm packages that compute them
        if (options%physics%landsurface>=2) then
            call check( nf90_inq_varid(ncid, "rlus",  temp_id) )
            varid(29)=temp_id
            call check( nf90_inq_varid(ncid, "shs",   temp_id) )
            varid(23)=temp_id
            call check( nf90_inq_varid(ncid, "hfls",  temp_id) )
            varid(24)=temp_id
            call check( nf90_inq_varid(ncid, "hfgs",  temp_id) )
            varid(25)=temp_id
            dimids(2)=soil_id
            call check( nf90_inq_varid(ncid, "soil_w", temp_id) )
            varid(26)=temp_id
            call check( nf90_inq_varid(ncid, "soil_t", temp_id) )
            varid(27)=temp_id
            call check( nf90_inq_varid(ncid, "ts",     temp_id) )
            varid(28)=temp_id
            call check( nf90_inq_varid(ncid, "snw",    temp_id) )
            varid(30)=temp_id
            call check( nf90_inq_varid(ncid, "canwat", temp_id) )
            varid(31)=temp_id
        endif
        
    end subroutine setup_varids
    
!   simple routine to write all domain data from this current time step to the output file. 
!   note these are instantaneous fields, precip etc are accumulated fluxes. 
!   u/v are destaggered first
!   We could accumulated multiple time periods per file at some point, but this routine would
!   still serve as a good restart file
    subroutine write_domain(domain,options,timestep,inputfilename)
        implicit none
        ! This is the name of the data file and variable we will read. 
        type(domain_type),intent(in)::domain
        type(options_type),intent(in)::options
        integer,intent(in)::timestep
        character(len=*),intent(in),optional :: inputfilename
        integer :: year, month, day, hour, minute, second
        
        ! note, set this on every call so that we can output surface variables most of the time, and 3D variables once / month for restarts
        surface_io_only = options%surface_io_only
        
        ! these are module level variables
        nx=size(domain%qv,1)
        nz=size(domain%qv,2)
        ny=size(domain%qv,3)
        if (.not.allocated(last_rain)) then
            allocate(last_rain(nx,ny))
            last_rain=0
        endif
        nsoil=size(domain%soil_t,2)
        
        current_step=1
        if (present(inputfilename)) then
            filename=inputfilename
        else
            if (timestep.eq.(-1)) then
                write(filename,"(A,A)") trim(options%output_file),"restart.nc"
            else
                call calendar_date(domain%model_time/86400.0+50000,year, month, day, hour, minute, second)
                if (output_frequency==DAILY_FREQUENCY) then
                    write(filename,'(A,i4,"_",i2.2"_"i2.2"_"i2.2"-"i2.2".nc")') trim(options%output_file),year,month,day,0,0
                    current_step=nint(((hour*60.0+minute)*60.0+second)/options%out_dt) + 1
                elseif (output_frequency==MONTHLY_FREQUENCY) then
                    write(filename,'(A,i4,"_",i2.2"_"i2.2"_"i2.2"-"i2.2".nc")') trim(options%output_file),year,month,1,0,0
                    current_step=nint((((day*24 + hour)*60.0 + minute)*60.0 + second)/options%out_dt) + 1
                else
                    write(filename,'(A,i4,"_",i2.2"_"i2.2"_"i2.2"-"i2.2".nc")') trim(options%output_file),year,month,day,hour,minute
                    current_step=1
                endif
                
            endif
        endif
        ! this is the time position to write to in the file
        start_three_D(4) = current_step
        start_two_D(3)   = current_step
        start_scalar(1)  = current_step
        
        
        call calendar_date(domain%model_time/86400.0+50000,year, month, day, hour, minute, second)
        write(*,'(A,i4,"/",i2.2"/"i2.2" "i2.2":"i2.2)') "Output Date:",year,month,day,hour,minute
        if (file_exists(filename)) then
            ! Open the file. NF90_WRITE tells netCDF we want write/append access to
            ! the file.
            call check( nf90_open(filename,NF90_WRITE,ncid))
            call setup_varids(ncid,options)
        else
            ! otherwise, create a new file
            call create_file(filename,options)
            ! and write constant (in time) variables
            call check( nf90_put_var(ncid, lat_id,    domain%lat), trim(filename)//":Latitude" )
            call check( nf90_put_var(ncid, lon_id,    domain%lon), trim(filename)//":Longitude" )
            call check( nf90_put_var(ncid, varid(20), domain%z) , trim(filename)//":Z")
        endif
        
        ! write the actual data
        call check( nf90_put_var(ncid, time_id,   domain%model_time/86400.0+50000, start_scalar ), trim(filename)//":Time" )
        
        if (.not.surface_io_only) then
            call check( nf90_put_var(ncid, varid(1),  domain%qv,    start_three_D),    trim(filename)//":qv" )
            call check( nf90_put_var(ncid, varid(2),  domain%cloud, start_three_D),    trim(filename)//":cloud" )
            call check( nf90_put_var(ncid, varid(3),  domain%ice,   start_three_D),    trim(filename)//":ice" )
            call check( nf90_put_var(ncid, varid(4),  domain%qrain, start_three_D),    trim(filename)//":qrain" )
            call check( nf90_put_var(ncid, varid(5),  domain%qsnow, start_three_D),    trim(filename)//":qsnow" )
            ! these should only be output for thompson microphysics
            if (options%physics%microphysics==1) then
                call check( nf90_put_var(ncid, varid(6),  domain%qgrau, start_three_D),trim(filename)//":qgraupel" )
                call check( nf90_put_var(ncid, varid(7),  domain%nrain, start_three_D),trim(filename)//":nrain" )
                call check( nf90_put_var(ncid, varid(8),  domain%nice,  start_three_D),trim(filename)//":nice" )
            endif
            call check( nf90_put_var(ncid, varid(9),  domain%u,     start_three_D),    trim(filename)//":u" )
            call check( nf90_put_var(ncid, varid(10), domain%v,     start_three_D),    trim(filename)//":v" )
            call check( nf90_put_var(ncid, varid(11), domain%w,     start_three_D),    trim(filename)//":w" )
            call check( nf90_put_var(ncid, varid(12), domain%p,     start_three_D),    trim(filename)//":p" )
            call check( nf90_put_var(ncid, varid(13), domain%th,    start_three_D),    trim(filename)//":th" )
            call check( nf90_put_var(ncid, varid(21), domain%rho,   start_three_D),    trim(filename)//":rho" )
            if (options%physics%windtype==1) then
                call check( nf90_put_var(ncid, varid(33), domain%nsquared, start_three_D), trim(filename)//":nsquared" )
            endif
        else
            call check( nf90_put_var(ncid, varid(1),  domain%qv(:,1,:), start_two_D),  trim(filename)//":qv" )
!             call check( nf90_put_var(ncid, varid(9),  domain%u(:,1,:),  start_two_D),  trim(filename)//":u" )
!             call check( nf90_put_var(ncid, varid(10), domain%v(:,1,:),  start_two_D),  trim(filename)//":v" )
            call check( nf90_put_var(ncid, varid(12), domain%p(:,1,:),  start_two_D),  trim(filename)//":p" )
            call check( nf90_put_var(ncid, varid(13), domain%th(:,1,:), start_two_D),  trim(filename)//":th" )
        endif
        call check( nf90_put_var(ncid, varid(14), domain%rain,           start_two_D), trim(filename)//":rain" )
        call check( nf90_put_var(ncid, varid(32), domain%rain-last_rain, start_two_D), trim(filename)//":rainrate" )
        call check( nf90_put_var(ncid, varid(15), domain%snow,           start_two_D), trim(filename)//":snow" )
        if (options%physics%microphysics==1) then
            call check( nf90_put_var(ncid, varid(16), domain%graupel,    start_two_D), trim(filename)//":graupel" )
        endif
        if (options%physics%convection>0) then
            call check( nf90_put_var(ncid, varid(17), domain%crain,      start_two_D), trim(filename)//":crain" )
        endif
        
        ! these should only be output for radiation packages that compute them
        if (options%physics%radiation>=2) then
            call check( nf90_put_var(ncid, varid(18), domain%swdown,     start_two_D ), trim(filename)//":swdown" )
            call check( nf90_put_var(ncid, varid(19), domain%lwdown,     start_two_D ), trim(filename)//":lwdown" )
            call check( nf90_put_var(ncid, varid(22), domain%cloudfrac,  start_two_D ), trim(filename)//":cloudfrac" )
        endif
        ! these should only be output for lsm packages that compute them
        if (options%physics%landsurface>=2) then
            call check( nf90_put_var(ncid, varid(23), domain%sensible_heat,start_two_D),  trim(filename)//":sensible_heat" )
            call check( nf90_put_var(ncid, varid(24), domain%latent_heat,  start_two_D),  trim(filename)//":latent_heat" )
            call check( nf90_put_var(ncid, varid(25), domain%ground_heat,  start_two_D),  trim(filename)//":ground_heat" )
            call check( nf90_put_var(ncid, varid(26), domain%soil_vwc,     start_three_D),trim(filename)//":soil_vwc" )
            call check( nf90_put_var(ncid, varid(27), domain%soil_t,       start_three_D),trim(filename)//":soil_t" )
            call check( nf90_put_var(ncid, varid(28), domain%skin_t,       start_two_D),  trim(filename)//":skin_t" )
            call check( nf90_put_var(ncid, varid(29), domain%lwup,         start_two_D),  trim(filename)//":lwup" )
            call check( nf90_put_var(ncid, varid(30), domain%snow_swe,     start_two_D),  trim(filename)//":snow_swe" )
            call check( nf90_put_var(ncid, varid(31), domain%canopy_water, start_two_D),  trim(filename)//":canopy_water" )
        endif
        
        last_rain=domain%rain
        ! Close the file, freeing all resources.
        call check( nf90_close(ncid), trim(filename) )
    end subroutine write_domain
end module output