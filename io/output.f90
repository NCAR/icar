!>------------------------------------------------------------
!!  Model Output
!!
!!  Writes all model data to a (mostly?) CF compliant netcdf file. 
!!  Eventually this should be updated to work with generic data structures rather
!!  than specifying / hard codeing every possible variable name and attribute. 
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
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
    public :: write_domain, output_init
    
    integer, parameter :: ndims = 4     !> number of dimensions in output (x,y,z,t)
    integer, parameter :: nvars = 41    !> current number of vars = 41
    !> This will be the netCDF ID for the file and data variable.
    integer :: ncid, temp_id
    !> dimension IDs
    integer :: t_id, x_id, y_id, xu_id, yv_id, soil_id
    !> arrays to store dims for a given variable
    integer :: dimids(ndims)
    integer :: dimtwo(2)
    integer :: dimtwo_time(3)
    !> variable IDs
    integer :: lat_id,lon_id,time_id
    !> array to store ALL var ids
    integer :: varid(nvars)
    
    !> the number of time steps to save in each file
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

    real, allocatable, dimension(:,:) :: last_rain, last_snow
    logical :: surface_io_only
    ! We are writing 3D data, a (ny x nz x nx) grid or (ny x nsoil x nx) grid
    integer :: nx,ny,nz,i,nsoil
    
contains
    
    !>------------------------------------------------------------
    !! Keeps accumulated precip variables in reasonable bounds with a "bucket"
    !!
    !! Check if precip variables have exceeded kPRECIP_BUCKET_SIZE and if so, "tip" into the bucket
    !! i.e. add one to the bucket and subtract kPRECIP_BUCKET_SIZE from the precip variable
    !! and repeat while precip variable > kPRECIP_BUCKET_SIZE
    !!
    !! @param domain only: rain, crain, snow, graupel and their respective buckets
    !!
    !!------------------------------------------------------------
    subroutine tip_precip_to_buckets(domain)
        implicit none
        type(domain_type), intent(inout) :: domain
        
        ! loop variables
        integer :: nx, ny, i, j
        
        nx=size(domain%rain,1)
        ny=size(domain%rain,2)
        ! 2 to n-2 so we don't process the edges 
        do j=2,ny-2
            do i=2,nx-1
                if (domain%rain(i,j)>kPRECIP_BUCKET_SIZE) then
                    do while (domain%rain(i,j)>kPRECIP_BUCKET_SIZE)
                        domain%rain(i,j) = domain%rain(i,j)-kPRECIP_BUCKET_SIZE
                        domain%rain_bucket(i,j) = domain%rain_bucket(i,j)+1
                        last_rain(i,j)=last_rain(i,j)-kPRECIP_BUCKET_SIZE
                    end do
                endif
                if (domain%crain(i,j)>kPRECIP_BUCKET_SIZE) then
                    do while (domain%crain(i,j)>kPRECIP_BUCKET_SIZE)
                        domain%crain(i,j) = domain%crain(i,j)-kPRECIP_BUCKET_SIZE
                        domain%crain_bucket(i,j) = domain%crain_bucket(i,j)+1
                    end do
                endif
                if (domain%snow(i,j)>kPRECIP_BUCKET_SIZE) then
                    do while (domain%snow(i,j)>kPRECIP_BUCKET_SIZE)
                        domain%snow(i,j) = domain%snow(i,j)-kPRECIP_BUCKET_SIZE
                        domain%snow_bucket(i,j) = domain%snow_bucket(i,j)+1
                        last_snow(i,j)=last_snow(i,j)-kPRECIP_BUCKET_SIZE
                    end do
                endif
                if (domain%graupel(i,j)>kPRECIP_BUCKET_SIZE) then
                    do while (domain%graupel(i,j)>kPRECIP_BUCKET_SIZE)
                        domain%graupel(i,j) = domain%graupel(i,j)-kPRECIP_BUCKET_SIZE
                        domain%graupel_bucket(i,j) = domain%graupel_bucket(i,j)+1
                    end do
                endif
            end do
        end do
    end subroutine tip_precip_to_buckets
    
    !>------------------------------------------------------------
    !! Creates a new NetCDF output file
    !!
    !! Creates a new file including all variables, dimensions, and attributes
    !!
    !! @param filename  name of NetCDF output file to be created
    !! @param options   model wide options so the correct output can be created 
    !!
    !!------------------------------------------------------------
    subroutine create_file(filename,options)
        character(len=255), intent(in) :: filename
        type(options_type), intent(in) :: options
        
        ! store real (not model) timestamp
        character(len=19) :: todays_date_time
        integer,dimension(8) :: date_time
        character(len=49) :: date_format
        character(len=5) :: UTCoffset
        character(len=MAXFILELENGTH) :: err
        
        
        call date_and_time(values=date_time,zone=UTCoffset)
        date_format='(I4,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2,":",I2.2)'
        write(todays_date_time,date_format) date_time(1:3),date_time(5:7)
    
        ! create the file (clobbering any existing files!)
        err="Creating:"//trim(filename)
        call check( nf90_create(filename, NF90_CLOBBER, ncid), trim(err))
        
        err="Creating global attributes"
        call check( nf90_put_att(ncid,NF90_GLOBAL,"Conventions","CF-1.6"), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"title","Intermediate Complexity Atmospheric Research Model output"), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"institution","National Center for Atmospheric Research"), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"source","Intermediate Complexity Atmospheric Model version:"//trim(options%version)), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"history","Created:"//todays_date_time//UTCoffset), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"references", &
                    "Gutmann et al. 2016: The Intermediate Complexity Atmospheric Model (ICAR). J.Hydrometeor. doi:10.1175/JHM-D-15-0155.1, in press."), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"comment",trim(options%comment)), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"contact","Ethan Gutmann : gutmann@ucar.edu"), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"git",VERSION), trim(err))

        call check( nf90_put_att(ncid,NF90_GLOBAL,"bucket_size",kPRECIP_BUCKET_SIZE), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"dx",options%dx), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"wind_smoothing",options%smooth_wind_distance), trim(err))
        
        ! only output linear wind parameters if ICAR is using the linear wind calculations
        if (options%physics%windtype>0) then
            call check( nf90_put_att(ncid,NF90_GLOBAL,"vert_smooth",options%lt_options%vert_smooth), trim(err))
            call check( nf90_put_att(ncid,NF90_GLOBAL,"linear_contribution",options%lt_options%linear_contribution), trim(err))
            if (options%lt_options%variable_N) then
                call check( nf90_put_att(ncid,NF90_GLOBAL,"variable_N","time varying N"), trim(err))
            else
                call check( nf90_put_att(ncid,NF90_GLOBAL,"fixed_N",options%lt_options%N_squared), trim(err))
            endif
            if (options%lt_options%remove_lowres_linear) then
                if (.not.options%lt_options%variable_N) then
                    call check( nf90_put_att(ncid,NF90_GLOBAL,"rm_lin_frac_fixed_N",options%lt_options%rm_N_squared), trim(err))
                endif
                call check( nf90_put_att(ncid,NF90_GLOBAL,"remove_lin_fraction",options%lt_options%rm_linear_contribution), trim(err))
            endif
            if (options%lt_options%spatial_linear_fields) then
                call check( nf90_put_att(ncid,NF90_GLOBAL,"spatially_variable_N","spatially varying N"), trim(err))
            endif
        endif
        ! general physics options
        call check( nf90_put_att(ncid,NF90_GLOBAL,"microphysics",options%physics%microphysics), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"advection",options%physics%advection), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"boundarylayer",options%physics%boundarylayer), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"watersurface",options%physics%watersurface), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"landsurface",options%physics%landsurface), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"radiation",options%physics%radiation), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"convection",options%physics%convection), trim(err))
        call check( nf90_put_att(ncid,NF90_GLOBAL,"windtype",options%physics%windtype), trim(err))
    
        if (options%ideal) then
            call check( nf90_put_att(ncid,NF90_GLOBAL,"ideal","True"), trim(err))
        endif
        if (options%readz) then
            call check( nf90_put_att(ncid,NF90_GLOBAL,"readz","True"), trim(err))
        endif
        if (options%readdz) then
            call check( nf90_put_att(ncid,NF90_GLOBAL,"readdz","True"), trim(err))
        endif
        if (options%debug) then
            call check( nf90_put_att(ncid,NF90_GLOBAL,"debug","True"), trim(err))
        endif
        if (options%external_winds) then
            call check( nf90_put_att(ncid,NF90_GLOBAL,"external_winds","True"), trim(err))
        endif
        if (options%lt_options%remove_lowres_linear) then
            call check( nf90_put_att(ncid,NF90_GLOBAL,"remove_lowres_linear","True"), trim(err))
        endif
        if (options%mean_winds) then
            call check( nf90_put_att(ncid,NF90_GLOBAL,"mean_winds","True"), trim(err))
        endif
        if (options%mean_fields) then
            call check( nf90_put_att(ncid,NF90_GLOBAL,"mean_fields","True"), trim(err))
        endif
        if (options%restart) then
            call check( nf90_put_att(ncid,NF90_GLOBAL,"restart","True"), trim(err))
        endif
        if (options%advect_density) then
            call check( nf90_put_att(ncid,NF90_GLOBAL,"advect_density","True"), trim(err))
        endif
        if (options%lt_options%variable_N) then
            call check( nf90_put_att(ncid,NF90_GLOBAL,"variable_N","True"), trim(err))
        endif
        if (options%use_agl_height) then
            call check( nf90_put_att(ncid,NF90_GLOBAL,"use_agl_height","True"), trim(err))
        endif
    
    
    
        ! define the dimensions
        err="Creating Dimension:"
        call check( nf90_def_dim(ncid, "time", NF90_UNLIMITED, t_id), trim(err)//"time" )
        call check( nf90_def_dim(ncid, "lon", nx, x_id), trim(err)//"lon" )
        call check( nf90_def_dim(ncid, "lat", ny, y_id), trim(err)//"lat" )
        call check( nf90_def_dim(ncid, "lev", nz, temp_id), trim(err)//"lev" )
        ! standard 3D + time dimensions
        dimids(1)=x_id
        dimids(2)=y_id
        dimids(3)=temp_id !=z
        dimids(4)=t_id
        ! 2D dimensions
        dimtwo(1)=x_id
        dimtwo(2)=y_id
        ! 2D + time dimensions
        dimtwo_time(1)=x_id
        dimtwo_time(2)=y_id
        dimtwo_time(3)=t_id
        ! dimensions for the staggered grids
        call check( nf90_def_dim(ncid, "lon_u", nx+1, xu_id), trim(err)//"lon_u" )
        call check( nf90_def_dim(ncid, "lat_v", ny+1, yv_id), trim(err)//"lat_v" )
        
        ! If ICAR is using Noah, then create a dimension for soil depth too. 
        if (options%physics%landsurface==kLSM_NOAH) then
            call check( nf90_def_dim(ncid, "depth", nsoil, soil_id), trim(err)//"depth" )
        endif
    
        ! Create the variable returns varid of the data variable
        err="Defining Variable: "
        call check( nf90_def_var(ncid, "lat", NF90_REAL, dimtwo, lat_id), trim(err)//"lat" )
        call check( nf90_put_att(ncid,lat_id,"standard_name","latitude"))
        call check( nf90_put_att(ncid,lat_id,"long_name","latitude"))
        call check( nf90_put_att(ncid,lat_id,"units","degree_north"))
    
        call check( nf90_def_var(ncid, "lon", NF90_REAL, dimtwo, lon_id), trim(err)//"lon" )
        call check( nf90_put_att(ncid,lon_id,"standard_name","longitude"))
        call check( nf90_put_att(ncid,lon_id,"long_name","longitude"))
        call check( nf90_put_att(ncid,lon_id,"units","degree_east"))
    
        call check( nf90_def_var(ncid, "time", NF90_DOUBLE, t_id, time_id), trim(err)//"time" )
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
        
        ! if we are outputting 3D fields (i.e. not surface fields only)
        if (.not.surface_io_only) then
            call check( nf90_def_var(ncid, "qv", NF90_REAL, dimids, temp_id), trim(err)//"qv" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","water_vapor_mixing_ratio"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Water Wapor Mixing Ratio"))
            call check( nf90_put_att(ncid,temp_id,"units","kg kg-1"))
            varid(1)=temp_id
    
            call check( nf90_def_var(ncid, "qc", NF90_REAL, dimids, temp_id), trim(err)//"qc" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","cloud_liquid_water_mixing_ratio"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Cloud liquid water content"))
            call check( nf90_put_att(ncid,temp_id,"units","kg kg-1"))
            varid(2)=temp_id
    
            call check( nf90_def_var(ncid, "qi", NF90_REAL, dimids, temp_id), trim(err)//"qi" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","cloud_ice_mixing_ratio"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Cloud ice content"))
            call check( nf90_put_att(ncid,temp_id,"units","kg kg-1"))
            varid(3)=temp_id
    
            call check( nf90_def_var(ncid, "qr", NF90_REAL, dimids, temp_id), trim(err)//"qr" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","mass_fraction_of_rain_with_air"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Rain water content"))
            call check( nf90_put_att(ncid,temp_id,"WARNING","Could be mixing ratio, not mass fraction w/thompson scheme"))
            call check( nf90_put_att(ncid,temp_id,"units","kg kg-1"))
            varid(4)=temp_id
    
            call check( nf90_def_var(ncid, "qs", NF90_REAL, dimids, temp_id), trim(err)//"qs" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","mass_fraction_of_snow_with_air"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Snow ice content"))
            call check( nf90_put_att(ncid,temp_id,"WARNING","Could be mixing ratio, not mass fraction w/thompson scheme"))
            call check( nf90_put_att(ncid,temp_id,"units","kg kg-1"))
            varid(5)=temp_id
    
            ! these should only be output for thompson microphysics
            if (options%physics%microphysics==kMP_THOMPSON) then
                call check( nf90_def_var(ncid, "qg", NF90_REAL, dimids, temp_id), trim(err)//"qg" )
                call check( nf90_put_att(ncid,temp_id,"standard_name","mass_fraction_of_graupel_in_air"))
                call check( nf90_put_att(ncid,temp_id,"long_name","Graupel ice content"))
                call check( nf90_put_att(ncid,temp_id,"WARNING","Could be mixing ratio, not mass fraction w/thompson scheme"))
                call check( nf90_put_att(ncid,temp_id,"units","kg kg-1"))
                varid(6)=temp_id
    
                call check( nf90_def_var(ncid, "nr", NF90_REAL, dimids, temp_id), trim(err)//"nr" )
                call check( nf90_put_att(ncid,temp_id,"standard_name","number_concentration_of_rain_particles_in_air"))
                call check( nf90_put_att(ncid,temp_id,"long_name","Rain number concentration"))
                call check( nf90_put_att(ncid,temp_id,"units","cm-3"))
                varid(7)=temp_id
    
                call check( nf90_def_var(ncid, "ni", NF90_REAL, dimids, temp_id), trim(err)//"ni" )
                call check( nf90_put_att(ncid,temp_id,"standard_name","number_concentration_of_ice_crystals_in_air"))
                call check( nf90_put_att(ncid,temp_id,"long_name","Cloud ice number concentration"))
                call check( nf90_put_att(ncid,temp_id,"units","cm-3"))
                varid(8)=temp_id
            elseif (options%physics%microphysics==kMP_MORRISON) then
                call check( nf90_def_var(ncid, "qg", NF90_REAL, dimids, temp_id), trim(err)//"qg" )
                call check( nf90_put_att(ncid,temp_id,"standard_name","mass_fraction_of_graupel_in_air"))
                call check( nf90_put_att(ncid,temp_id,"long_name","Graupel ice content"))
                call check( nf90_put_att(ncid,temp_id,"WARNING","Could be mixing ratio, not mass fraction w/thompson scheme"))
                call check( nf90_put_att(ncid,temp_id,"units","kg kg-1"))
                varid(6)=temp_id
    
                call check( nf90_def_var(ncid, "nr", NF90_REAL, dimids, temp_id), trim(err)//"nr" )
                call check( nf90_put_att(ncid,temp_id,"standard_name","number_concentration_of_rain_particles_in_air"))
                call check( nf90_put_att(ncid,temp_id,"long_name","Rain number concentration"))
                call check( nf90_put_att(ncid,temp_id,"units","cm-3"))
                varid(7)=temp_id
    
                call check( nf90_def_var(ncid, "ni", NF90_REAL, dimids, temp_id), trim(err)//"ni" )
                call check( nf90_put_att(ncid,temp_id,"standard_name","number_concentration_of_ice_crystals_in_air"))
                call check( nf90_put_att(ncid,temp_id,"long_name","Cloud ice number concentration"))
                call check( nf90_put_att(ncid,temp_id,"units","cm-3"))
                varid(8)=temp_id

                call check( nf90_def_var(ncid, "ngraupel", NF90_REAL, dimids, temp_id), trim(err)//"ngraupel" )
                call check( nf90_put_att(ncid,temp_id,"standard_name","number_concentration_of_graupel_particles_in_air"))
                call check( nf90_put_att(ncid,temp_id,"long_name","Graupel number concentration"))
                call check( nf90_put_att(ncid,temp_id,"units","cm-3"))
                varid(38)=temp_id
                
                call check( nf90_def_var(ncid, "nsnow", NF90_REAL, dimids, temp_id), trim(err)//"nsnow" )
                call check( nf90_put_att(ncid,temp_id,"standard_name","number_concentration_of_snow_particles_in_air"))
                call check( nf90_put_att(ncid,temp_id,"long_name","Snow number concentration"))
                call check( nf90_put_att(ncid,temp_id,"units","cm-3"))
                varid(39)=temp_id
            elseif (options%physics%microphysics==kMP_WSM6) then
                call check( nf90_def_var(ncid, "qg", NF90_REAL, dimids, temp_id), trim(err)//"qg" )
                call check( nf90_put_att(ncid,temp_id,"standard_name","mass_fraction_of_graupel_in_air"))
                call check( nf90_put_att(ncid,temp_id,"long_name","Graupel ice content"))
                call check( nf90_put_att(ncid,temp_id,"units","kg kg-1"))
                varid(6)=temp_id
            endif
            ! need to modify dimids for staggered grids. 
            dimids(1)=xu_id
            call check( nf90_def_var(ncid, "u", NF90_REAL, dimids, temp_id), trim(err)//"u" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","grid_eastward_wind"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Grid relative eastward wind"))
            call check( nf90_put_att(ncid,temp_id,"units","m s-1"))
            varid(9)=temp_id
            dimids(1)=x_id
            
            dimids(2)=yv_id
            call check( nf90_def_var(ncid, "v", NF90_REAL, dimids, temp_id), trim(err)//"v" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","grid_northward_wind"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Grid relative northward wind"))
            call check( nf90_put_att(ncid,temp_id,"units","m s-1"))
            varid(10)=temp_id
            dimids(2)=y_id
            
            call check( nf90_def_var(ncid, "w", NF90_REAL, dimids, temp_id), trim(err)//"w" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","upward_air_velocity"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Vertical wind"))
!             call check( nf90_put_att(ncid,temp_id,"WARNING","Grid relative (i.e. add u*dz/dx) and scaled by dx/dz"))
            call check( nf90_put_att(ncid,temp_id,"units","m s-1"))
            varid(11)=temp_id
            
            call check( nf90_def_var(ncid, "p", NF90_REAL, dimids, temp_id), trim(err)//"p")
            call check( nf90_put_att(ncid,temp_id,"standard_name","air_pressure"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Pressure"))
            call check( nf90_put_att(ncid,temp_id,"units","Pa"))
            varid(12)=temp_id
            
            call check( nf90_def_var(ncid, "th", NF90_REAL, dimids, temp_id), trim(err)//"th" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","air_potential_temperature"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Potential temperature"))
            call check( nf90_put_att(ncid,temp_id,"units","K"))
            varid(13)=temp_id
            
            call check( nf90_def_var(ncid, "z", NF90_REAL, dimids(1:3), temp_id), trim(err)//"z" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","height"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Model level height (AGL)"))
            call check( nf90_put_att(ncid,temp_id,"units","m"))
            varid(20)=temp_id
    
            call check( nf90_def_var(ncid, "rho", NF90_REAL, dimids, temp_id), trim(err)//"rho" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","air_density"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Density of dry air"))
            call check( nf90_put_att(ncid,temp_id,"units","kg m-3"))
            varid(21)=temp_id
        
            if (options%physics%windtype==kWIND_LINEAR) then
                call check( nf90_def_var(ncid, "nsq", NF90_REAL, dimids, temp_id), trim(err)//"nsq" )
                call check( nf90_put_att(ncid,temp_id,"standard_name","square_of_brunt_vaisala_frequency_in_air"))
                call check( nf90_put_att(ncid,temp_id,"long_name","Brunt Vaisala Frequency (squared)"))
                call check( nf90_put_att(ncid,temp_id,"description", "Frequency is the number of oscillations of a wave per unit time."))
                call check( nf90_put_att(ncid,temp_id,"units","s-2"))
                varid(33)=temp_id
            endif

        ! else we are just storing surface fields, create 2D fields for a lot of the 3D variables, but only store the lowest model layer
        else
            call check( nf90_def_var(ncid, "qv", NF90_REAL, dimtwo_time, temp_id), trim(err)//"qv" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","water_vapor_mixing_ratio"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Water Wapor Mixing Ratio"))
            call check( nf90_put_att(ncid,temp_id,"units","kg kg-1"))
            varid(1)=temp_id
            
            dimtwo(1)=xu_id
            call check( nf90_def_var(ncid, "u", NF90_REAL, dimtwo_time, temp_id), trim(err)//"u" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","grid_eastward_wind"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Grid relative eastward wind"))
            call check( nf90_put_att(ncid,temp_id,"units","m s-1"))
            varid(9)=temp_id
            dimtwo(1)=x_id
            
            dimtwo(2)=yv_id
            call check( nf90_def_var(ncid, "v", NF90_REAL, dimtwo_time, temp_id), trim(err)//"v" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","grid_northward_wind"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Grid relative northward wind"))
            call check( nf90_put_att(ncid,temp_id,"units","m s-1"))
            varid(10)=temp_id
            dimtwo(2)=y_id
            
            call check( nf90_def_var(ncid, "p", NF90_REAL, dimtwo_time, temp_id), trim(err)//"p" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","air_pressure"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Pressure"))
            call check( nf90_put_att(ncid,temp_id,"units","Pa"))
            varid(12)=temp_id
            
            call check( nf90_def_var(ncid, "th", NF90_REAL, dimtwo_time, temp_id), trim(err)//"th" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","air_potential_temperature"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Potential temperature"))
            call check( nf90_put_att(ncid,temp_id,"units","K"))
            varid(13)=temp_id
            
            call check( nf90_def_var(ncid, "z",  NF90_REAL, dimtwo, temp_id), trim(err)//"z" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","height"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Model level height (AGL)"))
            call check( nf90_put_att(ncid,temp_id,"units","m"))
            varid(20)=temp_id
        endif

        ! surface pressure
        call check( nf90_def_var(ncid, "ps", NF90_REAL, dimtwo_time, temp_id), trim(err)//"ps")
        call check( nf90_put_att(ncid,temp_id,"standard_name","surface_air_pressure"))
        call check( nf90_put_att(ncid,temp_id,"long_name","Surface Air Pressure"))
        call check( nf90_put_att(ncid,temp_id,"units","Pa"))
        varid(40)=temp_id

    
        ! surface precip fluxes
        call check( nf90_def_var(ncid, "rain", NF90_REAL, dimtwo_time, temp_id), trim(err)//"rain" )
        call check( nf90_put_att(ncid,temp_id,"standard_name","precipitation_amount"))
        call check( nf90_put_att(ncid,temp_id,"long_name","Combined large scale and convective rain, snow and graupel (accumulated)"))
        call check( nf90_put_att(ncid,temp_id,"units","kg m-2"))
        varid(14)=temp_id

        call check( nf90_def_var(ncid, "rain_rate", NF90_REAL, dimtwo_time, temp_id), trim(err)//"rain_rate" )
        call check( nf90_put_att(ncid,temp_id,"standard_name","precipitation_amount"))
        call check( nf90_put_att(ncid,temp_id,"long_name","Time step total combined rain, snow and graupel"))
        call check( nf90_put_att(ncid,temp_id,"units","kg m-2"))
        varid(32)=temp_id
    
        call check( nf90_def_var(ncid, "snow_rate", NF90_REAL, dimtwo_time, temp_id), trim(err)//"snow_rate" )
        call check( nf90_put_att(ncid,temp_id,"standard_name","snowfall_amount"))
        call check( nf90_put_att(ncid,temp_id,"long_name","Time step total snowfall (liquid equivalent)"))
        call check( nf90_put_att(ncid,temp_id,"units","kg m-2"))
        varid(41)=temp_id

        call check( nf90_def_var(ncid, "snow", NF90_REAL, dimtwo_time, temp_id), trim(err)//"snow" )
        call check( nf90_put_att(ncid,temp_id,"standard_name","snowfall_amount"))
        call check( nf90_put_att(ncid,temp_id,"long_name","Combined large scale and convective snow (accumulated)"))
        call check( nf90_put_att(ncid,temp_id,"units","kg m-2"))
        varid(15)=temp_id
        
        if     ((options%physics%microphysics==kMP_THOMPSON)        &
            .or.(options%physics%microphysics==kMP_MORRISON)        &
            .or.(options%physics%microphysics==kMP_WSM6))           then
            call check( nf90_def_var(ncid, "graupel", NF90_REAL, dimtwo_time, temp_id), trim(err)//"graupel" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","graupel_amount"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Combined large scale and convective graupel (accumulated)"))
            call check( nf90_put_att(ncid,temp_id,"units","kg m-2"))
            varid(16)=temp_id
        endif
        
        if (options%physics%convection>0) then
            call check( nf90_def_var(ncid, "crain", NF90_REAL, dimtwo_time, temp_id), trim(err)//"crain" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","convective_rainfall_amount"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Convective rain (accumulated)"))
            call check( nf90_put_att(ncid,temp_id,"units","kg m-2"))
            varid(17)=temp_id
        endif
        
        ! surface fluxes
        ! these should only be output for radiation packages that compute them
        if (options%physics%radiation==kRA_SIMPLE) then
            call check( nf90_def_var(ncid, "clt", NF90_REAL, dimtwo_time, temp_id), trim(err)//"clt" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","cloud_area_fraction"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Fractional cloud cover"))
            call check( nf90_put_att(ncid,temp_id,"units","[0-1]"))
            varid(22)=temp_id
        endif
        if (options%physics%radiation>0) then
            call check( nf90_def_var(ncid, "rsds", NF90_REAL, dimtwo_time, temp_id), trim(err)//"rsds" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","surface_downwelling_shortwave_flux_in_air"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Shortwave downward radiation energy flux at the surface"))
            call check( nf90_put_att(ncid,temp_id,"positive","down"))
            call check( nf90_put_att(ncid,temp_id,"units","W m-2"))
            varid(18)=temp_id

            call check( nf90_def_var(ncid, "rlds", NF90_REAL, dimtwo_time, temp_id), trim(err)//"rlds" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","surface_downwelling_longwave_flux_in_air"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Longwave downward radiation energy flux at the surface"))
            call check( nf90_put_att(ncid,temp_id,"positive","down"))
            call check( nf90_put_att(ncid,temp_id,"units","W m-2"))
            varid(19)=temp_id
        endif

        call check( nf90_def_var(ncid, "ts", NF90_REAL, dimtwo_time, temp_id), trim(err)//"ts" )
        call check( nf90_put_att(ncid,temp_id,"standard_name","surface_temperature"))
        call check( nf90_put_att(ncid,temp_id,"long_name","Land surface skin temperature"))
        call check( nf90_put_att(ncid,temp_id,"units","K"))
        varid(28)=temp_id

        ! diagnosed surface fields (e.g. T2m, U10)
        if (options%physics%landsurface>kLSM_BASIC) then
            call check( nf90_def_var(ncid, "ta2m", NF90_REAL, dimtwo_time, temp_id), trim(err)//"ta2m" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","air_temperature"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Bulk air temperature at 2m"))
            call check( nf90_put_att(ncid,temp_id,"units","K"))
            varid(34)=temp_id

            call check( nf90_def_var(ncid, "hus2m", NF90_REAL, dimtwo_time, temp_id), trim(err)//"hus2m" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","specific_humidity"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Bulk air specific humidity at 2m"))
            call check( nf90_put_att(ncid,temp_id,"units","kg kg-1"))
            varid(35)=temp_id
        endif
        
        call check( nf90_def_var(ncid, "u10m", NF90_REAL, dimtwo_time, temp_id), trim(err)//"u10m" )
        call check( nf90_put_att(ncid,temp_id,"standard_name","eastward_10m_wind_speed"))
        call check( nf90_put_att(ncid,temp_id,"long_name","Eastward wind speed at 2m"))
        call check( nf90_put_att(ncid,temp_id,"units","m s-1"))
        varid(36)=temp_id
        
        call check( nf90_def_var(ncid, "v10m", NF90_REAL, dimtwo_time, temp_id), trim(err)//"v10m" )
        call check( nf90_put_att(ncid,temp_id,"standard_name","northward_10m_wind_speed"))
        call check( nf90_put_att(ncid,temp_id,"long_name","Northward wind speed at 2m"))
        call check( nf90_put_att(ncid,temp_id,"units","m s-1"))
        varid(37)=temp_id
        
        call check( nf90_def_var(ncid, "hfss", NF90_REAL, dimtwo_time, temp_id), trim(err)//"hfss" )
        call check( nf90_put_att(ncid,temp_id,"standard_name","surface_upward_sensible_heat_flux"))
        call check( nf90_put_att(ncid,temp_id,"long_name","Sensible heat flux"))
        call check( nf90_put_att(ncid,temp_id,"positive","up"))
        call check( nf90_put_att(ncid,temp_id,"units","W m-2"))
        varid(23)=temp_id

        call check( nf90_def_var(ncid, "hfls", NF90_REAL, dimtwo_time, temp_id), trim(err)//"hfls" )
        call check( nf90_put_att(ncid,temp_id,"standard_name","surface_upward_latent_heat_flux"))
        call check( nf90_put_att(ncid,temp_id,"long_name","Latent Heat Flux"))
        call check( nf90_put_att(ncid,temp_id,"positive","up"))
        call check( nf90_put_att(ncid,temp_id,"units","W m-2"))
        varid(24)=temp_id
        
        
        ! these should only be output for lsm packages that compute them
        if (options%physics%landsurface==kLSM_NOAH) then
            call check( nf90_def_var(ncid, "rlus", NF90_REAL, dimtwo_time, temp_id), trim(err)//"rlus" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","surface_upwelling_longwave_flux_in_air"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Longwave upward radiative energy flux from the surface"))
            call check( nf90_put_att(ncid,temp_id,"positive","up"))
            call check( nf90_put_att(ncid,temp_id,"units","W m-2"))
            varid(29)=temp_id
    
            call check( nf90_def_var(ncid, "hfgs", NF90_REAL, dimtwo_time, temp_id), trim(err)//"hfgs" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","upward_heat_flux_at_ground_level_in_soil"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Ground Heat Flux"))
            call check( nf90_put_att(ncid,temp_id,"positive","up"))
            call check( nf90_put_att(ncid,temp_id,"units","W m-2"))
            varid(25)=temp_id
    
            dimids(3)=soil_id
            call check( nf90_def_var(ncid, "soil_w", NF90_REAL, dimids, temp_id), trim(err)//"soil_w" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","soil_moisture_content"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Column Soil Moisture"))
            call check( nf90_put_att(ncid,temp_id,"units","kg m-2"))
            varid(26)=temp_id
    
            call check( nf90_def_var(ncid, "soil_t", NF90_REAL, dimids, temp_id), trim(err)//"soil_t" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","soil_temperature"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Column Soil Temperature"))
            call check( nf90_put_att(ncid,temp_id,"units","K"))
            varid(27)=temp_id
    
            call check( nf90_def_var(ncid, "snw", NF90_REAL, dimtwo_time, temp_id), trim(err)//"snw" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","surface_snow_amount"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Snow water equivalent"))
            call check( nf90_put_att(ncid,temp_id,"units","kg m-2"))
            varid(30)=temp_id

            call check( nf90_def_var(ncid, "canwat", NF90_REAL, dimtwo_time, temp_id), trim(err)//"canwat" )
            call check( nf90_put_att(ncid,temp_id,"standard_name","canopy_water_amount"))
            call check( nf90_put_att(ncid,temp_id,"long_name","Canopy water content"))
            call check( nf90_put_att(ncid,temp_id,"units","kg m-2"))
            varid(31)=temp_id
            
        endif
        
        ! End define mode. This tells netCDF we are done defining metadata.
        call check( nf90_enddef(ncid) )
        
    end subroutine create_file
    
    
    !>------------------------------------------------------------
    !! Set up the varids for an existing NetCDF file by searching for the variable name
    !!
    !! Querys an existing NetCDF file, and setups of the variable ids for all variables
    !! in the varids array. In some ways it would be better to have it create variables if
    !! they don't exist, but this way there is an implicit check that a restarted model 
    !! has consistent physics options. 
    !!
    !! Stops if a variable is missing from the NetCDF file
    !!
    !! @param ncid      Specifies the existing file by NetCDF ID
    !! @param options   Model options to determine which varids to look for setup
    !!
    !!------------------------------------------------------------
    subroutine setup_varids(ncid,options)
        integer, intent(in) :: ncid
        type(options_type),intent(in)::options
        integer :: temp_id
        character(len=MAXFILELENGTH) :: err
        
        err="setup_varids: Searching for variable in existing file: "
        call check( nf90_inq_varid(ncid, "lat", lat_id), trim(err)//"lat" )
        call check( nf90_inq_varid(ncid, "lon", lon_id), trim(err)//"lon" )
        call check( nf90_inq_varid(ncid, "time", time_id), trim(err)//"time" )
        call check( nf90_inq_varid(ncid, "qv", temp_id), trim(err)//"qv" )
        varid(1)=temp_id
        if (.not.surface_io_only) then
            call check( nf90_inq_varid(ncid, "qc", temp_id), trim(err)//"qc" )
            varid(2)=temp_id
            call check( nf90_inq_varid(ncid, "qi", temp_id), trim(err)//"qi" )
            varid(3)=temp_id
            call check( nf90_inq_varid(ncid, "qr", temp_id), trim(err)//"qr" )
            varid(4)=temp_id
            call check( nf90_inq_varid(ncid, "qs", temp_id), trim(err)//"qs" )
            varid(5)=temp_id
            ! these should only be output for thompson microphysics
            if (options%physics%microphysics==kMP_THOMPSON) then
                call check( nf90_inq_varid(ncid, "qg", temp_id), trim(err)//"qg" )
                varid(6)=temp_id
                call check( nf90_inq_varid(ncid, "nr", temp_id), trim(err)//"nr" )
                varid(7)=temp_id
                call check( nf90_inq_varid(ncid, "ni", temp_id), trim(err)//"ni" )
                varid(8)=temp_id
            elseif (options%physics%microphysics==kMP_MORRISON) then
                call check( nf90_inq_varid(ncid, "qg", temp_id), trim(err)//"qg" )
                varid(6)=temp_id
                call check( nf90_inq_varid(ncid, "nr", temp_id), trim(err)//"nr" )
                varid(7)=temp_id
                call check( nf90_inq_varid(ncid, "ni", temp_id), trim(err)//"ni" )
                varid(8)=temp_id
                call check( nf90_inq_varid(ncid, "ngraupel", temp_id), trim(err)//"ngraupel" )
                varid(38)=temp_id
                call check( nf90_inq_varid(ncid, "nsnow", temp_id), trim(err)//"nsnow" )
                varid(39)=temp_id
            elseif (options%physics%microphysics==kMP_WSM6) then
                call check( nf90_inq_varid(ncid, "qg", temp_id), trim(err)//"qg" )
                varid(6)=temp_id
            endif
            call check( nf90_inq_varid(ncid, "w",  temp_id), trim(err)//"w" )
            varid(11)=temp_id
            call check( nf90_inq_varid(ncid, "rho", temp_id), trim(err)//"rho" )
            varid(21)=temp_id
            if (options%physics%windtype==kWIND_LINEAR) then
                call check( nf90_inq_varid(ncid, "nsq", temp_id), trim(err)//"nsq" )
                varid(33)=temp_id
            endif
        endif
        call check( nf90_inq_varid(ncid, "u",  temp_id), trim(err)//"u" )
        varid(9)=temp_id
        call check( nf90_inq_varid(ncid, "v",  temp_id), trim(err)//"v" )
        varid(10)=temp_id
        call check( nf90_inq_varid(ncid, "p",  temp_id), trim(err)//"p" )
        varid(12)=temp_id
        call check( nf90_inq_varid(ncid, "th", temp_id), trim(err)//"th" )
        varid(13)=temp_id
        call check( nf90_inq_varid(ncid, "z",  temp_id), trim(err)//"z" )
        varid(20)=temp_id
        call check( nf90_inq_varid(ncid, "ps",temp_id), trim(err)//"ps" )
        varid(40)=temp_id
        ! surface precip fluxes
        call check( nf90_inq_varid(ncid, "rain",temp_id), trim(err)//"rain" )
        varid(14)=temp_id
        call check( nf90_inq_varid(ncid, "rain_rate",temp_id), trim(err)//"rain_rate" )
        varid(32)=temp_id
        call check( nf90_inq_varid(ncid, "snow_rate",temp_id), trim(err)//"snow_rate" )
        varid(41)=temp_id
        call check( nf90_inq_varid(ncid, "snow",temp_id), trim(err)//"snow" )
        varid(15)=temp_id
        
        if ((options%physics%microphysics==kMP_THOMPSON)    &
        .or.(options%physics%microphysics==kMP_MORRISON)    &
        .or.(options%physics%microphysics==kMP_WSM6))       then
            call check( nf90_inq_varid(ncid, "graupel", temp_id), trim(err)//"graupel" )
            varid(16)=temp_id
        endif
        if (options%physics%convection>0) then
            call check( nf90_inq_varid(ncid, "crain", temp_id), trim(err)//"crain" )
            varid(17)=temp_id
        endif
    
        ! surface fluxes
        ! these should only be output for radiation packages that compute them
        if (options%physics%radiation==kRA_SIMPLE) then
            call check( nf90_inq_varid(ncid, "clt", temp_id), trim(err)//"clt" )
            varid(22)=temp_id
        endif
        if (options%physics%radiation>0) then
            call check( nf90_inq_varid(ncid, "rsds", temp_id), trim(err)//"rsds" )
            varid(18)=temp_id
            call check( nf90_inq_varid(ncid, "rlds", temp_id), trim(err)//"rlds" )
            varid(19)=temp_id
        endif
        
        call check( nf90_inq_varid(ncid, "ts",     temp_id), trim(err)//"ts" )
        varid(28)=temp_id
        if (options%physics%landsurface>kLSM_BASIC) then
            call check( nf90_inq_varid(ncid, "ta2m", temp_id), trim(err)//"ta2m" )
            varid(34)=temp_id
            call check( nf90_inq_varid(ncid, "hus2m", temp_id), trim(err)//"hus2m" )
            varid(35)=temp_id
        endif
        call check( nf90_inq_varid(ncid, "u10m", temp_id), trim(err)//"u10m" )
        varid(36)=temp_id
        call check( nf90_inq_varid(ncid, "v10m", temp_id), trim(err)//"v10m" )
        varid(37)=temp_id
        call check( nf90_inq_varid(ncid, "hfss",   temp_id), trim(err)//"hfss" )
        varid(23)=temp_id
        call check( nf90_inq_varid(ncid, "hfls",  temp_id), trim(err)//"hfls" )
        varid(24)=temp_id
        
        ! these should only be output for lsm packages that compute them
        if (options%physics%landsurface==kLSM_NOAH) then
            call check( nf90_inq_varid(ncid, "rlus",  temp_id), trim(err)//"rlus" )
            varid(29)=temp_id
            call check( nf90_inq_varid(ncid, "hfgs",  temp_id), trim(err)//"hfgs" )
            varid(25)=temp_id
            call check( nf90_inq_varid(ncid, "soil_w", temp_id), trim(err)//"soil_w" )
            varid(26)=temp_id
            call check( nf90_inq_varid(ncid, "soil_t", temp_id), trim(err)//"soil_t" )
            varid(27)=temp_id
            call check( nf90_inq_varid(ncid, "snw",    temp_id), trim(err)//"snw" )
            varid(30)=temp_id
            call check( nf90_inq_varid(ncid, "canwat", temp_id), trim(err)//"canwat" )
            varid(31)=temp_id
        endif
        
    end subroutine setup_varids
    
    !>------------------------------------------------------------
    !!  Initialize the output module
    !!
    !!  Set up module level variables nx,ny,nz,nsoil and last_rain
    !!  Sets output frequency (daily, monthly, or every step)
    !!
    !!  @param[in]  domain  model domain datatype, fields must be initialized
    !!  @param[in]  options model options datatype, fields must be initialized
    !!
    !!------------------------------------------------------------
    subroutine output_init(domain, options)
        implicit none
        ! This is the name of the data file and variable we will read. 
        type(domain_type), intent(in)::domain
        type(options_type),intent(in)::options
        
        ! these are module level variables
        nx = size(domain%qv,1)
        nz = size(domain%qv,2)
        ny = size(domain%qv,3)
        nsoil = size(domain%soil_t,2)
        
        !! Sets output frequency (daily, monthly, or every step)
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

        
        
        if (.not.allocated(last_rain)) then
            allocate(last_rain(nx,ny))
            allocate(last_snow(nx,ny))
            if (options%restart) then
                last_rain=domain%rain
                last_snow=domain%snow
            else
                last_rain=0
                last_snow=0
            endif
        endif
    end subroutine output_init

    
    !>------------------------------------------------------------
    !!  Simple routine to write all domain data from this current time step to the output file. 
    !!
    !!  Checks if the file needs to be created or simply checked for required variables. 
    !!  Writes the actual data to the file (calls helper routines to create or set up the file as necessary). 
    !!  Note these are mostly instantaneous output fields, precip etc are accumulated fluxes. 
    !!
    !! @param domain        Full model domain to write to the file
    !! @param options       Model options to decipher what to write
    !! @param timestep      Only used to permit -1 to be specified to write a file with the restarted conditions
    !! @param inputfilename OPTIONAL: specify the output filename (otherwise options%output_file + date_time)
    !!
    !!------------------------------------------------------------
    subroutine write_domain(domain,options,timestep,inputfilename)
        implicit none
        ! This is the name of the data file and variable we will read. 
        type(domain_type),intent(inout)::domain
        type(options_type),intent(in)::options
        integer,intent(in)::timestep
        character(len=*),intent(in),optional :: inputfilename
        integer :: year, month, day, hour, minute, second
        integer :: output_shape(3), zlast(3)
        logical :: output_rain_rate
        
        output_rain_rate=.True.
        
        ! note, set this on every call so that we can output surface variables most of the time, and 3D variables once / month for restarts
        surface_io_only = options%surface_io_only
        
        output_shape = [nx,ny,nz]
        zlast = [1,3,2]
        
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
            call check( nf90_put_var(ncid, varid(20), reshape(domain%z, output_shape, order=zlast)) , trim(filename)//":Z")
        endif
        
        ! write the actual data
        call check( nf90_put_var(ncid, time_id,   domain%model_time/86400.0+50000, start_scalar ), trim(filename)//":Time" )
        
        if (.not.surface_io_only) then
            call check( nf90_put_var(ncid, varid(1),  reshape(domain%qv,    output_shape, order=zlast), start_three_D),    trim(filename)//":qv" )
            call check( nf90_put_var(ncid, varid(2),  reshape(domain%cloud, output_shape, order=zlast), start_three_D),    trim(filename)//":cloud" )
            call check( nf90_put_var(ncid, varid(3),  reshape(domain%ice,   output_shape, order=zlast), start_three_D),    trim(filename)//":ice" )
            call check( nf90_put_var(ncid, varid(4),  reshape(domain%qrain, output_shape, order=zlast), start_three_D),    trim(filename)//":qrain" )
            call check( nf90_put_var(ncid, varid(5),  reshape(domain%qsnow, output_shape, order=zlast), start_three_D),    trim(filename)//":qsnow" )
            ! these should only be output for thompson microphysics
            if (options%physics%microphysics==kMP_THOMPSON) then
                call check( nf90_put_var(ncid, varid(6),  reshape(domain%qgrau, output_shape, order=zlast), start_three_D),trim(filename)//":qgraupel" )
                call check( nf90_put_var(ncid, varid(7),  reshape(domain%nrain, output_shape, order=zlast), start_three_D),trim(filename)//":nrain" )
                call check( nf90_put_var(ncid, varid(8),  reshape(domain%nice,  output_shape, order=zlast), start_three_D),trim(filename)//":nice" )
            elseif (options%physics%microphysics==kMP_MORRISON) then
                call check( nf90_put_var(ncid, varid(6),  reshape(domain%qgrau, output_shape, order=zlast), start_three_D),trim(filename)//":qgraupel" )
                call check( nf90_put_var(ncid, varid(7),  reshape(domain%nrain, output_shape, order=zlast), start_three_D),trim(filename)//":nrain" )
                call check( nf90_put_var(ncid, varid(8),  reshape(domain%nice,  output_shape, order=zlast), start_three_D),trim(filename)//":nice" )
                call check( nf90_put_var(ncid, varid(38), reshape(domain%ngraupel, output_shape, order=zlast), start_three_D),trim(filename)//":ngraupel" )
                call check( nf90_put_var(ncid, varid(39), reshape(domain%nsnow,  output_shape, order=zlast), start_three_D),trim(filename)//":nsnow" )
            elseif (options%physics%microphysics==kMP_WSM6) then
                call check( nf90_put_var(ncid, varid(6),  reshape(domain%qgrau, output_shape, order=zlast), start_three_D),trim(filename)//":qgraupel" )
            endif
            
            ! for u add one to the x shape then revert
            output_shape(1)=output_shape(1)+1
            call check( nf90_put_var(ncid, varid(9),  reshape(domain%u,     output_shape, order=zlast), start_three_D),    trim(filename)//":u" )
            output_shape(1)=output_shape(1)-1
            
            ! for v add one to the y shape then revert
            output_shape(2)=output_shape(2)+1
            call check( nf90_put_var(ncid, varid(10), reshape(domain%v,     output_shape, order=zlast), start_three_D),    trim(filename)//":v" )
            output_shape(2)=output_shape(2)-1
            
            call check( nf90_put_var(ncid, varid(11), reshape(domain%w_real,output_shape, order=zlast), start_three_D),    trim(filename)//":w" )
            call check( nf90_put_var(ncid, varid(12), reshape(domain%p,     output_shape, order=zlast), start_three_D),    trim(filename)//":p" )
            call check( nf90_put_var(ncid, varid(13), reshape(domain%th,    output_shape, order=zlast), start_three_D),    trim(filename)//":th" )
            call check( nf90_put_var(ncid, varid(21), reshape(domain%rho,   output_shape, order=zlast), start_three_D),    trim(filename)//":rho" )
            
            if (options%physics%windtype==kWIND_LINEAR) then
                call check( nf90_put_var(ncid, varid(33), reshape(domain%nsquared, output_shape, order=zlast), start_three_D), trim(filename)//":nsquared" )
            endif
        else
            call check( nf90_put_var(ncid, varid(1),  domain%qv(:,1,:), start_two_D),  trim(filename)//":qv" )
            call check( nf90_put_var(ncid, varid(12), domain%p(:,1,:),  start_two_D),  trim(filename)//":p" )
            call check( nf90_put_var(ncid, varid(13), domain%th(:,1,:), start_two_D),  trim(filename)//":th" )
        endif
        
        ! surface pressure
        call check( nf90_put_var(ncid, varid(40), domain%psfc, start_two_D), trim(filename)//":ps" )
        
        ! Write precip variables (rain, snow, graupel, crain) adjusting for the internal precip bucket
        call tip_precip_to_buckets(domain)
        call check( nf90_put_var(ncid, varid(14), &
                                 domain%rain + domain%rain_bucket*kPRECIP_BUCKET_SIZE, &
                                 start_two_D), trim(filename)//":rain" )
        if (output_rain_rate) then ! don't bother outputing if this is the first step in a restart run
            call check( nf90_put_var(ncid, varid(32), &
                                     domain%rain-last_rain, &
                                     start_two_D), trim(filename)//":rainrate" )

            call check( nf90_put_var(ncid, varid(41), &
                                     domain%snow-last_snow, &
                                     start_two_D), trim(filename)//":snowrate" )
        endif
        call check( nf90_put_var(ncid, varid(15), &
                                 domain%snow + domain%snow_bucket*kPRECIP_BUCKET_SIZE, &
                                 start_two_D), trim(filename)//":snow" )
        
        if ((options%physics%microphysics==kMP_THOMPSON)    &
        .or.(options%physics%microphysics==kMP_MORRISON)    &
        .or.(options%physics%microphysics==kMP_WSM6))       then
            call check( nf90_put_var(ncid, varid(16), &
                                     domain%graupel + domain%graupel_bucket*kPRECIP_BUCKET_SIZE, &
                                     start_two_D), trim(filename)//":graupel" )
        endif
        if (options%physics%convection>0) then
            call check( nf90_put_var(ncid, varid(17), &
                                     domain%crain + domain%crain_bucket*kPRECIP_BUCKET_SIZE, &
                                     start_two_D), trim(filename)//":crain" )
        endif
        
        ! these should only be output for radiation packages that compute them
        if (options%physics%radiation>0) then
            call check( nf90_put_var(ncid, varid(18), domain%swdown,     start_two_D ), trim(filename)//":swdown" )
            call check( nf90_put_var(ncid, varid(19), domain%lwdown,     start_two_D ), trim(filename)//":lwdown" )
        endif
        if (options%physics%radiation==kRA_SIMPLE) then
            call check( nf90_put_var(ncid, varid(22), domain%cloudfrac,  start_two_D ), trim(filename)//":cloudfrac" )
        endif
        
        call check( nf90_put_var(ncid, varid(28), domain%skin_t,       start_two_D),  trim(filename)//":skin_t" )
        if (output_rain_rate) then ! don't bother outputing if this is the first step in a restart run
            if (options%physics%landsurface>kLSM_BASIC) then
                call check( nf90_put_var(ncid, varid(34), domain%T2m, start_two_D),  trim(filename)//":t2m" )
                call check( nf90_put_var(ncid, varid(35), domain%Q2m, start_two_D),  trim(filename)//":hus2m" )
            endif
            call check( nf90_put_var(ncid, varid(36), domain%u10, start_two_D),  trim(filename)//":u10m" )
            call check( nf90_put_var(ncid, varid(37), domain%v10, start_two_D),  trim(filename)//":v10m" )
            call check( nf90_put_var(ncid, varid(23), domain%sensible_heat,start_two_D),  trim(filename)//":sensible_heat" )
            call check( nf90_put_var(ncid, varid(24), domain%latent_heat,  start_two_D),  trim(filename)//":latent_heat" )
        endif
        ! these should only be output for lsm packages that compute them
        if (options%physics%landsurface==kLSM_NOAH) then
            output_shape(3) = nsoil
            call check( nf90_put_var(ncid, varid(26), reshape(domain%soil_vwc, output_shape, order=zlast), start_three_D),trim(filename)//":soil_vwc" )
            call check( nf90_put_var(ncid, varid(27), reshape(domain%soil_t,   output_shape, order=zlast), start_three_D),trim(filename)//":soil_t" )
            
            call check( nf90_put_var(ncid, varid(30), domain%snow_swe,     start_two_D),  trim(filename)//":snow_swe" )
            call check( nf90_put_var(ncid, varid(31), domain%canopy_water, start_two_D),  trim(filename)//":canopy_water" )
            
            if (output_rain_rate) then ! don't bother outputing if this is the first step in a restart run
                call check( nf90_put_var(ncid, varid(25), domain%ground_heat,  start_two_D),  trim(filename)//":ground_heat" )
                call check( nf90_put_var(ncid, varid(29), domain%lwup,         start_two_D),  trim(filename)//":lwup" )
            endif
        endif
        
        last_rain=domain%rain
        last_snow=domain%snow
        ! Close the file, freeing all resources.
        call check( nf90_close(ncid), trim(filename) )
    end subroutine write_domain
end module output
