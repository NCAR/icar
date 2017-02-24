!> ----------------------------------------------------------------------------
!!  Model Initialization includes allocating memory for boundary and domain 
!!      data structures.  It reads all of the options from the namelist
!!      file (or files).  It also reads in Lat/Lon and Terrain data.  This module
!!      also sets up geographic (and vertical) look uptables for the forcing data
!!      Finally, there is a driver routine to initialize all model physics packages
!! 
!!   The module has been updated to allow arbitrary named variables
!!       this allows the use of e.g. ERAi, but still is not as flexible as it could be
!!
!!   The use of various python wrapper scripts in helpers/ makes it easy to add new
!!       datasets, and make them conform to the expectations of the current system. 
!!      For now there are no plans to near term plans to substantially modify this. 
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module init
    use data_structures
    use io_routines,                only : io_read, &
                                           io_write3d,io_write3di, io_write
    use geo,                        only : geo_LUT, geo_interp, geo_interp2d
    use vertical_interpolation,     only : vLUT, vinterp
    use microphysics,               only : mp_init
    use advection,                  only : adv_init
    use convection,                 only : init_convection
    use planetary_boundary_layer,   only : pbl_init
    use radiation,                  only : radiation_init
    use land_surface,               only : lsm_init
    use wind,                       only : init_winds
    use initialize_options,         only : init_options
    
    implicit none
    private
    public::init_model, init_physics
    
contains
    subroutine init_model(options,domain,boundary)
        implicit none
        type(options_type), intent(inout) :: options
        type(domain_type), intent(inout):: domain
        type(bc_type), intent(inout):: boundary

        call welcome_message()
        ! read in options file
        write(*,*) "Initializing Options"
        call init_options(options)

        ! allocate and initialize the domain
        write(*,*) ""
        write(*,*) "Initializing Domain"
        call init_domain(options,domain)
        ! allocate and initialize the boundary conditions structure (includes 3D grids too...)
        ! this might be more apropriately though of as a forcing data structure (for low res model)
        write(*,*) ""
        write(*,*) "Initializing Boundaries"
        call init_bc(options,domain,boundary)
        write(*,*) ""
        call setup_bias_correction(options,domain)
        write(*,'(/ A)') "Finished basic initialization"
        write(*,'(A /)') "---------------------------------------"
        
    end subroutine init_model
    
    subroutine init_physics(options,domain)
        implicit none
        type(options_type), intent(in) :: options
        type(domain_type), intent(inout) :: domain

        ! initialize microphysics code (e.g. compute look up tables in Thompson et al)
        call mp_init(options) !this could easily be moved to init_model...
        
        call init_convection(domain,options)
        
        call pbl_init(domain,options)
        
        call radiation_init(domain,options)
        
        call lsm_init(domain,options)
        
        call adv_init(domain,options)
        
    end subroutine init_physics
    
! ------------------------------------------------------------------------------------------
!-==== Model Domain Section ====
!
! Begining of section focused on allocating and initializing the model domain data structures
!
! ------------------------------------------------------------------------------------------

    !> ----------------------------------------------------------------------------
    !!  Conver longitudes into a common format for geographic interpolation
    !!
    !!  First convert longitudes that may be -180-180 into 0-360 range first. 
    !!  If data straddle the 0 degree longitude, convert back to -180-180 so interpolation will be smooth
    !!
    !!  @param[inout]   long_data  2D array of longitudes (either -180 - 180 or 0 - 360)
    !!
    !! ----------------------------------------------------------------------------
    subroutine convert_longitudes(long_data)
        implicit none
        real, dimension(:,:), intent(inout) :: long_data
        
        where(long_data<0) long_data = 360+long_data
        
        if ((minval(long_data)<10).and.(maxval(long_data)>350)) then
            where(long_data>180) long_data = long_data-360
        endif
        
    end subroutine convert_longitudes

!   Allow running over a sub-domain, by removing the outer N grid cells from all sides of the domain (lat,lon,terrain)
    subroutine remove_edges(domain,edgesize)
        implicit none
        type(domain_type), intent(inout) :: domain
        integer, intent(in)::edgesize
        integer::nx1,ny1,nx2,ny2,nz
        real,allocatable,dimension(:,:)::temp_data
        
        nx1 = size(domain%lat,1)
        ny1 = size(domain%lat,2)
        nx2 = nx1-(edgesize*2)
        ny2 = ny1-(edgesize*2)
        
        allocate(temp_data(nx1,ny1))
        
        temp_data = domain%lat
        deallocate(domain%lat)
        allocate(domain%lat(nx2,ny2))
        domain%lat = temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)
        
        temp_data = domain%lon
        deallocate(domain%lon)
        allocate(domain%lon(nx2,ny2))

        domain%lon = temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)

        temp_data = domain%terrain
        deallocate(domain%terrain)
        allocate(domain%terrain(nx2,ny2))
        domain%terrain = temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)

        temp_data = domain%landmask
        deallocate(domain%landmask)
        allocate(domain%landmask(nx2,ny2))
        domain%landmask = temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)

        
        deallocate(temp_data)
        
        nx1 = size(domain%u_geo%lat,1)
        ny1 = size(domain%u_geo%lat,2)
        nx2 = nx1-(edgesize*2)
        ny2 = ny1-(edgesize*2)
        
        allocate(temp_data(nx1,ny1))
        temp_data = domain%u_geo%lat
        deallocate(domain%u_geo%lat)
        allocate(domain%u_geo%lat(nx2,ny2))
        domain%u_geo%lat = temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)
        
        temp_data = domain%u_geo%lon
        deallocate(domain%u_geo%lon)
        allocate(domain%u_geo%lon(nx2,ny2))
        domain%u_geo%lon = temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)

        deallocate(temp_data)
        
        nx1 = size(domain%v_geo%lat,1)
        ny1 = size(domain%v_geo%lat,2)
        nx2 = nx1-(edgesize*2)
        ny2 = ny1-(edgesize*2)
        
        allocate(temp_data(nx1,ny1))
        temp_data = domain%v_geo%lat
        deallocate(domain%v_geo%lat)
        allocate(domain%v_geo%lat(nx2,ny2))
        domain%v_geo%lat = temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)
        
        temp_data = domain%v_geo%lon
        deallocate(domain%v_geo%lon)
        allocate(domain%v_geo%lon(nx2,ny2))
        domain%v_geo%lon = temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)
        
        deallocate(temp_data)
    end subroutine remove_edges
    
!   allocate all arrays in domain
    subroutine domain_allocation(domain,options,nx,nz,ny,nsoil)
        implicit none
        type(domain_type),  intent(inout)   :: domain
        type(options_type), intent(in)      :: options
        integer,            intent(in)      :: nx,nz,ny
        integer,            intent(in), optional   :: nsoil
        
        integer :: ns
        
        ns = 4
        if (present(nsoil)) then
            ns = nsoil
        endif
        
        ! atmosphere allocation
        allocate(domain%p(nx,nz,ny))                ! air pressure [Pa]
        domain%p = 100000
        allocate(domain%u(nx+1,nz,ny))              ! eastward wind [m/s]
        domain%u = 0
        allocate(domain%v(nx,nz,ny+1))              ! northward wind [m/s]
        domain%v = 0
        allocate(domain%w(nx,nz,ny))                ! vertical wind [grid/s]
        domain%w = 0
        allocate(domain%w_real(nx,nz,ny))           ! real vertical wind [m/s] including the U,V * dz/dx component
        domain%w_real = 0
        if (options%advect_density) then
            allocate(domain%ur(nx+1,nz,ny))         ! eastward wind * density [m/s kg/m^3]
            domain%ur = 0
            allocate(domain%vr(nx,nz,ny+1))         ! northward wind * density[m/s kg/m^3]
            domain%vr = 0
            allocate(domain%wr(nx,nz,ny))           ! vertical wind * density [grid/s kg/m^3]
            domain%wr = 0
        endif
        allocate(domain%th(nx,nz,ny))               ! potential temperature [K]
        domain%th = 280
        allocate(domain%qv(nx,nz,ny))               ! water vapor [kg/kg]
        domain%qv = 0.0002
        allocate(domain%cloud(nx,nz,ny))            ! liquid cloud water content mixing ratio [kg/kg]
        domain%cloud = 0
        allocate(domain%ice(nx,nz,ny))              ! frozen cloud water content mixing ratio [kg/kg]
        domain%ice = 0
        allocate(domain%nice(nx,nz,ny))             ! cloud ice number concentration [cm-3]
        domain%nice = 0
        allocate(domain%qrain(nx,nz,ny))            ! rain mixing ratio [kg/kg]
        domain%qrain = 0
        allocate(domain%nrain(nx,nz,ny))            ! rain drop number concentration [cm-3]
        domain%nrain = 0
        allocate(domain%qsnow(nx,nz,ny))            ! snow  mixing ratio [kg/kg]
        domain%qsnow = 0
        allocate(domain%nsnow(nx,nz,ny))            ! snow number concentration [cm-3]
        domain%nsnow = 0
        allocate(domain%qgrau(nx,nz,ny))            ! graupel mixing ratio [kg/kg]
        domain%qgrau = 0
        allocate(domain%ngraupel(nx,nz,ny))         ! graupel number concentration [cm-3]
        domain%ngraupel = 0
        allocate(domain%pii(nx,nz,ny))              ! exner function
        domain%pii = 1
        allocate(domain%rho(nx,nz,ny))              ! air density [kg/m^3]
        domain%rho = 1
        allocate(domain%cloudfrac(nx,ny))           ! cloud fraction
        domain%cloudfrac = 0

        allocate(domain%t(nx,nz,ny))                ! real air temperature [K]
        domain%t = domain%th*domain%pii
        allocate(domain%p_inter(nx,nz,ny))          ! air pressure on vertical interfaces [Pa]
        domain%p_inter = 100000
        allocate(domain%Um(nx,nz,ny))               ! eastward wind on mass grid [m/s]
        domain%Um = 0
        allocate(domain%Vm(nx,nz,ny))               ! northward wind on mass grid [m/s]
        domain%Vm = 0
        allocate(domain%mut(nx,nz,ny))              ! dry mass in each grid cell (p_inter[i] - p_inter[i+1])
        domain%mut = 0
        
        ! land-atm flux allocation
        allocate(domain%rain(nx,ny))                ! accumulated total rainfall [kg/m^2]
        domain%rain = 0
        allocate(domain%crain(nx,ny))               ! accumulated convective rainfall
        domain%crain = 0
        allocate(domain%snow(nx,ny))                ! accumulated snow fall
        domain%snow = 0
        allocate(domain%graupel(nx,ny))             ! accumulated graupel fall
        domain%graupel = 0
        allocate(domain%rain_bucket(nx,ny))         ! accumulated total rainfall [kg/m^2]
        domain%rain_bucket = 0
        allocate(domain%crain_bucket(nx,ny))        ! accumulated convective rainfall
        domain%crain_bucket = 0
        allocate(domain%snow_bucket(nx,ny))         ! accumulated snow fall
        domain%snow_bucket = 0
        allocate(domain%graupel_bucket(nx,ny))      ! accumulated graupel fall
        domain%graupel_bucket = 0
        allocate(domain%current_rain(nx,ny))        ! rain fall in current time step
        domain%current_rain = 0
        allocate(domain%current_snow(nx,ny))        ! snow fall in current time step
        domain%current_snow = 0
        allocate(domain%swdown(nx,ny))              ! shortwave down at surface
        domain%swdown = 0
        allocate(domain%lwdown(nx,ny))              ! longwave down at surface
        domain%lwdown = 0
        allocate(domain%lwup(nx,ny))                ! longwave up from surface
        domain%lwup = 0

        allocate(domain%sst(nx,ny))                 ! sea surface temperature
        domain%sst = 280

        allocate(domain%sensible_heat(nx,ny))       ! sensible heat flux from surface
        domain%sensible_heat = 0
        allocate(domain%latent_heat(nx,ny))         ! latent heat flux from surface
        domain%latent_heat = 0
        allocate(domain%ground_heat(nx,ny))         ! ground heat flux into ground
        domain%ground_heat = 0
        allocate(domain%pbl_height(nx,ny))          ! planetary boundary layer height (not always used)
        domain%pbl_height = 0
        
        ! land surface allocation
        allocate(domain%soil_t(nx,ns,ny))           ! 3D soil temperature
        domain%soil_t = 280
        allocate(domain%soil_vwc(nx,ns,ny))         ! 3D soil volumetric water content
        domain%soil_vwc = 0.25
        
        allocate(domain%soil_tdeep(nx,ny))          ! deep soil temperature
        domain%soil_tdeep = 280
        allocate(domain%soil_totalmoisture(nx,ny))  ! soil column total moisture content
        domain%soil_totalmoisture = 500 ! =2000mm * 0.25 (vwc)
        allocate(domain%skin_t(nx,ny))              ! skin temperature
        domain%skin_t = 280
        allocate(domain%snow_swe(nx,ny))            ! snow water equivalent
        domain%snow_swe = 0
        
        if (options%lsm_options%monthly_vegfrac) then
            allocate(domain%vegfrac(nx,ny,12))      ! vegetation cover fraction (%)
        else
            allocate(domain%vegfrac(nx,ny,1))       ! vegetation cover fraction (%)
        endif
        domain%vegfrac = 50 !% veg cover
        
        allocate(domain%canopy_water(nx,ny))        ! canopy water content
        domain%canopy_water = 0
        allocate(domain%soil_type(nx,ny))           ! USGS soil type
        domain%soil_type = 6 ! Loam
        allocate(domain%veg_type(nx,ny))            ! Vegetation type
        domain%veg_type = 7  ! grassland
        
        allocate(domain%u10(nx,ny))                 ! 10m height U wind
        domain%u10 = 0
        allocate(domain%v10(nx,ny))                 ! 10m height V wind
        domain%v10 = 0
        allocate(domain%t2m(nx,ny))                 ! 2m height air temperature
        domain%t2m = domain%t(:,1,:)
        allocate(domain%q2m(nx,ny))                 ! 2m height air mixing ratio
        domain%q2m = domain%qv(:,1,:)
        
        allocate(domain%znt(nx,ny))                 ! surface roughness
        domain%znt = 0.2
        allocate(domain%ustar(nx,ny))               ! surface shear stress (u*)
        domain%ustar = 0
        allocate(domain%ptop(nx,ny))                ! model top pressure
        domain%ptop = domain%p(:,nz,:)
        allocate(domain%psfc(nx,ny))                ! model surface pressure
        domain%psfc = domain%p_inter(:,1,:)

    end subroutine domain_allocation

    
! interpolate intput%z to output%z assuming that input has one less grid cell 
! in the interpolate_dim dimension
    subroutine copy_z(input,output,interpolate_dim)
        implicit none
        class(interpolable_type), intent(in) :: input
        class(interpolable_type), intent(inout) :: output
        integer,intent(in)::interpolate_dim
        
        integer::nxi,nyi,nzi,nxo,nyo,nzo
        
        ! dimensions of the input data
        nxi = size(input%z,1)
        nzi = size(input%z,2)
        nyi = size(input%z,3)
        ! dimensions of the output data
        nxo = size(output%lat,1)
        nzo = nzi
        nyo = size(output%lat,2)

        if (allocated(output%z)) then
            deallocate(output%z)
        endif
        allocate(output%z(nxo,nzo,nyo))
        if (interpolate_dim==1) then
            if ((nxo/=(nxi+1)).or.(nyo/=nyi).or.(nzo/=nzi)) then
                stop("Error copying z levels from mass grid to U grid, dimensions don't match")
            endif
            output%z(2:nxo-1,:,:) = (input%z(1:nxi-1,:,:) + input%z(2:nxi,:,:))/2
            output%z(1,:,:) = input%z(1,:,:)
            output%z(nxo,:,:) = input%z(nxi,:,:)
        else if (interpolate_dim==3) then
            if ((nyo/=(nyi+1)).or.(nxo/=nxi).or.(nzo/=nzi)) then
                stop("Error copying z levels from mass grid to V grid, dimensions don't match")
            endif
            output%z(:,:,2:nyo-1) = (input%z(:,:,1:nyi-1) + input%z(:,:,2:nyi))/2
            output%z(:,:,1) = input%z(:,:,1)
            output%z(:,:,nyo) = input%z(:,:,nyi)
        else
            write(*,*) "Can not interpolate z data over z dimension"
        endif
        
    end subroutine copy_z

    subroutine setup_bias_correction(options, domain)
        implicit none
        type(options_type), intent(in) :: options
        type(domain_type), intent(inout):: domain
        
        if (options%use_bias_correction) then
            call io_read(options%bias_options%filename, options%bias_options%rain_fraction_var, domain%rain_fraction)
        endif
    end subroutine setup_bias_correction

    subroutine init_domain_land(domain,options)
        implicit none
        type(options_type), intent(in) :: options
        type(domain_type), intent(inout):: domain
        ! these are temporary variables used for IO before storing data in the domain data structure
        integer,dimension(:,:),allocatable :: temp_idata
        real,dimension(:,:,:),allocatable :: temp_rdata
        real,dimension(:,:),allocatable :: temp_rdata_2d
        
        integer:: buffer,nx,ny,nz, i
        
        buffer = options%buffer
        nx = size(domain%veg_type,1)+buffer*2
        ny = size(domain%veg_type,2)+buffer*2
        nz = size(domain%soil_t,2)
        
        
        ! Veg cover fraction = 2D/3D real
        if (options%vegfrac_var.ne."") then
            ! if we are supposed to read a veg fraction for each month then it will be a 3D variable
            if (options%lsm_options%monthly_vegfrac) then
                call io_read(options%init_conditions_file,options%vegfrac_var,temp_rdata,1)
                domain%vegfrac = temp_rdata(1+buffer:nx-buffer,1+buffer:ny-buffer,:)  ! subset the data by buffer grid cells
                deallocate(temp_rdata)
            else
                ! otherwise we just read a 2D variable and store it in the first time entry in the 3D vegfrac field
                call io_read(options%init_conditions_file,options%vegfrac_var,temp_rdata_2d,1)
                domain%vegfrac(:,:,1) = temp_rdata_2d(1+buffer:nx-buffer,1+buffer:ny-buffer)  ! subset the data by buffer grid cells
                deallocate(temp_rdata_2d)
            endif
        endif
        ! if the input vegetation fraction is a fraction (0-1) 
        ! then convert it to a percent
        if (maxval(domain%vegfrac)<=1) then
            domain%vegfrac = domain%vegfrac*100
        endif

        ! Veg TYPE = 2D integer
        if (options%vegtype_var.ne."") then
            call io_read(options%init_conditions_file,options%vegtype_var,temp_idata,1)
            domain%veg_type = temp_idata(1+buffer:nx-buffer,1+buffer:ny-buffer)  ! subset the data by buffer grid cells
            deallocate(temp_idata)
        endif
        
        ! Soil TYPE = 2D integer
        if (options%soiltype_var.ne."") then
            call io_read(options%init_conditions_file,options%soiltype_var,temp_idata,1)
            domain%soil_type = temp_idata(1+buffer:nx-buffer,1+buffer:ny-buffer)  ! subset the data by buffer grid cells
            deallocate(temp_idata)
        endif
        
        ! Soil Volumetric Water Content = 3D real
        if (options%soil_vwc_var.ne."") then
            call io_read(options%init_conditions_file,options%soil_vwc_var,temp_rdata,1)  
            domain%soil_vwc = reshape(temp_rdata(1+buffer:nx-buffer,1+buffer:ny-buffer,1:nz) &    ! subset the data by buffer grid cells
                                    ,[nx-buffer*2,nz,ny-buffer*2],order=[1,3,2])                ! and reshape to move the z axis to the middle
            deallocate(temp_rdata)
        endif
        
        ! Soil Temperature = 3D real
        if (options%soil_t_var.ne."") then
            call io_read(options%init_conditions_file,options%soil_t_var,temp_rdata,1)
            domain%soil_t = reshape(temp_rdata(1+buffer:nx-buffer,1+buffer:ny-buffer,1:nz) &  ! subset the data by buffer grid cells
                                    ,[nx-buffer*2,nz,ny-buffer*2],order=[1,3,2])            ! and reshape to move the z axis to the middle
            deallocate(temp_rdata)
        endif
        
        ! Deep Soil Temperature = 2D real
        if (options%soil_deept_var.ne."") then
            call io_read(options%init_conditions_file,options%soil_deept_var,temp_rdata_2d,1)
            domain%soil_tdeep = temp_rdata_2d(1+buffer:nx-buffer,1+buffer:ny-buffer)  ! subset the data by buffer grid cells
            where(domain%soil_tdeep<200) domain%soil_tdeep = 200 ! mitigates zeros that can cause problems
            if (options%soil_t_var=="") then
                write(*,*) "Missing explicit soil T, using deep soil T for all depths"
                do i = 1,nz
                    domain%soil_t(:,i,:) = domain%soil_tdeep
                end do
            endif
            deallocate(temp_rdata_2d)
        endif
        
        
    end subroutine init_domain_land
    
    !> ----------------------------------------------------------------------------
    !!  Interpolate a 2D data array in the "x" (first) dimension, extrapolate out one on edges
    !!  
    !!  If data_out is not allocated, it will be allocated as (nx+1,ny). Then a simple average between
    !!  grid cells computes a bilinear interpolation. Edge cells are extrapolated outwards using a
    !!  similar method. 
    !!  
    !!  @param[in]  data_in     real 2D array to be interpolated
    !!  @param[out] data_out    real 2D allocatable array to store the interpolated output
    !!
    !> ----------------------------------------------------------------------------
    subroutine interpolate_in_x(data_in, data_out)
        implicit none
        real, dimension(:,:), intent(in) :: data_in
        real, dimension(:,:), allocatable, intent(out) :: data_out
        
        integer :: nx, ny
        
        nx = size(data_in,1)
        ny = size(data_in,2)
        if (.not.allocated(data_out)) then
            allocate(data_out(nx+1,ny))
        else
            if ((size(data_out,1)/=nx+1).or.(size(data_out,2)/=ny)) then
                deallocate(data_out)
                allocate(data_out(nx+1,ny))
            endif
        endif
        
        ! simple bilinear inteprolation (included extrapolation one grid cell out)
        data_out(2:nx,:) = (data_in(2:,:) + data_in(1:nx-1,:))/2
        data_out(1,:)    = data_in(1,:)*2 - data_in(2,:)
        data_out(nx+1,:) = data_in(nx,:)*2 - data_in(nx-1,:)
        ! result is data_out
    end subroutine interpolate_in_x

    !> ----------------------------------------------------------------------------
    !!  Interpolate a 2D data array in the "y" (second) dimension, extrapolate out one on edges
    !!  
    !!  If data_out is not allocated, it will be allocated as (nx,ny+1). Then a simple average between
    !!  grid cells computes a bilinear interpolation. Edge cells are extrapolated outwards using a
    !!  similar method. 
    !!  
    !!  @param[in]  data_in     real 2D array to be interpolated
    !!  @param[out] data_out    real 2D allocatable array to store the interpolated output
    !!
    !> ----------------------------------------------------------------------------
    subroutine interpolate_in_y(data_in, data_out)
        implicit none
        real, dimension(:,:), intent(in) :: data_in
        real, dimension(:,:), allocatable, intent(out) :: data_out
        
        integer :: nx, ny
        
        nx = size(data_in,1)
        ny = size(data_in,2)
        if (.not.allocated(data_out)) then
            allocate(data_out(nx,ny+1))
        else
            if ((size(data_out,1)/=nx).or.(size(data_out,2)/=ny+1)) then
                deallocate(data_out)
                allocate(data_out(nx,ny+1))
            endif
        endif
        
        ! simple bilinear inteprolation (included extrapolation one grid cell out)
        data_out(:,2:ny) = (data_in(:,2:) + data_in(:,1:ny-1))/2
        data_out(:,1)    = data_in(:,1)*2 - data_in(:,2)
        data_out(:,ny+1) = data_in(:,ny)*2 - data_in(:,ny-1)
        ! result is data_out
    end subroutine interpolate_in_y

!   initialize the domain e.g. lat,lon,terrain, 3D z coordinate
    subroutine init_domain(options, domain)
        implicit none
        type(options_type), intent(in) :: options
        type(domain_type), intent(inout):: domain
        real,dimension(:,:,:),allocatable::temporary_z
        integer:: ny,nz,nx,i,buf
        
!       these are the only required variables on a high-res grid, lat, lon, and terrain elevation
        call io_read(options%init_conditions_file,options%hgt_hi,domain%terrain,1)
        call io_read(options%init_conditions_file,options%lat_hi,domain%lat,1)
        call io_read(options%init_conditions_file,options%lon_hi,domain%lon,1)
        call convert_longitudes(domain%lon)
        
        !  because u/vlat/lon_hi are optional, we calculated them by interplating from the mass lat/lon if necessary
        if (options%ulat_hi/="") then
            call io_read(options%init_conditions_file,options%ulat_hi,domain%u_geo%lat,1)
        else
            call interpolate_in_x(domain%lat,domain%u_geo%lat)
        endif
        if (options%ulon_hi/="") then
            call io_read(options%init_conditions_file,options%ulon_hi,domain%u_geo%lon,1)
        else
            call interpolate_in_x(domain%lon,domain%u_geo%lon)
        endif
        call convert_longitudes(domain%u_geo%lon)
        
        if (options%vlat_hi/="") then
            call io_read(options%init_conditions_file,options%vlat_hi,domain%v_geo%lat,1)
        else
            call interpolate_in_y(domain%lat,domain%v_geo%lat)
        endif
        if (options%vlon_hi/="") then
            call io_read(options%init_conditions_file,options%vlon_hi,domain%v_geo%lon,1)
        else
            call interpolate_in_y(domain%lon,domain%v_geo%lon)
        endif
        call convert_longitudes(domain%v_geo%lon)
        
        if (options%landvar/="") then
            call io_read(options%init_conditions_file,options%landvar,domain%landmask,1)
            where(domain%landmask==0) domain%landmask = kLC_WATER
        else
            nx = size(domain%lat,1)
            ny = size(domain%lat,2)
            allocate(domain%landmask(nx,ny))
            domain%landmask = kLC_LAND !if we weren't supplied a landmask field, assume all is land (what we care about anyway)
        endif
        
        ! remove n grid cells from all sides of the domain if requested
        if(options%buffer>0) then
            write(*,*) "Removing ",options%buffer," gridcells from edges of domain"
            call remove_edges(domain,options%buffer)
        endif
        
        ! use the lat variable to define the x and y dimensions for all other variables
        nx = size(domain%lat,1)
        ny = size(domain%lat,2)
        ! assumes nz is defined in the options
        nz = options%nz
        
        ! if a 3d grid was also specified, then read those data in
        if ((options%readz).and.(options%ideal).and.(options%zvar.ne."")) then
            
            stop("Reading Z from an external file is not currently supported, use a fixed dz")
            
            call io_read(options%init_conditions_file,options%zvar, domain%z)
            ! dz also has to be calculated from the 3d z file
            buf = options%buffer
            allocate(domain%dz(nx,nz,ny))
            allocate(temporary_z(nx,nz,ny))
            do i = 1,nz-1
                domain%dz(:,i,:) = domain%z(buf+1:nx+buf,buf+1:ny+buf,i+1)-domain%z(buf+1:nx+buf,buf+1:ny+buf,i)
                temporary_z(:,i,:) = domain%z(buf+1:nx+buf,buf+1:ny+buf,i)
            enddo
            temporary_z(:,nz,:) = domain%z(buf+1:nx+buf,buf+1:ny+buf,nz)
            domain%dz(:,nz,:) = domain%dz(:,nz-1,:)
            deallocate(domain%z)
            allocate(domain%z(nx,nz,ny))
            domain%z = temporary_z
            deallocate(temporary_z)
            
        else
            ! otherwise, set up the z grid to be evenly spaced in z using the terrain +dz/2 for the base
            ! and z[i-1]+(dz[i-1]+dz[i])/2 for the rest
            allocate(domain%z(nx,nz,ny))
            allocate(domain%z_inter(nx,nz,ny))
            allocate(domain%dz(nx,nz,ny))
            allocate(domain%dz_inter(nx,nz,ny))
            ! lowest model level is half of the lowest dz above the land surface
            domain%z(:,1,:) = domain%terrain+options%dz_levels(1)/2
            domain%z_inter(:,1,:) = domain%terrain
            ! dz between the mass grid points is half of each dz
            domain%dz(:,1,:) = (options%dz_levels(1) + options%dz_levels(2))/2
            ! dz between the interface points = dz = thickness of mass grid cells
            domain%dz_inter(:,1,:) = options%dz_levels(1)
            do i = 2,nz
                domain%z(:,i,:)       = domain%z(:,i-1,:)       + (options%dz_levels(i)+options%dz_levels(i-1))/2
                domain%z_inter(:,i,:) = domain%z_inter(:,i-1,:) +  options%dz_levels(i)
                
                ! dz between the interface points = dz = thickness of mass grid cells
                domain%dz_inter(:,i,:) = options%dz_levels(i)
                ! dz between the mass grid points is half of each dz
                if (i<nz) then
                    domain%dz(:,i,:) = (options%dz_levels(i) + options%dz_levels(i+1))/2
                endif
            enddo
            domain%dz(:,nz,:) = options%dz_levels(nz)
            
        endif
        call copy_z(domain,domain%u_geo,interpolate_dim=1)
        call copy_z(domain,domain%v_geo,interpolate_dim=3)
        
        ! all other variables should be allocated and initialized to 0
        call domain_allocation(domain,options,nx,nz,ny)
        ! initializing land
        call init_domain_land(domain,options)
        
        ! store dx in domain as well as options, read as an option, but it is more appropriate in domain
        domain%dx = options%dx
        
        call init_winds(domain,options)
    end subroutine init_domain

! ------------------------------------------------------------------------------------------
!-==== Boundary Conditions Section ====
!
! Begining of section focused on allocating and initializeing the boundary data structures
!
! ------------------------------------------------------------------------------------------
!   allocate arrays in boundary condition data structure
    subroutine boundary_allocate(boundary,nx,nz,ny)
        implicit none
        type(bc_type), intent(inout) :: boundary
        integer,intent(in)::nx,nz,ny
        
        allocate(boundary%du_dt(nx+1,nz,ny))
        boundary%du_dt = 0
        allocate(boundary%dv_dt(nx,nz,ny+1))
        boundary%dv_dt = 0
        allocate(boundary%dp_dt(nx,nz,ny))
        boundary%dp_dt = 0
        allocate(boundary%dth_dt(nz,max(nx,ny),4))
        boundary%dth_dt = 0
        allocate(boundary%dqv_dt(nz,max(nx,ny),4))
        boundary%dqv_dt = 0
        allocate(boundary%dqc_dt(nz,max(nx,ny),4))
        boundary%dqc_dt = 0
        allocate(boundary%dlh_dt(nx,ny))
        boundary%dlh_dt = 0
        allocate(boundary%dsh_dt(nx,ny))
        boundary%dsh_dt = 0
        allocate(boundary%dlw_dt(nx,ny))
        boundary%dlw_dt = 0
        allocate(boundary%dsw_dt(nx,ny))
        boundary%dsw_dt = 0
        allocate(boundary%dsst_dt(nx,ny))
        boundary%dsst_dt = 0
        allocate(boundary%dpblh_dt(nx,ny))
        boundary%dpblh_dt = 0
        allocate(boundary%drain_dt(nx,ny))
        boundary%drain_dt = 0
    end subroutine boundary_allocate
    
! initialize the boundary condition data structure e.g. lat,lon,terrain,3D Z coord
    subroutine init_bc_data(options,boundary,domain)
        implicit none
        type(options_type), intent(in) :: options
        type(bc_type), intent(inout):: boundary
        type(domain_type), intent(in):: domain
        real,dimension(45)::fulldz
        real,dimension(:,:,:),allocatable::zbase
        integer::nx,ny,nz,i
        
        ! these variables are required for any boundary/forcing file type
        call io_read(options%boundary_files(1),options%latvar,boundary%lat)
        call io_read(options%boundary_files(1),options%lonvar,boundary%lon)
        call convert_longitudes(boundary%lon)
        call io_read(options%boundary_files(1),options%ulat,boundary%u_geo%lat)
        call io_read(options%boundary_files(1),options%ulon,boundary%u_geo%lon)
        call convert_longitudes(boundary%u_geo%lon)
        call io_read(options%boundary_files(1),options%vlat,boundary%v_geo%lat)
        call io_read(options%boundary_files(1),options%vlon,boundary%v_geo%lon)
        call convert_longitudes(boundary%v_geo%lon)
        call io_read(options%boundary_files(1),options%hgtvar,boundary%terrain)
        
        ! read in the vertical coordinate
        call io_read(options%boundary_files(1),options%zvar,zbase)
        if (options%debug) write(*,*) "Raw input forcing z min=",minval(zbase), " z max=", maxval(zbase)
        nx = size(zbase,1)
        ny = size(zbase,2)
        nz = size(zbase,3)
        allocate(boundary%lowres_z(nx,nz,ny))
        boundary%lowres_z = reshape(zbase,[nx,nz,ny],order=[1,3,2])
        deallocate(zbase)

        if (trim(options%zbvar)/="") then
            call io_read(options%boundary_files(1),options%zbvar, zbase)
            boundary%lowres_z = boundary%lowres_z + reshape(zbase,[nx,nz,ny],order=[1,3,2])
            if (options%debug) write(*,*) "Raw forcing z_base min=",minval(zbase), " max=", maxval(zbase)
            deallocate(zbase)
        endif
        if (options%z_is_geopotential) then
             boundary%lowres_z = boundary%lowres_z/ gravity
        endif
        if (options%debug) write(*,*) "Pre-interpolation forcing z min=",minval(boundary%lowres_z), " z max=", maxval(boundary%lowres_z)
        
        if (options%z_is_on_interface) then
            write(*,*) "Interpreting geopotential height as residing between model layers"
            nz = nz - 1
            ! zbase is used as a temporary variable while we deallocate/reallocate bc%z
            if (allocated(zbase)) deallocate(zbase)
            allocate(zbase(nx,nz,ny))
            zbase = (boundary%lowres_z(:,1:nz,:) + boundary%lowres_z(:,2:nz+1,:))/2
            deallocate(boundary%lowres_z)
            allocate(boundary%lowres_z(nx,nz,ny))
            boundary%lowres_z = zbase
            deallocate(zbase)
        endif
        
        ! all other structures must be allocated and initialized, but will be set on a high-res grid
        ! u/v are seperate so we can read them on the low res grid and adjust/rm-linearwinds before interpolating
        ! this also makes it easier to change how these variables are read from various forcing model file structures
        nx = size(boundary%u_geo%lon,1)
        ny = size(boundary%u_geo%lon,2)
        allocate(boundary%u(nx,nz,ny))
        boundary%u = 0
        nx = size(boundary%v_geo%lon,1)
        ny = size(boundary%v_geo%lon,2)
        allocate(boundary%v(nx,nz,ny))
        boundary%v = 0
        
        nz = options%nz
        nx = size(domain%lat,1)
        ny = size(domain%lat,2)
        
        call boundary_allocate(boundary,nx,nz,ny)
    end subroutine init_bc_data
    
    
!   sets up the data used to rotate the wind field when using an external
!   high-res wind field from e.g. a high res WRF run
    subroutine setup_extwinds(domain)
        implicit none
        type(wind_type),intent(inout)::domain
        integer::nx,ny
        
        nx = size(domain%terrain,1)
        ny = size(domain%terrain,2)
        
        ! dzdx/y used in rotating windfield back to terrain following grid in a simple fashion
        allocate(domain%dzdx(nx-1,ny))
        allocate(domain%dzdy(nx,ny-1))
        domain%dzdx = sqrt((domain%terrain(2:nx,:)-domain%terrain(1:nx-1,:))**2+domain%dx**2)/domain%dx
        domain%dzdy = sqrt((domain%terrain(:,2:ny)-domain%terrain(:,1:ny-1))**2+domain%dx**2)/domain%dx
    end subroutine setup_extwinds
    
    
!   initialize the external wind system (primarily GEOLUTs)
    subroutine init_ext_winds(options,bc)
        implicit none
        type(options_type), intent(in) :: options
        type(bc_type),intent(inout) :: bc
            
        real, allocatable, dimension(:,:,:) :: u,v
        real, allocatable, dimension(:,:) :: lat,lon
        real, allocatable, dimension(:,:) :: temporary_terrain
        
        call io_read(options%ext_wind_files(1),options%hgt_hi,bc%ext_winds%terrain,1)
        call io_read(options%ext_wind_files(1),options%latvar,bc%ext_winds%lat)
        call io_read(options%ext_wind_files(1),options%lonvar,bc%ext_winds%lon)
        call convert_longitudes(bc%ext_winds%lon)
        
        ! ulat_hi and vlat_hi are optional, if they are not supplied interpolate the mass lat/lon grid
        if (options%ulat_hi/="") then
            call io_read(options%ext_wind_files(1),options%ulat_hi,bc%ext_winds%u_geo%lat,1)
        else
            call interpolate_in_x(bc%ext_winds%lat,bc%ext_winds%u_geo%lat)
        endif
        if (options%ulon_hi/="") then
            call io_read(options%ext_wind_files(1),options%ulon_hi,bc%ext_winds%u_geo%lon,1)
        else
            call interpolate_in_x(bc%ext_winds%lon,bc%ext_winds%u_geo%lon)
        endif
        call convert_longitudes(bc%ext_winds%u_geo%lon)
        
        if (options%vlat_hi/="") then
            call io_read(options%ext_wind_files(1),options%vlat_hi,bc%ext_winds%v_geo%lat,1)
        else
            call interpolate_in_y(bc%ext_winds%lat,bc%ext_winds%v_geo%lat)
        endif
        if (options%vlon_hi/="") then
            call io_read(options%ext_wind_files(1),options%vlon_hi,bc%ext_winds%v_geo%lon,1)
        else
            call interpolate_in_y(bc%ext_winds%lon,bc%ext_winds%v_geo%lon)
        endif
        call convert_longitudes(bc%ext_winds%v_geo%lon)

        
        write(*,*) "Setting up ext wind geoLUTs"
        call geo_LUT(bc%next_domain%u_geo, bc%ext_winds%u_geo)
        call geo_LUT(bc%next_domain%v_geo, bc%ext_winds%v_geo)
        
        ! if the external wind file has a different shape in either dimension, compute the GEOLUT and interpolate terrain
        ! this was particularly problematic for ideal WRF simluation comparisons
        if ((size(bc%ext_winds%terrain,1)/=size(bc%next_domain%terrain,1)).or. &
            (size(bc%ext_winds%terrain,2)/=size(bc%next_domain%terrain,2))) then
            call geo_LUT(bc%next_domain, bc%ext_winds)
            allocate(temporary_terrain(size(bc%ext_winds%terrain,1),size(bc%ext_winds%terrain,2)))
            temporary_terrain = bc%ext_winds%terrain
            deallocate(bc%ext_winds%terrain)
            allocate(bc%ext_winds%terrain(size(bc%next_domain%terrain,1),size(bc%next_domain%terrain,2)))
            call geo_interp2d(bc%ext_winds%terrain,temporary_terrain,bc%ext_winds%geolut)
            deallocate(temporary_terrain)
        endif
        
        ! force all weight to be on the first x,y pair...
        ! this assumes the "external winds" file is on the exact same grid as the high res model grid
        ! probably better not to do this incase they are not on the same grid. This also assumes the first x,y pair
        ! is the closes grid cell.  This should be the case for the current geolut algorithm, but it is not required. 
        ! bc%ext_winds%u_geo%geolut%w(2:,:,:) = 0
        ! bc%ext_winds%u_geo%geolut%w(1,:,:) = 1
        ! bc%ext_winds%v_geo%geolut%w(2:,:,:) = 0
        ! bc%ext_winds%v_geo%geolut%w(1,:,:) = 1
        bc%ext_winds%dx = bc%next_domain%dx
        call setup_extwinds(bc%ext_winds)
    end subroutine init_ext_winds
    
    subroutine swap_z(bc)
        ! swap the lowres_z and z variables in bc
        type(bc_type), intent(inout) :: bc
        real,allocatable,dimension(:,:,:) :: tempz
        integer::nx,nz,ny, i
        nx = size(bc%lowres_z,1)
        nz = size(bc%lowres_z,2)
        ny = size(bc%lowres_z,3)
        allocate(tempz(nx,nz,ny))
        tempz = bc%lowres_z
        
        nx = size(bc%z,1)
        nz = size(bc%z,2)
        ny = size(bc%z,3)
        deallocate(bc%lowres_z)
        allocate(bc%lowres_z(nx,nz,ny))
        bc%lowres_z = bc%z

        nx = size(tempz,1)
        nz = size(tempz,2)
        ny = size(tempz,3)
        deallocate(bc%z)
        allocate(bc%z(nx,ny,nz))
        do i = 1,nz
            bc%z(:,:,i) = tempz(:,i,:)
        end do
        deallocate(tempz)
        
    end subroutine swap_z

    subroutine move_lut(inputgeo,outputgeo)
        ! move the contents of one geographic look up table into another
        ! moving implies deletion of the initial data, so inputgeo is destroyed
        type(geo_look_up_table), intent(inout) :: inputgeo,outputgeo
        integer::nx,ny,nz
        nx = size(inputgeo%x,1)
        ny = size(inputgeo%x,2)
        nz = size(inputgeo%x,3)
        
        allocate(outputgeo%x(nx,ny,nz))
        allocate(outputgeo%y(nx,ny,nz))
        allocate(outputgeo%w(nx,ny,nz))
        
        outputgeo%x = inputgeo%x
        outputgeo%y = inputgeo%y
        outputgeo%w = inputgeo%w
        
        call destroy_lut(inputgeo)
    end subroutine move_lut
    
    subroutine destroy_lut(geolut)
        ! deallocate all memory associated with a geographic look up table
        type(geo_look_up_table), intent(inout) :: geolut
        deallocate(geolut%x)
        deallocate(geolut%y)
        deallocate(geolut%w)
    end subroutine destroy_lut
    
!   initialize the boundary condiditions (init data structures and GEOLUT)
!   initializes external winds if necessary, adds low-res terrain to the high-res domain if desired
    subroutine init_bc(options,domain,boundary)
        implicit none
        type(options_type), intent(in) :: options
        type(domain_type), intent(inout):: domain
        type(bc_type), intent(inout):: boundary
        type(geo_look_up_table) :: u_temp_geo,v_temp_geo
        real, dimension(:,:,:), allocatable :: tempz
        integer::i,nx,ny,nz
            
        boundary%dx = options%dxlow
        ! set up base data
        call init_bc_data(options,boundary,domain)
        call init_domain(options,boundary%next_domain) !set up a domain to hold the forcing for the next time step
        
        ! create the geographic look up table used to calculate boundary forcing data
        write(*,*) "Setting up domain geographic Look Up Tables"
        ! NOTE: these first two geoLUTs are for translating from the mass grid to the U/V grids
        ! These are only used once to translate the terrain to those grids. 
        ! set up a look up table from low-res grid center to high-res u-offset coordinates
        call geo_LUT(domain%u_geo,boundary)
        call move_lut(boundary%geolut,u_temp_geo)
        ! set up a look up table from low-res grid center to high-res v-offset coordinates
        call geo_LUT(domain%v_geo,boundary)
        call move_lut(boundary%geolut,v_temp_geo)
        ! main geoLUTs
        call geo_LUT(domain,boundary)
        call geo_LUT(domain%u_geo,boundary%u_geo)
        call geo_LUT(domain%v_geo,boundary%v_geo)
        
        if (options%external_winds) then
            call init_ext_winds(options,boundary)
        endif
        
        nz = size(boundary%lowres_z,2)
        if (maxval(boundary%terrain)>maxval(boundary%lowres_z(:,1,:))) then
            write(*,*) "WARNING : Forcing Z levels are below the terrain provided. "
            write(*,*) "          For pressure level forcing data this is expected, but not well tested. "
        endif
        
        ! if we want vertical interpolations between forcing and model grid to be done from 
        ! height above ground level (and we should for wind?), then we need to remove the 
        ! topography from both grids first
        if (options%use_agl_height) then
            nz = size(boundary%lowres_z,2)
            
            ! First set the forcing data z into AGL coordinates
            do i = 1,nz
                boundary%lowres_z(:,i,:) = boundary%lowres_z(:,i,:)-boundary%terrain
            enddo
            ! Then put the model domain into AGL coordinates
            nz = size(domain%z,2)
            do i = 1,nz
                domain%z(:,i,:) = domain%z(:,i,:)-domain%terrain
            enddo

            ! finally, get the model terrain onto the u and v staggered grids
            call copy_z(domain,domain%u_geo,interpolate_dim=1)
            call copy_z(domain,domain%v_geo,interpolate_dim=3)

            ! this is the only place we need the domain z to be in AGL space, so revert back immediately
            do i = 1,nz
                domain%z(:,i,:) = domain%z(:,i,:)+domain%terrain
            enddo
        endif
        
        ! This section all happens in AGL space if use_agl_height=True
        ! interpolate lowres_z into the u_geo and v_geo data structures (on their staggered grids)
        nx = size(domain%u_geo%lat,1)
        ny = size(domain%u_geo%lat,2)
        nz = size(boundary%lowres_z,2)
        allocate(boundary%u_geo%z(nx,nz,ny))
        call geo_interp(boundary%u_geo%z,boundary%lowres_z,u_temp_geo,.false.)

        nx = size(domain%v_geo%lat,1)
        ny = size(domain%v_geo%lat,2)
        allocate(boundary%v_geo%z(nx,nz,ny))
        call geo_interp(boundary%v_geo%z,boundary%lowres_z,v_temp_geo,.false.)
        
        ! destroy the geographic look up tables to save space
        ! geoLUTs for forcing data are separate
        call destroy_lut(v_temp_geo)
        call destroy_lut(u_temp_geo)
        
        write(*,*) "Setting up vertical interpolation Look Up Tables"
        ! generate the vertical inteprolateion look up tables for the u and v grids
        call vLUT(domain%u_geo,boundary%u_geo)
        call vLUT(domain%v_geo,boundary%v_geo)
        
        ! restore the lowres_z variable to ASL space for the mass grid vLUT generation
        if (options%use_agl_height) then
            nz = size(boundary%lowres_z,2)
            do i = 1,nz
                boundary%lowres_z(:,i,:) = boundary%lowres_z(:,i,:)+boundary%terrain
            enddo
        endif

        ! Finally, take care of the mass grid in ASL space
        nx = size(domain%terrain,1)
        ny = size(domain%terrain,2)
        nz = size(boundary%lowres_z,2)
        allocate(boundary%z(nx,nz,ny))
        call geo_interp(boundary%z,boundary%lowres_z,boundary%geolut,.false.)

        if (options%debug) write(*,*) "Terrain z min=",minval(domain%terrain), "Terrain z max=", maxval(domain%terrain)
        if (options%debug) write(*,*) "Domain z min=", minval(domain%z), "Domain z max=", maxval(domain%z)
        if (options%debug) write(*,*) "Forcing z min=",minval(boundary%z), "Forcing z max=", maxval(boundary%z)
        call vLUT(domain,boundary)

        ! these are not longer needed and (without adjustments) are potentially unreliable (AGL vs ASL)
        deallocate(boundary%v_geo%z, boundary%u_geo%z,domain%u_geo%z,domain%v_geo%z)
        ! swaps z and lowres_z (one of the cases where pointers would make life a lot easier)
        call swap_z(boundary)
        
        ! after swap_z call, boundary%z is on the forcing data grid, boundary%lowres_z is interpolated to the ICAR grid
        ! Finally apply the vertical interpolation here as well
        ! this is a little annoying because we have to allocate a temporary. 
        nx = size(domain%terrain,1)
        ny = size(domain%terrain,2)
        nz = size(domain%p,2)
        allocate(tempz(nx,nz,ny))
        
        call vinterp(tempz,    & 
                     boundary%lowres_z, &
                     boundary%vert_lut)
        
        deallocate(boundary%lowres_z)
        allocate(boundary%lowres_z(nx,nz,ny))
        boundary%lowres_z = tempz
        deallocate(tempz)
        
    end subroutine init_bc
    
    subroutine welcome_message()
        implicit none
        
        write(*,*) ""
        write(*,*) "============================================================"
        write(*,*) "|                                                          |"
        write(*,*) "|  The Intermediate Complexity Atmospheric Research Model  |"
        write(*,*) "|                          (ICAR)                          |"
        write(*,*) "|                                                          |"
        write(*,*) "|   Developed at NCAR:                                     |"
        write(*,*) "|     The National Center for Atmospheric Research         |"
        write(*,*) "|     NCAR is sponsored by the National Science Foundation |"
        write(*,*) "|                                                          |"
        write(*,*) "|   Version: ",kVERSION_STRING,"                                         |"
        write(*,*) "|                                                          |"
        write(*,*) "============================================================"
        write(*,*) ""        
        
    end subroutine welcome_message
end module
