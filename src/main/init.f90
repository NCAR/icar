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
module initialization
    use data_structures
    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use boundary_interface, only : boundary_t
    use microphysics,               only : mp_init
    use advection,                  only : adv_init
    use radiation,                  only : radiation_init
    use convection,                 only : init_convection
    use planetary_boundary_layer,   only : pbl_init
    use land_surface,               only : lsm_init

    use mod_atm_utilities,          only : init_atm_utilities
    use wind,                       only : update_winds

    ! use io_routines,                only : io_read, &
    !                                        io_write3d,io_write3di, io_write
    ! use geo,                        only : geo_LUT, geo_interp, geo_interp2d, standardize_coordinates
    ! use vertical_interpolation,     only : vLUT, vinterp
    ! use wind,                       only : init_winds
    ! use initialize_options,         only : init_options
    ! use string,                     only : str


    implicit none
    private
    public::init_model, init_physics

contains
    subroutine init_model(options,domain,boundary)
        implicit none
        type(options_t), intent(inout) :: options
        type(domain_t),  intent(inout) :: domain
        type(boundary_t),intent(inout) :: boundary

        integer :: omp_get_max_threads, num_threads

#if defined(_OPENMP)
        num_threads = omp_get_max_threads()
#else
        num_threads = 1
#endif
        if (this_image()==1) call welcome_message()

        if (this_image()==1) then
            write(*,*) "  Number of coarray image:",num_images()
            write(*,*) "  Max number of OpenMP Threads:",num_threads
        endif

        ! read in options file
        if (this_image()==1) write(*,*) "Initializing Options"
        call options%init()

        if (this_image()==1) write(*,*) "Initializing Domain"
        call domain%init(options)

        if (this_image()==1) write(*,*) "Initializing boundary condition data structure"
        call boundary%init(options)

        if (this_image()==1) write(*,*) "Reading Initial conditions from boundary dataset"
        call domain%get_initial_conditions(boundary, options)

        if (this_image()==1) write(*,*) "Updating initial winds"
        call update_winds(domain, options)

        ! initialize the atmospheric helper utilities
        call init_atm_utilities(options)

        call init_physics(options, domain)

        ! call setup_bias_correction(options,domain)
        if (this_image()==1) write(*,'(/ A)') "Finished basic initialization"
        if (this_image()==1) write(*,'(A /)') "---------------------------------------"

    end subroutine init_model

    subroutine init_physics(options, domain)
        implicit none
        type(options_t), intent(inout) :: options
        type(domain_t),  intent(inout) :: domain

        ! initialize microphysics code (e.g. compute look up tables in Thompson et al)
        call mp_init(options) !this could easily be moved to init_model...

        call init_convection(domain,options)

        call pbl_init(domain,options)

        call radiation_init(domain,options)

        call lsm_init(domain,options)

        call adv_init(domain,options)

    end subroutine init_physics

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


! ------------------------------------------------------------------------------------------
!-==== Model Domain Section ====
!
! Begining of section focused on allocating and initializing the model domain data structures
!
! ------------------------------------------------------------------------------------------

!   Allow running over a sub-domain, by removing the outer N grid cells from all sides of the domain (lat,lon,terrain)
    ! subroutine remove_edges(domain,edgesize)
    !     implicit none
    !     type(domain_t), intent(inout) :: domain
    !     integer, intent(in)::edgesize
    !     integer::nx1,ny1,nx2,ny2,nz
    !     real,allocatable,dimension(:,:)::temp_data
    !
    !     nx1 = size(domain%lat,1)
    !     ny1 = size(domain%lat,2)
    !     nx2 = nx1-(edgesize*2)
    !     ny2 = ny1-(edgesize*2)
    !
    !     allocate(temp_data(nx1,ny1))
    !
    !     temp_data = domain%lat
    !     deallocate(domain%lat)
    !     allocate(domain%lat(nx2,ny2))
    !     domain%lat = temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)
    !
    !     temp_data = domain%lon
    !     deallocate(domain%lon)
    !     allocate(domain%lon(nx2,ny2))
    !
    !     domain%lon = temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)
    !
    !     temp_data = domain%terrain
    !     deallocate(domain%terrain)
    !     allocate(domain%terrain(nx2,ny2))
    !     domain%terrain = temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)
    !
    !     temp_data = domain%landmask
    !     deallocate(domain%landmask)
    !     allocate(domain%landmask(nx2,ny2))
    !     domain%landmask = temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)
    !
    !
    !     deallocate(temp_data)
    !
    !     nx1 = size(domain%u_geo%lat,1)
    !     ny1 = size(domain%u_geo%lat,2)
    !     nx2 = nx1-(edgesize*2)
    !     ny2 = ny1-(edgesize*2)
    !
    !     allocate(temp_data(nx1,ny1))
    !     temp_data = domain%u_geo%lat
    !     deallocate(domain%u_geo%lat)
    !     allocate(domain%u_geo%lat(nx2,ny2))
    !     domain%u_geo%lat = temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)
    !
    !     temp_data = domain%u_geo%lon
    !     deallocate(domain%u_geo%lon)
    !     allocate(domain%u_geo%lon(nx2,ny2))
    !     domain%u_geo%lon = temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)
    !
    !     deallocate(temp_data)
    !
    !     nx1 = size(domain%v_geo%lat,1)
    !     ny1 = size(domain%v_geo%lat,2)
    !     nx2 = nx1-(edgesize*2)
    !     ny2 = ny1-(edgesize*2)
    !
    !     allocate(temp_data(nx1,ny1))
    !     temp_data = domain%v_geo%lat
    !     deallocate(domain%v_geo%lat)
    !     allocate(domain%v_geo%lat(nx2,ny2))
    !     domain%v_geo%lat = temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)
    !
    !     temp_data = domain%v_geo%lon
    !     deallocate(domain%v_geo%lon)
    !     allocate(domain%v_geo%lon(nx2,ny2))
    !     domain%v_geo%lon = temp_data(1+edgesize:nx1-edgesize,1+edgesize:ny1-edgesize)
    !
    !     deallocate(temp_data)
    ! end subroutine remove_edges

!   allocate all arrays in domain

! interpolate intput%z to output%z assuming that input has one less grid cell
! in the interpolate_dim dimension
    ! subroutine copy_z(input,output,interpolate_dim)
    !     implicit none
    !     class(interpolable_type), intent(in) :: input
    !     class(interpolable_type), intent(inout) :: output
    !     integer,intent(in)::interpolate_dim
    !
    !     integer::nxi,nyi,nzi,nxo,nyo,nzo
    !
    !     ! dimensions of the input data
    !     nxi = size(input%z,1)
    !     nzi = size(input%z,2)
    !     nyi = size(input%z,3)
    !     ! dimensions of the output data
    !     nxo = size(output%lat,1)
    !     nzo = nzi
    !     nyo = size(output%lat,2)
    !
    !     if (allocated(output%z)) then
    !         deallocate(output%z)
    !     endif
    !     allocate(output%z(nxo,nzo,nyo))
    !     if (interpolate_dim==1) then
    !         if ((nxo/=(nxi+1)).or.(nyo/=nyi).or.(nzo/=nzi)) then
    !             stop "Error copying z levels from mass grid to U grid, dimensions don't match"
    !         endif
    !         output%z(2:nxo-1,:,:) = (input%z(1:nxi-1,:,:) + input%z(2:nxi,:,:))/2
    !         output%z(1,:,:) = input%z(1,:,:)
    !         output%z(nxo,:,:) = input%z(nxi,:,:)
    !     else if (interpolate_dim==3) then
    !         if ((nyo/=(nyi+1)).or.(nxo/=nxi).or.(nzo/=nzi)) then
    !             stop "Error copying z levels from mass grid to V grid, dimensions don't match"
    !         endif
    !         output%z(:,:,2:nyo-1) = (input%z(:,:,1:nyi-1) + input%z(:,:,2:nyi))/2
    !         output%z(:,:,1) = input%z(:,:,1)
    !         output%z(:,:,nyo) = input%z(:,:,nyi)
    !     else
    !         write(*,*) "Can not interpolate z data over z dimension"
    !     endif
    !
    ! end subroutine copy_z

    ! subroutine setup_bias_correction(options, domain)
    !     implicit none
    !     type(options_t), intent(in) :: options
    !     type(domain_t), intent(inout):: domain
    !
    !     if (options%use_bias_correction) then
    !         call io_read(options%bias_options%filename, options%bias_options%rain_fraction_var, domain%rain_fraction)
    !     endif
    ! end subroutine setup_bias_correction
    !
    ! subroutine init_domain_land(domain,options)
    !     implicit none
    !     type(options_t), intent(in) :: options
    !     type(domain_t), intent(inout):: domain
    !     ! these are temporary variables used for IO before storing data in the domain data structure
    !     integer,dimension(:,:),allocatable :: temp_idata
    !     real,dimension(:,:,:),allocatable :: temp_rdata
    !     real,dimension(:,:),allocatable :: temp_rdata_2d
    !
    !     integer:: buffer,nx,ny,nz, i
    !
    !     buffer = options%buffer
    !     nx = size(domain%veg_type,1)+buffer*2
    !     ny = size(domain%veg_type,2)+buffer*2
    !     nz = size(domain%soil_t,2)
    !
    !
    !     ! Veg cover fraction = 2D/3D real
    !     if (options%vegfrac_var.ne."") then
    !         ! if we are supposed to read a veg fraction for each month then it will be a 3D variable
    !         if (options%lsm_options%monthly_vegfrac) then
    !             call io_read(options%init_conditions_file,options%vegfrac_var,temp_rdata,1)
    !             domain%vegfrac = temp_rdata(1+buffer:nx-buffer,1+buffer:ny-buffer,:)  ! subset the data by buffer grid cells
    !             deallocate(temp_rdata)
    !         else
    !             ! otherwise we just read a 2D variable and store it in the first time entry in the 3D vegfrac field
    !             call io_read(options%init_conditions_file,options%vegfrac_var,temp_rdata_2d,1)
    !             domain%vegfrac(:,:,1) = temp_rdata_2d(1+buffer:nx-buffer,1+buffer:ny-buffer)  ! subset the data by buffer grid cells
    !             deallocate(temp_rdata_2d)
    !         endif
    !     endif
    !     ! if the input vegetation fraction is a fraction (0-1)
    !     ! then convert it to a percent
    !     if (maxval(domain%vegfrac)<=1) then
    !         domain%vegfrac = domain%vegfrac*100
    !     endif
    !
    !     ! Veg TYPE = 2D integer
    !     if (options%vegtype_var.ne."") then
    !         call io_read(options%init_conditions_file,options%vegtype_var,temp_idata,1)
    !         domain%veg_type = temp_idata(1+buffer:nx-buffer,1+buffer:ny-buffer)  ! subset the data by buffer grid cells
    !         deallocate(temp_idata)
    !     endif
    !
    !     ! Soil TYPE = 2D integer
    !     if (options%soiltype_var.ne."") then
    !         call io_read(options%init_conditions_file,options%soiltype_var,temp_idata,1)
    !         domain%soil_type = temp_idata(1+buffer:nx-buffer,1+buffer:ny-buffer)  ! subset the data by buffer grid cells
    !         deallocate(temp_idata)
    !     endif
    !
    !     ! Soil Volumetric Water Content = 3D real
    !     if (options%soil_vwc_var.ne."") then
    !         call io_read(options%init_conditions_file,options%soil_vwc_var,temp_rdata,1)
    !         domain%soil_vwc = reshape(temp_rdata(1+buffer:nx-buffer,1+buffer:ny-buffer,1:nz) &    ! subset the data by buffer grid cells
    !                                 ,[nx-buffer*2,nz,ny-buffer*2],order=[1,3,2])                ! and reshape to move the z axis to the middle
    !         deallocate(temp_rdata)
    !     endif
    !
    !     ! Soil Temperature = 3D real
    !     if (options%soil_t_var.ne."") then
    !         call io_read(options%init_conditions_file,options%soil_t_var,temp_rdata,1)
    !         domain%soil_t = reshape(temp_rdata(1+buffer:nx-buffer,1+buffer:ny-buffer,1:nz) &  ! subset the data by buffer grid cells
    !                                 ,[nx-buffer*2,nz,ny-buffer*2],order=[1,3,2])            ! and reshape to move the z axis to the middle
    !         deallocate(temp_rdata)
    !     endif
    !
    !     ! Deep Soil Temperature = 2D real
    !     if (options%soil_deept_var.ne."") then
    !         call io_read(options%init_conditions_file,options%soil_deept_var,temp_rdata_2d,1)
    !         domain%soil_tdeep = temp_rdata_2d(1+buffer:nx-buffer,1+buffer:ny-buffer)  ! subset the data by buffer grid cells
    !         where(domain%soil_tdeep<200) domain%soil_tdeep = 200 ! mitigates zeros that can cause problems
    !         if (options%soil_t_var=="") then
    !             write(*,*) "Missing explicit soil T, using deep soil T for all depths"
    !             do i = 1,nz
    !                 domain%soil_t(:,i,:) = domain%soil_tdeep
    !             end do
    !         endif
    !         deallocate(temp_rdata_2d)
    !     endif
    !
    ! end subroutine init_domain_land

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
    ! subroutine interpolate_in_x(data_in, data_out)
    !     implicit none
    !     real, dimension(:,:), intent(in) :: data_in
    !     real, dimension(:,:), allocatable, intent(out) :: data_out
    !
    !     integer :: nx, ny
    !
    !     nx = size(data_in,1)
    !     ny = size(data_in,2)
    !     if (.not.allocated(data_out)) then
    !         allocate(data_out(nx+1,ny))
    !     else
    !         if ((size(data_out,1)/=nx+1).or.(size(data_out,2)/=ny)) then
    !             deallocate(data_out)
    !             allocate(data_out(nx+1,ny))
    !         endif
    !     endif
    !
    !     ! simple bilinear inteprolation (included extrapolation one grid cell out)
    !     data_out(2:nx,:) = (data_in(2:,:) + data_in(1:nx-1,:))/2
    !     data_out(1,:)    = data_in(1,:)*2 - data_in(2,:)
    !     data_out(nx+1,:) = data_in(nx,:)*2 - data_in(nx-1,:)
    !     ! result is data_out
    ! end subroutine interpolate_in_x

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
    ! subroutine interpolate_in_y(data_in, data_out)
    !     implicit none
    !     real, dimension(:,:), intent(in) :: data_in
    !     real, dimension(:,:), allocatable, intent(out) :: data_out
    !
    !     integer :: nx, ny
    !
    !     nx = size(data_in,1)
    !     ny = size(data_in,2)
    !     if (.not.allocated(data_out)) then
    !         allocate(data_out(nx,ny+1))
    !     else
    !         if ((size(data_out,1)/=nx).or.(size(data_out,2)/=ny+1)) then
    !             deallocate(data_out)
    !             allocate(data_out(nx,ny+1))
    !         endif
    !     endif
    !
    !     ! simple bilinear inteprolation (included extrapolation one grid cell out)
    !     data_out(:,2:ny) = (data_in(:,2:) + data_in(:,1:ny-1))/2
    !     data_out(:,1)    = data_in(:,1)*2 - data_in(:,2)
    !     data_out(:,ny+1) = data_in(:,ny)*2 - data_in(:,ny-1)
    !     ! result is data_out
    ! end subroutine interpolate_in_y

!   initialize the domain e.g. lat,lon,terrain, 3D z coordinate
!     subroutine init_domain(options, domain)
!         implicit none
!         type(options_t), intent(in) :: options
!         type(domain_t), intent(inout):: domain
!         real,dimension(:,:,:),allocatable::temporary_z
!         integer:: ny,nz,nx,i,buf
!
! !       these are the only required variables on a high-res grid, lat, lon, and terrain elevation
!         call io_read(options%init_conditions_file,options%hgt_hi,domain%terrain,1)
!         call io_read(options%init_conditions_file,options%lat_hi,domain%lat,1)
!         call io_read(options%init_conditions_file,options%lon_hi,domain%lon,1)
!         call standardize_coordinates(domain)
!
!         !  because u/vlat/lon_hi are optional, we calculated them by interplating from the mass lat/lon if necessary
!         if (options%ulat_hi/="") then
!             call io_read(options%init_conditions_file,options%ulat_hi,domain%u_geo%lat,1)
!         else
!             call interpolate_in_x(domain%lat,domain%u_geo%lat)
!         endif
!         if (options%ulon_hi/="") then
!             call io_read(options%init_conditions_file,options%ulon_hi,domain%u_geo%lon,1)
!         else
!             call interpolate_in_x(domain%lon,domain%u_geo%lon)
!         endif
!
!         call standardize_coordinates(domain%u_geo)
!
!         if (options%vlat_hi/="") then
!             call io_read(options%init_conditions_file,options%vlat_hi,domain%v_geo%lat,1)
!         else
!             call interpolate_in_y(domain%lat,domain%v_geo%lat)
!         endif
!         if (options%vlon_hi/="") then
!             call io_read(options%init_conditions_file,options%vlon_hi,domain%v_geo%lon,1)
!         else
!             call interpolate_in_y(domain%lon,domain%v_geo%lon)
!         endif
!         call standardize_coordinates(domain%v_geo)
!
!         if (options%landvar/="") then
!             call io_read(options%init_conditions_file,options%landvar,domain%landmask,1)
!             where(domain%landmask==0) domain%landmask = kLC_WATER
!         else
!             nx = size(domain%lat,1)
!             ny = size(domain%lat,2)
!             allocate(domain%landmask(nx,ny))
!             domain%landmask = kLC_LAND !if we weren't supplied a landmask field, assume all is land (what we care about anyway)
!         endif
!
!         ! remove n grid cells from all sides of the domain if requested
!         if(options%buffer>0) then
!             write(*,*) "  Removing ",trim(str(options%buffer))," gridcells from edges of domain"
!             call remove_edges(domain,options%buffer)
!         endif
!
!         ! use the lat variable to define the x and y dimensions for all other variables
!         nx = size(domain%lat,1)
!         ny = size(domain%lat,2)
!         ! assumes nz is defined in the options
!         nz = options%nz
!
!         ! if a 3d grid was also specified, then read those data in
!         if ((options%readz).and.(options%ideal).and.(options%zvar.ne."")) then
!
!             stop "Reading Z from an external file is not currently supported, use a fixed dz"
!
!             call io_read(options%init_conditions_file,options%zvar, domain%z)
!             ! dz also has to be calculated from the 3d z file
!             buf = options%buffer
!             allocate(domain%dz(nx,nz,ny))
!             allocate(temporary_z(nx,nz,ny))
!             do i = 1,nz-1
!                 domain%dz(:,i,:) = domain%z(buf+1:nx+buf,buf+1:ny+buf,i+1)-domain%z(buf+1:nx+buf,buf+1:ny+buf,i)
!                 temporary_z(:,i,:) = domain%z(buf+1:nx+buf,buf+1:ny+buf,i)
!             enddo
!             temporary_z(:,nz,:) = domain%z(buf+1:nx+buf,buf+1:ny+buf,nz)
!             domain%dz(:,nz,:) = domain%dz(:,nz-1,:)
!             deallocate(domain%z)
!             allocate(domain%z(nx,nz,ny))
!             domain%z = temporary_z
!             deallocate(temporary_z)
!
!         else
!             ! otherwise, set up the z grid to be evenly spaced in z using the terrain +dz/2 for the base
!             ! and z[i-1]+(dz[i-1]+dz[i])/2 for the rest
!             allocate(domain%z(nx,nz,ny))
!             allocate(domain%z_inter(nx,nz,ny))
!             allocate(domain%dz(nx,nz,ny))
!             allocate(domain%dz_inter(nx,nz,ny))
!             allocate(domain%z_layers(nz))
!             allocate(domain%z_interface_layers(nz+1))
!             ! lowest model level is half of the lowest dz above the land surface
!             domain%z(:,1,:) = domain%terrain+options%dz_levels(1)/2
!             domain%z_inter(:,1,:) = domain%terrain
!             ! dz between the mass grid points is half of each dz
!             domain%dz(:,1,:) = (options%dz_levels(1) + options%dz_levels(2))/2
!             ! dz between the interface points = dz = thickness of mass grid cells
!             domain%dz_inter(:,1,:) = options%dz_levels(1)
!             ! create 1D arrays to store the height of each level and interface above the terrain too
!             domain%z_layers(1) = options%dz_levels(1)/2
!             domain%z_interface_layers(1) = 0
!             domain%z_interface_layers(2) = options%dz_levels(1)/2
!             do i = 2,nz
!                 domain%z_interface_layers(i+1) = domain%z_interface_layers(i) + options%dz_levels(i)
!                 domain%z_layers(i)    = domain%z_interface_layers(i) + options%dz_levels(i) / 2
!
!                 ! although these are almost identical, they have topography added to the first layer
!                 domain%z(:,i,:)       = domain%z(:,i-1,:)       + (options%dz_levels(i)+options%dz_levels(i-1))/2
!                 domain%z_inter(:,i,:) = domain%z_inter(:,i-1,:) +  options%dz_levels(i)
!
!                 ! dz between the interface points = dz = thickness of mass grid cells
!                 domain%dz_inter(:,i,:) = options%dz_levels(i)
!                 ! dz between the mass grid points is half of each dz
!                 if (i<nz) then
!                     domain%dz(:,i,:) = (options%dz_levels(i) + options%dz_levels(i+1))/2
!                 endif
!             enddo
!             domain%dz(:,nz,:) = options%dz_levels(nz)
!
!         endif
!         call copy_z(domain,domain%u_geo,interpolate_dim=1)
!         call copy_z(domain,domain%v_geo,interpolate_dim=3)
!
!         ! all other variables should be allocated and initialized to 0
!         call domain_allocation(domain,options,nx,nz,ny)
!         ! initializing land
!         call init_domain_land(domain,options)
!
!         ! store dx in domain as well as options, read as an option, but it is more appropriate in domain
!         domain%dx = options%dx
!
!         call init_winds(domain,options)
!     end subroutine init_domain

end module
