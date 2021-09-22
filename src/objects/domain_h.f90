module domain_interface
  use iso_c_binding
  use options_interface,        only : options_t
  use boundary_interface,       only : boundary_t
  use exchangeable_interface,   only : exchangeable_t
  use grid_interface,           only : grid_t
  use variable_interface,       only : variable_t
  use variable_dict_interface,  only : var_dict_t
  use meta_data_interface,      only : meta_data_t
  use time_object,              only : Time_type
  use time_delta_object,        only : time_delta_t
  use data_structures,          only : interpolable_type, tendencies_type
  implicit none

  private
  public :: domain_t

  type domain_t
    type(meta_data_t)    :: info
    type(grid_t)         :: grid,   u_grid,   v_grid
    type(grid_t)         :: grid2d, u_grid2d, v_grid2d
    type(grid_t)         :: u_grid2d_ext, v_grid2d_ext!, grid2d_ext ! extended grids for u and v fields pre smoothing (grid_2d_ext is for SLEVE topography smoothing)
    type(grid_t)         :: grid_monthly, grid_soil

    type(Time_type) :: model_time

    ! note that not all variables are allocated at runtime, physics packages must request a variable be created
    ! though variables considered "required" are requested by the domain object itself (e.g. terrain)
    ! core model species to be advected

    ! wind field to control advection
    type(exchangeable_t) :: u
    type(exchangeable_t) :: v
    type(exchangeable_t) :: w

    type(exchangeable_t) :: water_vapor
    type(exchangeable_t) :: potential_temperature
    type(exchangeable_t) :: cloud_water_mass
    type(exchangeable_t) :: cloud_number
    type(exchangeable_t) :: cloud_ice_mass
    type(exchangeable_t) :: cloud_ice_number
    type(exchangeable_t) :: rain_mass
    type(exchangeable_t) :: rain_number
    type(exchangeable_t) :: snow_mass
    type(exchangeable_t) :: snow_number
    type(exchangeable_t) :: graupel_mass
    type(exchangeable_t) :: graupel_number

    ! other model variables (not advected)
    type(variable_t) :: exner
    type(variable_t) :: density
    type(variable_t) :: pressure
    type(variable_t) :: pressure_interface
    type(variable_t) :: temperature
    type(variable_t) :: z
    type(variable_t) :: dz_interface
    type(variable_t) :: z_interface
    type(variable_t) :: dz_mass
    type(variable_t) :: nsquared
    type(variable_t) :: graupel
    type(variable_t) :: accumulated_precipitation
    integer,allocatable :: precipitation_bucket(:,:)
    type(variable_t) :: accumulated_convective_pcp
    integer,allocatable :: cu_precipitation_bucket(:,:)
    type(variable_t) :: accumulated_snowfall
    integer,allocatable :: snowfall_bucket(:,:)
    type(variable_t) :: cloud_fraction
    type(variable_t) :: longwave
    type(variable_t) :: shortwave
    type(variable_t) :: terrain
    type(variable_t) :: forcing_terrain  ! BK 05/2020: The forcing terrain interpolated 2d to the hi-res grid. In order to calculate difference in slope
        type(variable_t) :: forcing_terrain2 ! test 9-6-2020
        ! type(variable_t) :: forcing_terrain_u1 ! test 9-6-2020
    type(variable_t) :: u_10m
    type(variable_t) :: v_10m
    type(variable_t) :: temperature_2m
    type(variable_t) :: humidity_2m
    type(variable_t) :: surface_pressure
    type(variable_t) :: longwave_up
    type(variable_t) :: ground_heat_flux
    type(variable_t) :: sensible_heat
    type(variable_t) :: latent_heat
    integer,allocatable :: veg_type(:,:)
    integer,allocatable :: soil_type(:,:)
    type(variable_t) :: roughness_z0
    type(variable_t) :: vegetation_fraction
    type(variable_t) :: lai
    type(variable_t) :: canopy_water
    type(variable_t) :: snow_water_equivalent
    type(variable_t) :: snow_height
    type(variable_t) :: skin_temperature
    type(variable_t) :: sst
    type(variable_t) :: soil_water_content
    type(variable_t) :: soil_temperature
    type(variable_t) :: soil_totalmoisture
    type(variable_t) :: soil_deep_temperature

    integer,allocatable :: land_mask(:,:)
    type(variable_t) :: latitude
    type(variable_t) :: longitude
    type(variable_t) :: u_latitude
    type(variable_t) :: u_longitude
    type(variable_t) :: v_latitude
    type(variable_t) :: v_longitude

    type(variable_t) :: w_real
    type(variable_t) :: u_mass
    type(variable_t) :: v_mass

    type(tendencies_type) :: tend

    type(var_dict_t) :: variables_to_force

    type(interpolable_type) :: geo
    type(interpolable_type) :: geo_u
    type(interpolable_type) :: geo_v

    real :: smooth_height, dx
    integer :: nsmooth

    complex(C_DOUBLE_COMPLEX),  allocatable :: terrain_frequency(:,:) ! FFT(terrain)
    double precision,           allocatable :: costheta(:,:)
    double precision,           allocatable :: sintheta(:,:)
    real,                       allocatable :: advection_dz(:,:,:)
    ! store the ratio between the average dz and each grid cells topographically modified dz (for space varying dz only)
    real,                       allocatable :: jacobian(:,:,:)
    real,                       allocatable :: jacobian_u(:,:,:)
    real,                       allocatable :: jacobian_v(:,:,:)
    real,                       allocatable :: jacobian_w(:,:,:)
    real,                       allocatable :: zr_u(:,:,:)
    real,                       allocatable :: zr_v(:,:,:)
    real,                       allocatable :: dzdx(:,:,:) ! change in height with change in x/y position (used to calculate w_real vertical motions)
    real,                       allocatable :: dzdy(:,:,:) ! change in height with change in x/y position (used to calculate w_real vertical motions)
    ! BK 2020/05
    real,                       allocatable :: delta_dzdx(:,:,:) ! change in height difference (between hi and lo-res data) with change in x/y position (used to calculate w_real vertical motions)
    real,                       allocatable :: delta_dzdy(:,:,:) ! change in height difference (between hi and lo-res data) with change in x/y position (used to calculate w_real vertical motions)
    real,                       allocatable :: zfr_u(:,:,:)     ! ratio between z levels (on grid)
    real,                       allocatable :: zfr_v(:,:,:)
    real,                       allocatable :: terrain_u(:,:)
    real,                       allocatable :: terrain_v(:,:)
    real,                       allocatable :: forcing_terrain_u(:,:)
    real,                       allocatable :: forcing_terrain_v(:,:)
    real,                       allocatable :: forc(:,:)
    real,                       allocatable :: h1(:,:)     ! the large-scale terrain (h1) for the SLEVE coordinate (achieved by smoothin the org terrain)
    real,                       allocatable :: h2(:,:)     ! the small-scale terrain (h2) for the SLEVE coordinate (difference org and h1 terrain)
    real,                       allocatable :: h1_u(:,:)     ! the large-scale terrain (h1) on the u grid
    real,                       allocatable :: h1_v(:,:)     ! the large-scale terrain (h1) on the v grid
    real,                       allocatable :: h2_u(:,:)     ! the small-scale terrain (h2) on the u grid
    real,                       allocatable :: h2_v(:,:)     ! the small-scale terrain (h2) on the v grid

    real,                       allocatable :: dz_scl(:)  ! the scaled dz levels, required for delta terrain calculation

    real,                       allocatable :: ustar(:,:)
    real,                       allocatable :: znu(:)
    real,                       allocatable :: znw(:)

    ! these data are stored on the domain wide grid even if this process is only looking at a subgrid
    ! these variables are necessary with linear winds, especially with spatially variable dz, to compute the LUT
    real,                       allocatable :: global_terrain(:,:)
    real,                       allocatable :: global_z_interface(:,:,:)
    real,                       allocatable :: global_dz_interface(:,:,:)
    real,                       allocatable :: global_jacobian(:,:,:)


    ! these coarrays are used to send all data to/from a master image for IO... ?
    ! For now this will be taken care of in the boundary conditions object
    ! real, allocatable :: transfer_3d(:,:,:)[:]
    ! real, allocatable :: transfer_2d(:,:)[:]

    ! Array listing variables to advect with pointers to local data
    type(variable_t), allocatable :: advected_species(:)

    ! contains the size of the domain (or the local tile?)
    integer :: nx, ny, nz, nx_global, ny_global
    integer :: ximg, ximages, yimg, yimages

    logical :: north_boundary = .True.
    logical :: south_boundary = .True.
    logical :: east_boundary = .True.
    logical :: west_boundary = .True.

    ! store the start (s) and end (e) for the i,j,k dimensions
    integer ::  ids,ide, jds,jde, kds,kde, & ! for the entire model domain    (d)
                ims,ime, jms,jme, kms,kme, & ! for the memory in these arrays (m)
                its,ite, jts,jte, kts,kte    ! for the data tile to process   (t)

  contains
    procedure :: init
    procedure :: var_request

    procedure :: halo_send
    procedure :: halo_retrieve
    procedure :: halo_exchange
    procedure :: enforce_limits

    procedure :: get_initial_conditions
    procedure :: interpolate_forcing
    procedure :: interpolate_external
    procedure :: update_delta_fields
    procedure :: apply_forcing

    procedure :: calculate_delta_terrain

  end type

  integer, parameter :: space_dimension=3

  interface

    ! Set default component values
    module subroutine init(this, options)
        implicit none
        class(domain_t), intent(inout) :: this
        type(options_t), intent(inout) :: options
    end subroutine

    ! read initial atmospheric conditions from forcing data
    module subroutine get_initial_conditions(this, forcing, options, external_conditions)
        implicit none
        class(domain_t),  intent(inout) :: this
        type(boundary_t), intent(inout) :: forcing
        type(boundary_t), intent(inout), optional :: external_conditions  ! external data such as SWE
        type(options_t),  intent(in)    :: options
    end subroutine

    module subroutine interpolate_external(this, external_conditions, options)
        implicit none
        class(domain_t),  intent(inout) :: this
        type(boundary_t), intent(in)    :: external_conditions
        type(options_t),  intent(in)    :: options
    end subroutine

    module subroutine interpolate_forcing(this, forcing, update)
        implicit none
        class(domain_t),  intent(inout) :: this
        type(boundary_t), intent(in)    :: forcing
        logical,          intent(in),   optional :: update
    end subroutine

    module subroutine var_request(this, options)
        implicit none
        class(domain_t), intent(inout) :: this
        type(options_t), intent(inout)  :: options
    end subroutine

    module subroutine halo_send(this)
        implicit none
        class(domain_t), intent(inout) :: this
    end subroutine

    module subroutine halo_retrieve(this)
        implicit none
        class(domain_t), intent(inout) :: this
    end subroutine


    ! Exchange subdomain boundary information
    module subroutine halo_exchange(this)
        implicit none
        class(domain_t), intent(inout) :: this
    end subroutine

    ! Make sure no hydrometeors are getting below 0
    module subroutine enforce_limits(this)
        implicit none
        class(domain_t), intent(inout) :: this
    end subroutine

    module subroutine update_delta_fields(this, dt)
        implicit none
        class(domain_t),    intent(inout) :: this
        type(time_delta_t), intent(in)    :: dt
    end subroutine

    module subroutine apply_forcing(this, dt)
        implicit none
        class(domain_t),    intent(inout) :: this
        type(time_delta_t), intent(in)    :: dt
    end subroutine

    module subroutine calculate_delta_terrain(this, forcing, options)
        implicit none
        class(domain_t), intent(inout) :: this
        type(boundary_t), intent(in)    :: forcing
        type(options_t), intent(in) :: options
    end subroutine

  end interface

end module
