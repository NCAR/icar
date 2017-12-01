module domain_interface
  use options_interface,        only : options_t
  use exchangeable_interface,   only : exchangeable_t
  use grid_interface,           only : grid_t
  use variable_interface,       only : variable_t
  use meta_data_interface,      only : meta_data_t
  use io_routines,              only : io_read
  use time_object,              only : Time_type
  implicit none

  private
  public :: domain_t

  type domain_t
    type(meta_data_t)    :: info
    type(grid_t)         :: grid, u_grid, v_grid, grid2d

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
    real,   allocatable :: exner                (:,:,:)
    real,   allocatable :: density              (:,:,:)
    real,   allocatable :: pressure             (:,:,:)
    real,   allocatable :: pressure_interface   (:,:,:)
    real,   allocatable :: temperature          (:,:,:)
    real,   allocatable :: z                    (:,:,:)
    real,   allocatable :: dz_interface         (:,:,:)
    real,   allocatable :: z_interface          (:,:,:)
    real,   allocatable :: dz_mass              (:,:,:)
    real,   allocatable :: graupel                  (:,:)
    real,   allocatable :: accumulated_precipitation(:,:)
    integer,allocatable :: precipitation_bucket     (:,:)
    real,   allocatable :: accumulated_snowfall     (:,:)
    integer,allocatable :: snowfall_bucket          (:,:)
    real,   allocatable :: longwave             (:,:)
    real,   allocatable :: shortwave            (:,:)
    real,   allocatable :: terrain              (:,:)
    integer,allocatable :: land_cover_type      (:,:)
    real,   allocatable :: vegetation_fraction  (:,:,:)
    real,   allocatable :: lai                  (:,:)
    real,   allocatable :: canopy_water         (:,:)
    real,   allocatable :: snow_water_equivalent(:,:)
    real,   allocatable :: skin_temperature     (:,:)
    real,   allocatable :: soil_water_content   (:,:,:)
    real,   allocatable :: soil_temperature     (:,:,:)
    integer,allocatable :: land_mask            (:,:)
    real,   allocatable :: latitude             (:,:)
    real,   allocatable :: longitude            (:,:)
    real,   allocatable :: u_latitude           (:,:)
    real,   allocatable :: u_longitude          (:,:)
    real,   allocatable :: v_latitude           (:,:)
    real,   allocatable :: v_longitude          (:,:)


    ! these coarrays are used to send all data to/from a master image for IO... ?
    ! For now this will be taken care of in the boundary conditions object
    ! real, allocatable :: transfer_3d(:,:,:)[:]
    ! real, allocatable :: transfer_2d(:,:)[:]

    ! Array listing variables to advect with pointers to local data
    type(variable_t), allocatable :: advected_species(:)

    ! contains the size of the domain (or the local tile?)
    integer :: nx, ny, nz, nx_global, ny_global
    integer :: ximg, ximages, yimg, yimages

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

  end type

  integer, parameter :: space_dimension=3

  interface

    ! Set default component values
    module subroutine init(this, options)
      implicit none
      class(domain_t), intent(inout) :: this
      class(options_t),intent(inout) :: options
    end subroutine

    module subroutine var_request(this, options)
        implicit none
        class(domain_t), intent(inout) :: this
        class(options_t),intent(inout) :: options
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

  end interface

end module
