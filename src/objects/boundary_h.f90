module boundary_interface
  use data_structures,          only : linearizable_type, domain_type
  use options_interface,        only : options_t
  ! use exchangeable_interface,   only : exchangeable_t
  use grid_interface,           only : grid_t
  ! use variable_interface,       only : variable_t
  use meta_data_interface,      only : meta_data_t
  use time_object,              only : Time_type
  use time_delta_object,        only : time_delta_t

  implicit none

  private
  public :: boundary_t

  ! ------------------------------------------------
  ! boundary conditions type, must be linearizable so we can remove low res linear wind field
  ! ------------------------------------------------
  type, extends(linearizable_type) :: boundary_t
      type(meta_data_t)    :: info
      type(grid_t)         :: grid
      type(grid_t)         :: grid2d

      ! store the full high-res 3D grid for the next time step to compute dXdt fields
      ! includes high res versions of low res terrain and z
      ! type(domain_t) :: next_domain
      ! temporarily use old style domain type so it doesn't store coarrays... have to think about this
      type(domain_type) :: next_domain
      type(domain_type) :: current_domain


      ! store the timestep to the next input
      type(time_delta_t) :: dt

      type(Time_type) :: forcing_time

      ! dX_dt variables are the change in variable X between two forcing time steps
      ! wind and pressure dX_dt fields applied to full 3d grid, others applied only to boundaries
      real, allocatable, dimension(:,:,:) :: du_dt, dv_dt, dp_dt, dth_dt, dqv_dt, dqc_dt

      real, allocatable, dimension(:,:)   :: drain_dt
      ! change in shortwave and longwave at surface if read from forcing
      real, allocatable, dimension(:,:)   :: dsw_dt, dlw_dt
      ! change in sst if read from forcing file
      real, allocatable, dimension(:,:)   :: dsst_dt

      ! store the low resolution version of terrain and atmospheric elevations
      ! real,allocatable,dimension(:,:)     :: lowres_terrain
      real,allocatable,dimension(:,:,:)   :: lowres_z

      ! contains the size of the domain (or the local tile?)
      integer :: nx, ny, nz, nx_global, ny_global
      integer :: ximg, ximages, yimg, yimages

      ! store the start (s) and end (e) for the i,j,k dimensions
      integer ::  ids,ide, jds,jde, kds,kde, & ! for the entire model domain    (d)
                  ims,ime, jms,jme, kms,kme, & ! for the memory in these arrays (m)
                  its,ite, jts,jte, kts,kte    ! for the data tile to process   (t)

  contains
    procedure :: init

    procedure :: update_forcing

    procedure :: distribute_update
    procedure :: distribute_initial_conditions

  end type boundary_t

  interface

    ! Set default component values
    module subroutine init(this, options)
      implicit none
      class(boundary_t), intent(inout) :: this
      class(options_t),  intent(inout) :: options
    end subroutine

    module subroutine update_forcing(this, options)
        implicit none
        class(boundary_t), intent(inout) :: this
        class(options_t),  intent(inout) :: options
    end subroutine

    module subroutine distribute_update(this)
      implicit none
      class(boundary_t), intent(inout) :: this
    end subroutine

    module subroutine distribute_initial_conditions(this)
      implicit none
      class(boundary_t), intent(inout) :: this
    end subroutine

  end interface

end module
