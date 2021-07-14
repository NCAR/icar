module exchangeable_interface
  use assertions_mod,          only : assert, assertions
  use icar_constants,          only : kMAX_NAME_LENGTH
  use grid_interface,          only : grid_t
  use variable_interface,      only : variable_t
  implicit none

  private
  public :: exchangeable_t

  type :: exchangeable_t
    private
    ! real, pointer,     public :: local(:,:,:)   => null()
    real, pointer,     public :: data_3d(:,:,:) => null() ! provide the local data
    type(variable_t),  public :: meta_data

    real, allocatable :: halo_south_in(:,:,:)[:]
    real, allocatable :: halo_north_in(:,:,:)[:]
    real, allocatable :: halo_west_in(:,:,:)[:]
    real, allocatable :: halo_east_in(:,:,:)[:]

    logical :: north_boundary=.false.
    logical :: south_boundary=.false.
    logical :: east_boundary=.false.
    logical :: west_boundary=.false.

  contains
    private
    procedure, public :: const
    procedure, public :: set_neighbors
    procedure, public :: send
    procedure, public :: retrieve
    procedure, public :: exchange
    procedure, public :: exchange_u
    procedure, public :: exchange_v
    procedure, public :: set_outputdata
    generic,   public :: initialize=>const

    procedure :: put_north
    procedure :: put_south
    procedure :: put_west
    procedure :: put_east
    procedure :: retrieve_north_halo
    procedure :: retrieve_south_halo
    procedure :: retrieve_west_halo
    procedure :: retrieve_east_halo
  end type

  integer, parameter :: space_dim=3

  interface

    module subroutine const(this, grid, metadata, forcing_var)
      implicit none
      class(exchangeable_t),           intent(inout) :: this
      type(grid_t),                    intent(in)    :: grid
      type(variable_t),                intent(in),     optional :: metadata
      character(len=kMAX_NAME_LENGTH), intent(in),    optional :: forcing_var
    end subroutine

    module subroutine set_neighbors(this, grid)
        implicit none
        class(exchangeable_t), intent(inout) :: this
        type(grid_t),          intent(in)    :: grid
    end subroutine

    module subroutine set_outputdata(this, metadata)
      implicit none
      class(exchangeable_t), intent(inout)  :: this
      type(variable_t),      intent(in),    optional :: metadata
    end subroutine


    module subroutine send(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine retrieve(this, no_sync)
      implicit none
      class(exchangeable_t), intent(inout) :: this
      logical,               intent(in),   optional :: no_sync
    end subroutine

    module subroutine exchange(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine exchange_u(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine exchange_v(this)
      implicit none
      class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine put_north(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine put_south(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine retrieve_north_halo(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine retrieve_south_halo(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine


    module subroutine put_east(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine put_west(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine retrieve_east_halo(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine

    module subroutine retrieve_west_halo(this)
        implicit none
        class(exchangeable_t), intent(inout) :: this
    end subroutine


  end interface

end module
