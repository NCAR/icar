module grid_interface
  implicit none

  private
  public :: grid_t

  type grid_t
      integer :: yimg,    ximg
      integer :: yimages, ximages
      integer :: ims, ime
      integer :: jms, jme
      integer :: kms, kme
      integer :: ns_halo_nx, ew_halo_ny, halo_nz

  contains
      procedure :: get_dims
      
  end type

  interface
      module function get_dims(this) result(dims)
          implicit none
          class(grid_t) :: this
          integer :: dims(3)
      end function
  end interface
end module
