submodule(grid_interface) grid_implementation
    use assertions_mod, only : assert, assertions

    implicit none

contains

    !> -------------------------------
    !! Return the dimensions of this grid as an n element array
    !!
    !! -------------------------------
    module function get_dims(this) result(dims)
        class(grid_t), intent(in) :: this
        integer, allocatable :: dims(:)

        if (this%is2d) then
            allocate(dims(2))
            dims(1) = this%ime - this%ims + 1
            dims(2) = this%jme - this%jms + 1
        endif
        if (this%is3d) then
            allocate(dims(3))
            dims(1) = this%ime - this%ims + 1
            dims(2) = this%kme - this%kms + 1
            dims(3) = this%jme - this%jms + 1
        endif
    end function

    !> -------------------------------
    !! Decompose the domain into as even a set of tiles as possible in two dimensions
    !!
    !! Searches through possible numbers of x and y tiles that multiple evenly to
    !! give the total number of images requested.
    !!
    !! For each x/y split compute the number of grid cells in both dimensions in each tile
    !! return the split that provides the closest match between the number of x and y grid cells
    !!
    !! -------------------------------
    module subroutine domain_decomposition(this, nx, ny, nimages, ratio, for_image)
        class(grid_t),  intent(inout) :: this
        integer,        intent(in)    :: nx, ny, nimages
        real,           intent(in), optional :: ratio
        integer,        intent(in), optional :: for_image
        real    :: multiplier
        integer :: ysplit, xsplit, xs, ys, i
        real    :: best, current, x, y
        integer :: image

        multiplier=1
        if (present(ratio)) multiplier = ratio

        image = this_image()
        if (present(for_image)) image = for_image

        xsplit = 1
        ysplit = nimages
        xs = xsplit
        ys = ysplit

        x = (nx/real(xsplit))
        y = (ny/real(ysplit))

        if (y > (multiplier*x)) then
            best = abs(1 - ( y / (multiplier*x) ))
        else
            best = abs(1 - ( (multiplier*x) / y ))
        endif

        do i=nimages,1,-1
            if (mod(nimages,i)==0) then
                ysplit = i
                xsplit = nimages / i

                x = (nx/float(xsplit))
                y = (ny/float(ysplit))

                if (y > (multiplier*x)) then
                    current = abs(1 - ( y / (multiplier*x) ))
                else
                    current = abs(1 - ( (multiplier*x) / y ))
                endif

                if (current < best) then
                    best = current
                    xs = xsplit
                    ys = ysplit
                endif
            endif
        enddo

        this%ximages = xs
        this%yimages = ys

        this%ximg = mod(image-1,  this%ximages)+1
        this%yimg = floor(real(image-1) / this%ximages)+1

        x = (nx/float(xs))
        y = (ny/float(ys))

        if (assertions) call assert((xs*ys) == nimages, "Number of tiles does not sum to number of images")
        ! if (image==1) print*, "ximgs=",xs, "yimgs=",ys

    end subroutine domain_decomposition

    !> -------------------------------
    !! Compute the number of grid cells in the current image along a dimension
    !!
    !! This takes care of the fact that generally the number of images will not evenly divide the number of grid cells
    !! In this case the extra grid cells need to be evenly distributed among all images
    !!
    !! n_global should be the full domain size of the dimension
    !! me should be this images image number along this dimension
    !! nimg should be the number of images this dimension will be divided into
    !!
    !! -------------------------------
    function my_n(n_global, me, nimg) result(n_local)
       integer, intent(in) :: n_global, me, nimg
       integer :: n_local

       ! add 1 if this image is less than the remainder that need an extra grid cell
       n_local = n_global / nimg + merge(1,0,me <= mod(n_global,nimg)  )
    end function

    !> -------------------------------
    !! Return the starting coordinate in the global domain coordinate system for a given image (me)
    !!
    !! -------------------------------
    function my_start(n_global, me, nimg) result(memory_start)
        implicit none
        integer, intent(in) :: n_global, me, nimg
        integer :: memory_start
        integer :: base_n

        base_n = n_global / nimg

        memory_start = (me-1)*(base_n) + min(me-1,mod(n_global,nimg)) + 1

    end function my_start

    !> -------------------------------
    !! Generate the domain decomposition mapping and compute the indicies for local memory
    !!
    !! -------------------------------
    module subroutine set_grid_dimensions(this, nx, ny, nz, nx_extra, ny_extra, halo_width, for_image)
      class(grid_t),   intent(inout) :: this
      integer,         intent(in)    :: nx, ny, nz
      integer,         intent(in), optional :: nx_extra, ny_extra, halo_width, for_image

      integer :: nx_e, ny_e, halo_size
      integer :: image

      halo_size = kDEFAULT_HALO_SIZE
      if (present(halo_width)) halo_size = halo_width

      image = this_image()
      if (present(for_image)) image = for_image

      nx_e = 0
      ny_e = 0
      if (present(nx_extra)) nx_e = nx_extra ! used to add 1 to the u-field staggered grid
      if (present(ny_extra)) ny_e = ny_extra ! used to add 1 to the v-field staggered grid

      call this%domain_decomposition(nx, ny, num_images(), for_image=image)

      if (nz<1) then
          this%is2d = .True.
          this%is3d = .False.
          if (allocated(this%dimensions)) deallocate(this%dimensions)
          allocate(this%dimensions(2))
          this%dimensions(1) = "lat"
          this%dimensions(2) = "lon"
      else
          this%is2d = .False.
          this%is3d = .True.
          if (allocated(this%dimensions)) deallocate(this%dimensions)
          allocate(this%dimensions(3))
          this%dimensions(1) = "lat"
          this%dimensions(2) = "height"
          this%dimensions(3) = "lon"
      endif

      this%ny_global  = ny + ny_e                                     ! global model domain grid size
      this%nx_global  = nx + nx_e                                     ! global model domain grid size
      this%nz         = nz                                            ! note nz is both global and local
      this%nx         = my_n(this%nx_global-nx_e, this%ximg, this%ximages) ! local grid size
      this%ny         = my_n(this%ny_global-ny_e, this%yimg, this%yimages) ! local grid size

      ! define the bounds needed for memory to store the data local to this image
      this%ims        = my_start(this%nx_global-nx_e, this%ximg, this%ximages)
      this%ime        = this%ims + this%nx + nx_e - 1

      this%jms        = my_start(this%ny_global-ny_e, this%yimg, this%yimages)
      this%jme        = this%jms + this%ny + ny_e - 1

      this%kms        = 1
      this%kme        = this%nz

      ! Now define the tile of data to process in physics routines
      this%kts = this%kms
      this%kte = this%kme

      ! The entire model domain begins at 1 and ends at nx,y,z
      this%ids = 1
      this%jds = 1
      this%kds = 1
      this%ide = this%nx_global
      this%jde = this%ny_global
      this%kde = this%nz

      this%halo_nz    = this%nz

      this%halo_size = halo_size
      call update_with_halos(this, halo_size)

      ! define the halo needed to manage communications between images
      ! perhaps this should be defined in exchangeable instead though?
      this%ns_halo_nx = this%nx !_global / this%ximages + 1 + nx_e  ! number of grid cells in x in the ns halo
      this%ew_halo_ny = this%ny !_global / this%yimages + 1 + ny_e  ! number of grid cells in y in the ew halo


  end subroutine

  !> -------------------------------
  !! updates the grid memory dimensions with halo sizes if necessary
  !!
  !! -------------------------------
  subroutine update_with_halos(grid, halo_size)
      type(grid_t), intent(inout)   :: grid
      integer,      intent(in)      :: halo_size

      logical :: north_boundary, south_boundary, &
                 east_boundary,  west_boundary

      north_boundary = (grid%yimg == grid%yimages)
      south_boundary = (grid%yimg == 1)
      east_boundary  = (grid%ximg == grid%ximages)
      west_boundary  = (grid%ximg == 1)

      ! if this is on a given boundary, then add 0, if it is not a boundary than add/subtract halo_size
      grid%ims = grid%ims - merge(0, halo_size, west_boundary)
      grid%ime = grid%ime + merge(0, halo_size, east_boundary)
      grid%jms = grid%jms - merge(0, halo_size, south_boundary)
      grid%jme = grid%jme + merge(0, halo_size, north_boundary)

      ! if this is on a boundary, we should skip 1 grid cell (the boundary conditions) else we should skip the halo
      grid%its = grid%ims + merge(1, halo_size, west_boundary)
      grid%ite = grid%ime - merge(1, halo_size, east_boundary)
      grid%jts = grid%jms + merge(1, halo_size, south_boundary)
      grid%jte = grid%jme - merge(1, halo_size, north_boundary)

      grid%nx = grid%ime - grid%ims + 1
      grid%ny = grid%jme - grid%jms + 1

  end subroutine

end submodule
