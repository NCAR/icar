!>------------------------------------------------------------
!!  Handles reading boundary conditions from the forcing file(s)
!!  Provides necessary interpolation on to the grid. 
!!
!! <pre>
!!  Primary entry points 
!!      bc_init      - first call only
!!      bc_update    - all successive calls
!!
!!  Contains options to use the 
!!      mean wind field
!!      mean boundary forcing
!!      wind field smoothing
!!      removal of low resolution model linear wind field
!!
!!  bc_init loads the restart file as necessary. 
!!
!!  Both init and update compute the exner function and density fields
!!  for the forcing step and update the wind field with linear perturbations
!!  as necessary
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module boundary_conditions
! ----------------------------------------------------------------------------
!   NOTE: This module attempts to be sufficiently general to work with a variety
!       of possible input files; however it may be necessary to modify it.
! ----------------------------------------------------------------------------
    use data_structures
    use io_routines,            only : io_getdims, io_read3d, io_maxDims, io_read2d, io_variable_is_present, &
                                       io_write
    use wind,                   only : update_winds,balance_uvw
    use linear_theory_winds,    only : linear_perturb
    use geo,                    only : geo_interp2d, geo_interp
    use vertical_interpolation, only : vinterp, vLUT_forcing
    use output,                 only : write_domain, output_init
    use string,                 only : str
    
    implicit none
    private
!   these could be better stored in bc_type and initialized in init_bc?
    character (len=255),dimension(:),allocatable :: file_list
!   manage file pointer and position in file for boundary conditions
    integer::curfile,curstep
    integer::steps_in_file,nfiles
!   manage file pointer and position in file for external winds
    character (len=255),dimension(:),allocatable :: ext_winds_file_list
    integer::ext_winds_curfile,ext_winds_curstep
    integer::ext_winds_steps_in_file,ext_winds_nfiles

    integer::smoothing_window=1 ! this will get updated in bc_init if it isn't an ideal run
    
    public :: bc_init
    public :: bc_update
    public :: bc_find_step
    public :: update_pressure
contains

    !>------------------------------------------------------------
    !! Find the time step in the input forcing to start the model on 
    !!
    !! The model start date (start_mjd) may not be the same as the first forcing
    !! date (initial_mjd).  Convert the difference between the two into forcing
    !! steps by dividing by the time delta between forcing steps (in_dt) after 
    !! converting in_dt from seconds to days. 
    !!
    !! @param  options  model options structure
    !! @retval step     integer number of steps into the forcing sequence
    !!
    !!------------------------------------------------------------
    function bc_find_step(options) result(step)
        implicit none
        type(options_type), intent(in) :: options
        integer :: step
        
        step = (options%start_mjd-options%initial_mjd)/(options%in_dt / 86400.0d+0)
        if (options%debug) write(*,*), "bc_find_step: First forcing time step = ",trim(str(step))
        
    end function bc_find_step
    
    !>------------------------------------------------------------
    !! Smooth an array (written for wind but will work for anything)
    !!
    !! Only smooths over the first (x) and second (y or z) or third (y or z) dimension
    !! ydim can be specified to allow working with (x,y,z) data or (x,z,y) data
    !! WARNING: this is a moderately complex setup to be efficient for the ydim=3 (typically large arrays, SLOW) case
    !! be careful when editing.  
    !! For the complex case it pre-computes the sum of all columns for a given row, 
    !! then to move from one column to the next it just has add the next column from the sums and subtracts the last one
    !! similarly, moving to the next row just means adding the next row to the sums, and subtracting the last one. 
    !! Each point also has to be divided by N, but this decreases the compution from O(windowsize^2) to O(constant) 
    !! Where O(constant) = 2 additions, 2 subtractions, and 1 divide regardless of windowsize!
    !!
    !! @param wind          3D array to be smoothed
    !! @param windowsize    size to smooth in both directions (i.e. the full window is this * 2 + 1)
    !! @param ysim          axis the y dimension is on (2 or 3) in the wind array
    !!
    !!------------------------------------------------------------
    subroutine smooth_wind(wind,windowsize,ydim)
        implicit none
        real, intent(inout), dimension(:,:,:):: wind    !> 3 dimensional wind field to be smoothed
        integer,intent(in)::windowsize                  !> halfwidth-1/2 of window to smooth over
                                                        ! Specified in grid cells, (+/- windowsize)
        integer,intent(in)::ydim                        !> the dimension to use for the y coordinate
                                                        ! It can be 2, or 3 (but not 1)
        real,allocatable,dimension(:,:,:)::inputwind    !> temporary array to store the input data in
        integer::i,j,k,nx,ny,nz,startx,endx,starty,endy ! various array indices/bounds
        ! intermediate sums to speed up the computation
        real,allocatable,dimension(:) :: rowsums,rowmeans
        real :: cursum
        integer :: cur_n,curcol,ncols,nrows
        
        ncols=windowsize*2+1
        nx=size(wind,1)
        ny=size(wind,2) !note, this is Z for the high-res domain (ydim=3)
        nz=size(wind,3) !note, this could be the Y or Z dimension depending on ydim
        ! this is ~20MB that has to be allocated deallocated every call...
        allocate(inputwind(nx,ny,nz)) ! Can't be module level because nx,ny,nz could change between calls, 
                                      ! could be part of a "smoothable" object to avoid allocate-deallocating constantly
!       nx=nx-1
!       if (ydim==3) then
!           nz=nz-1
!       endif
        inputwind=wind !make a copy so we always use the unsmoothed data when computing the smoothed data
        if (((windowsize*2+1)>nx).and.(ydim==3)) then
            write(*,*) "WARNING can not operate if windowsize*2+1 is larger than nx"
            write(*,*) "NX         = ", nx
            write(*,*) "windowsize = ", windowsize
            stop
        endif
        
        !parallelize over a slower dimension (not the slowest because it is MUCH easier this way)
        ! as long as the inner loops (the array operations) are over the fastest dimension we are mostly OK
        !$omp parallel firstprivate(windowsize,nx,ny,nz,ydim), &
        !$omp private(i,j,k,startx,endx,starty,endy, rowsums,rowmeans,nrows,ncols,cursum), &
        !$omp shared(wind,inputwind)
        allocate(rowsums(nx)) !this is only used when ydim=3, so nz is really ny
        allocate(rowmeans(nx)) !this is only used when ydim=3, so nz is really ny
        nrows=windowsize*2+1
        ncols=windowsize*2+1
        !$omp do schedule(static)
        do j=1,ny
            
            ! so we pre-compute the sum over rows for each column in the current window
            if (ydim==3) then
                rowsums=inputwind(1:nx,j,1)*(windowsize+1)
                do i=1,windowsize
                    rowsums=rowsums+inputwind(1:nx,j,i)
                enddo
            endif
            ! don't parallelize over this loop because it is much more efficient to be able to assume
            ! that you ran the previous row in serial for the slow (ydim=3) case
            do k=1,nz ! note this is y for ydim=3
                ! ydim=3 for the main model grid which is large and takes a long time
                if (ydim==3) then
                    ! if we are pinned to the top edge
                    if ((k-windowsize)<=1) then
                        starty=1
                        endy  =k+windowsize
                        rowsums=rowsums-inputwind(1:nx,j,starty)+inputwind(1:nx,j,endy)
                    ! if we are pinned to the bottom edge
                    else if ((k+windowsize)>nz) then
                        starty=k-windowsize
                        endy  =nz
                        rowsums=rowsums-inputwind(1:nx,j,starty-1)+inputwind(1:nx,j,endy)
                    ! if we are in the middle (this is the most common)
                    else
                        starty=k-windowsize
                        endy  =k+windowsize
                        rowsums=rowsums-inputwind(1:nx,j,starty-1)+inputwind(:,j,endy)
                    endif
                    rowmeans=rowsums/nrows
                    cursum=sum(rowmeans(1:windowsize))+rowmeans(1)*(windowsize+1)
                endif
                
                do i=1,nx
                    if (ydim==3) then
                        ! if we are pinned to the left edge
                        if ((i-windowsize)<=1) then
                            startx=1
                            endx  =i+windowsize
                            cursum=cursum-rowmeans(startx)+rowmeans(endx)
                        ! if we are pinned to the right edge
                        else if ((i+windowsize)>nx) then
                            startx=i-windowsize
                            endx  =nx
                            cursum=cursum-rowmeans(startx-1)+rowmeans(endx)
                        ! if we are in the middle (this is the most common)
                        else
                            startx=i-windowsize
                            endx  =i+windowsize
                            cursum=cursum-rowmeans(startx-1)+rowmeans(endx)
                        endif

                        wind(i,j,k)=cursum/ncols

                        ! old slower way, also creates artifacts windowsize distance from the borders
!                       startx=max(1, i-windowsize)
!                       endx  =min(nx,i+windowsize)
!                       starty=max(1, k-windowsize)
!                       endy  =min(nz,k+windowsize)
!                       wind(i,j,k)=sum(inputwind(startx:endx,j,starty:endy)) &
!                                   / ((endx-startx+1)*(endy-starty+1))
                    else ! ydim==2
                        ! ydim=2 for the input data which is a small grid, thus cheap so we still use the slow method
                        !first find the current window bounds
                        startx=max(1, i-windowsize)
                        endx  =min(nx,i+windowsize)
                        starty=max(1, j-windowsize)
                        endy  =min(ny,j+windowsize)
                        ! then compute the mean within that window (sum/n)
                        ! note, artifacts near the borders (mentioned in ydim==3) don't affect this
                        ! because the borders *should* be well away from the domain
                        ! if the forcing data are not much larger than the model domain this *could* create issues
                        wind(i,j,k)=sum(inputwind(startx:endx,starty:endy,k)) &
                                    / ((endx-startx+1)*(endy-starty+1))
                    endif
                enddo
            enddo
        enddo
        !$omp end do
        deallocate(rowmeans,rowsums)
        !$omp end parallel
        
        deallocate(inputwind)
    end subroutine smooth_wind
    
    !>------------------------------------------------------------
    !!  Generic routine to read a forcing variable (varname) from a netcdf file (filename) at a time step (curstep)
    !!  
    !!  Data are interpolated to the high res grid either in 3D or at the boundaries only (boundary_only)
    !!  Applies modifications specificaly for U,V,T, and P variables
    !!
    !! @param highres       Allocated input array to store the result.
    !! @param filename      Name of the NetCDF file to read.
    !! @param varname       Name of the variable to read from <filename>.
    !! @param geolut        Geographic Look Up Table that defines the 2D horizontal interpolation.
    !! @param vlut          Vertical Look Up table that defines the vertical interpolation.
    !! @param curstep       The time step in <filename> to read. 
    !! @param boundary_only Logical:interpolate to all points or just to boundaries
    !! @param options       Model options structure
    !! @param z_lo          Low-resolution (input) 3D vertical coordinate for pressure adjustment. [meters]
    !! @param z_hi          High-resolution (output) 3D vertical coordinate for pressure adjustment. [meters]
    !! @param time_varying_zlut Logical: if true, first adjust the low-resolution data to a common vertical coordinate
    !! @param interp_vertical Logical: if true, perform vertical interpolation (false for pressure)
    !!
    !!------------------------------------------------------------
    subroutine read_var(highres, filename, varname, geolut, vlut, curstep, boundary_only, options, &
                        z_lo, z_hi, time_varying_zlut, interp_vertical)
        implicit none
        real,dimension(:,:,:),   intent(inout):: highres
        character(len=*),        intent(in)   :: filename,varname
        type(geo_look_up_table), intent(in)   :: geolut
        type(vert_look_up_table),intent(in)   :: vlut
        integer,                 intent(in)   :: curstep
        logical,                 intent(in)   :: boundary_only
        type(options_type),      intent(in)   :: options
        real,dimension(:,:,:),   intent(in), optional :: z_lo, z_hi
        type(vert_look_up_table),intent(in), optional :: time_varying_zlut
        logical,                 intent(in), optional :: interp_vertical
        
        ! local variables
        real,dimension(:,:,:),allocatable :: inputdata,extra_data
        integer :: nx,ny,nz,i
        logical :: apply_vertical_interpolation
        
        ! the special case is to skip vertical interpolation (pressure and one pass of temperature)
        if (present(interp_vertical)) then
            apply_vertical_interpolation=interp_vertical
        else
            ! so default to applying interpolation
            apply_vertical_interpolation=.True.
        endif
        
        ! Read the data in, should be relatively fast because we are reading a low resolution forcing file
        call io_read3d(filename,varname,inputdata,curstep)
        
        nx=size(inputdata,1)
        ny=size(inputdata,2)
        nz=size(inputdata,3)
        
        ! Variable specific options
        ! For wind variables run a first pass of smoothing over the low res data
        if (((varname==options%vvar).or.(varname==options%uvar)).and.(.not.options%ideal)) then
            call smooth_wind(inputdata,1,2)
            
        ! For Temperature, we may need to add an offset
        else if ((varname==options%tvar).and.(options%t_offset/=0)) then
            inputdata=inputdata+options%t_offset
        
        ! For pressure, we may need to add a base pressure offset read from pbvar
        else if ((varname==options%pvar).and.(options%pbvar/='')) then
            call io_read3d(filename,options%pbvar,extra_data,curstep)
            inputdata=inputdata+extra_data
            deallocate(extra_data)
        endif
        
        ! if the z-axis of the input data varies over time, we need to first interpolate
        ! to the "standard" z-axis so that the hi-res vlut doesn't need to change
        if ((options%time_varying_z).and.present(time_varying_zlut)) then
            allocate(extra_data(nx,ny,nz))
            extra_data=inputdata
            call vinterp(inputdata, extra_data, time_varying_zlut, axis=3)
            deallocate(extra_data)
        endif
            
        
        ! just read the low res version with out interpolating for e.g. external wind data
        if ( (nx==size(highres,1)).and.(ny==size(highres,3)).and.(nz==size(highres,2)) ) then
            highres=reshape(inputdata,[nx,nz,ny],order=[1,3,2])
            deallocate(inputdata)
        else
            ! interpolate data onto the high resolution grid after re-arranging the dimensions. 
            allocate(extra_data(nx,nz,ny))
            extra_data=reshape(inputdata,[nx,nz,ny],order=[1,3,2])
        
            ! first interpolate to a high res grid (temporarily stored in inputdata)
            deallocate(inputdata)
            allocate(inputdata(size(highres,1),nz,size(highres,3)))
            call geo_interp(inputdata, &
                            extra_data, &
                            geolut,boundary_only)
            ! Then apply vertical interpolation on that grid
            if (apply_vertical_interpolation) then
                call vinterp(highres, inputdata, &
                             vlut,boundary_only)
            else
                ! if we aren't interpolating, just copy over the output
                highres=inputdata(:,:size(highres,2),:)
            endif
            deallocate(extra_data)
            deallocate(inputdata)
        endif
        
        ! highres is the useful output of the subroutine
    end subroutine read_var

    !>------------------------------------------------------------
    !!  Same as read_var but for 2-dimensional data
    !!
    !!  Data are read and horizontally interpolated. This routine is simpler because it is used For
    !!  a more limited set of variables, and they are simpler 2D instead of 3D. 
    !!  Primarily used for surface variables: Sensible and latent heat fluxes, PBL height, skin temperature, radiation
    !!
    !! @param highres       Allocated input array to store the result.
    !! @param filename      Name of the NetCDF file to read.
    !! @param varname       Name of the variable to read from <filename>.
    !! @param geolut        Geographic Look Up Table that defines the 2D horizontal interpolation.
    !! @param curstep       The time step in <filename> to read. 
    !!
    !!------------------------------------------------------------
    subroutine read_2dvar(highres,filename,varname,geolut,curstep)
        implicit none
        real,dimension(:,:),intent(inout)::highres
        character(len=*),intent(in) :: filename,varname
        type(geo_look_up_table),intent(in) :: geolut
        integer,intent(in)::curstep
    
        real,dimension(:,:),allocatable :: inputdata
    
!       Read the data in
        call io_read2d(filename,varname,inputdata,curstep)
!       interpolate data onto the high resolution grid
        call geo_interp2d(highres,inputdata,geolut)
        deallocate(inputdata)
                    
    end subroutine read_2dvar

    
    !>------------------------------------------------------------
    !!  Rotate winds from real space back to terrain following grid (approximately)
    !!
    !!  Only used for external winds. 
    !!  Assumes a simple slope transform in u and v independantly constant w/height
    !!
    !! @param domain    Full model domain
    !! @param ext_winds External wind data structure 
    !!
    !!------------------------------------------------------------
    subroutine rotate_ext_wind_field(domain,ext_winds)
        implicit none
        type(domain_type),intent(inout)::domain
        type(wind_type),intent(inout)::ext_winds
        integer :: nx,ny,nz,i
    
        nx=size(domain%z,1)
        nz=size(domain%z,2)
        ny=size(domain%z,3)
        do i=1,nz
            domain%u(1:nx-1,i,:)=domain%u(1:nx-1,i,:)*ext_winds%dzdx
            domain%v(:,i,1:ny-1)=domain%v(:,i,1:ny-1)*ext_winds%dzdy
        end do
    
    end subroutine rotate_ext_wind_field

    
    !>------------------------------------------------------------
    !!  Initialize the external winds information (filenames, nfiles, etc) and read the initial conditions
    !!
    !! @param domain    Model Domain
    !! @param bc        Model forcing boundary conditions data structure
    !! @param options   Model Options structure 
    !!
    !!------------------------------------------------------------
    subroutine ext_winds_init(domain,bc,options)
        implicit none
        type(domain_type),intent(inout)::domain
        type(bc_type),intent(inout)::bc
        type(options_type),intent(in)::options
        integer,dimension(io_maxDims)::dims !note, io_maxDims is included from io_routines.
        ! MODULE variables : ext_winds_ curstep, curfile, nfiles, steps_in_file, file_list
        ext_winds_curfile=1
        if (options%restart) then
            ext_winds_curstep=options%restart_step
        else
            ext_winds_curstep=1
        endif
        ext_winds_nfiles=options%ext_winds_nfiles
        allocate(ext_winds_file_list(ext_winds_nfiles))
        ext_winds_file_list=options%ext_wind_files
        call io_getdims(ext_winds_file_list(ext_winds_curfile),options%uvar, dims)
        if (dims(1)==3) then
            ext_winds_steps_in_file=1
        else
            ext_winds_steps_in_file=dims(dims(1)+1) !dims(1) = ndims
        endif
        
        do while (ext_winds_curstep>ext_winds_steps_in_file)
            ext_winds_curfile=ext_winds_curfile+1
            if (ext_winds_curfile>ext_winds_nfiles) then
                stop "Ran out of files to process!"
            endif
            ext_winds_curstep=ext_winds_curstep-ext_winds_steps_in_file 
            !instead of setting=1, this way we can set an arbitrary starting point multiple files in
            call io_getdims(ext_winds_file_list(ext_winds_curfile),options%uvar, dims)
            if (dims(1)==3) then
                ext_winds_steps_in_file=1
            else
                ext_winds_steps_in_file=dims(dims(1)+1) !dims(1) = ndims; dims(ndims+1)=ntimesteps
            endif
        enddo
        
        write(*,*) "Initial external wind file= ",ext_winds_curfile," : step= ",ext_winds_curstep
        call read_var(domain%u, ext_winds_file_list(ext_winds_curfile), options%uvar,  &
              bc%ext_winds%u_geo%geolut, bc%u_geo%vert_lut, ext_winds_curstep, .FALSE., options)
              
        call read_var(domain%v,    ext_winds_file_list(ext_winds_curfile),  options%vvar,  &
              bc%ext_winds%v_geo%geolut,bc%v_geo%vert_lut,ext_winds_curstep,.FALSE.,options)
        call rotate_ext_wind_field(domain,bc%ext_winds)
    end subroutine ext_winds_init
    
    !>------------------------------------------------------------
    !!  Remove linear theory topographic winds perturbations from the low resolution wind field. 
    !!
    !! NOTE: this may not work at the moment (and certainly won't with the Linear wind look up table option)
    !!
    !! @param domain    Model Domain
    !! @param bc        Model forcing boundary conditions data structure
    !! @param options   Model Options structure 
    !! @param filename  Name of the file to read the current boundary conditions forcing from (as in read_var)
    !! @param curstep   Time step in forcing file to read 
    !!
    !!------------------------------------------------------------
    subroutine remove_linear_winds(domain,bc,options,filename,curstep)
        implicit none
        type(domain_type), intent(inout) :: domain
        type(bc_type), intent(inout) :: bc
        type(options_type),intent(in) :: options
        character(len=*),intent(in) :: filename
        integer,intent(in)::curstep
        real,allocatable,dimension(:,:,:)::inputdata,extra_data
        logical :: reverse_winds=.TRUE.
        integer :: nx,ny,nz
        character(len=255) :: outputfilename
        
        ! first read in the low-res U and V data directly
        ! load low-res U data
        call io_read3d(filename,options%uvar,extra_data,curstep)
        nx=size(extra_data,1)
        ny=size(extra_data,2)
        nz=size(extra_data,3)
        call smooth_wind(extra_data,1,2)
        bc%u=reshape(extra_data,[nx,nz,ny],order=[1,3,2])
        deallocate(extra_data)
        
        ! load low-res V data
        call io_read3d(filename,options%vvar,extra_data,curstep)
        nx=size(extra_data,1)
        ny=size(extra_data,2)
        nz=size(extra_data,3)
        call smooth_wind(extra_data,1,2)
        bc%v=reshape(extra_data,[nx,nz,ny],order=[1,3,2])
        deallocate(extra_data)
        
        ! remove the low-res linear wind contribution effect
        call linear_perturb(bc,options,options%lt_options%vert_smooth,reverse_winds,options%advect_density)
        
        ! finally interpolate low res winds to the high resolutions grid
        nx=size(domain%u,1)
        nz=size(bc%u,2)
        ny=size(domain%u,3)
        allocate(extra_data(nx,nz,ny))
        call geo_interp(extra_data, bc%u,bc%u_geo%geolut,.FALSE.)
        call vinterp(domain%u,extra_data,bc%u_geo%vert_lut)
        deallocate(extra_data)
        
        nx=size(domain%v,1)
        nz=size(bc%v,2)
        ny=size(domain%v,3)
        allocate(extra_data(nx,nz,ny))
        call geo_interp(extra_data, bc%v,bc%v_geo%geolut,.FALSE.)
        call vinterp(domain%v,extra_data,bc%v_geo%vert_lut)
        deallocate(extra_data)
    end subroutine remove_linear_winds
    
    !>------------------------------------------------------------
    !! Make the forcing boundary conditions into average values
    !!
    !! Averages are computed independantly for each level and each side of the model
    !!
    !! @param inputdata     3D Array to compute means over the boundaries
    !!
    !!------------------------------------------------------------
    subroutine mean_boundaries(inputdata)
        implicit none
        real, dimension(:,:,:), intent(inout) :: inputdata
        integer:: nx, ny, nz, i

        nx=size(inputdata,1)
        nz=size(inputdata,2)
        ny=size(inputdata,3)
        do i=1,nz
            inputdata(1,i,:)  = sum(inputdata(1,i,:))  / ny
            inputdata(nx,i,:) = sum(inputdata(nx,i,:)) / ny
            inputdata(:,i,1)  = sum(inputdata(:,i,1))  / nx
            inputdata(:,i,ny) = sum(inputdata(:,i,ny)) / nx
        end do
        
    end subroutine mean_boundaries

    !>------------------------------------------------------------
    !!  Compute the mean wind field and make them constant everywhere...
    !!
    !!  Averages are computed over the entire domain. U and V are read in hear to avoid the 
    !!  cost of first interpolating them to the high-res grid. BUT this means that a much larger domain
    !!  then the one ICAR is actually running on is used, so it should be changed to use the high-res grid...
    !!
    !! @param domain    Model Domain
    !! @param filename  Name of the file to read the current boundary conditions forcing from (as in read_var)
    !! @param curstep   Time step in forcing file to read 
    !! @param options   Model Options structure 
    !!
    !!------------------------------------------------------------
    subroutine mean_winds(domain,filename,curstep,options)
        implicit none
        type(domain_type), intent(inout) :: domain
        character(len=*),intent(in)::filename
        integer,intent(in)::curstep
        type(options_type):: options
        
        real,allocatable,dimension(:,:,:)::extra_data
        integer::nx,ny,nz
        
        nz=size(domain%u,2)
        
!       load low-res U data
        call io_read3d(filename,options%uvar,extra_data,curstep)
        domain%u=sum(extra_data(:,:,:nz))/size(extra_data(:,:,:nz))
        deallocate(extra_data)

!       load low-res V data
        call io_read3d(filename,options%vvar,extra_data,curstep)
        domain%v=sum(extra_data(:,:,:nz))/size(extra_data(:,:,:nz))
        deallocate(extra_data)
                
    end subroutine mean_winds
    
    !>------------------------------------------------------------
    !! Check that two model grids have the same shape
    !!
    !! If the size of all dimensions in data1 and data2 are not exactly the same the model will stop
    !!
    !! @param data1     First 3D array to check dimension sizes
    !! @param data2     Second 3D array to check dimension sizes
    !!
    !!------------------------------------------------------------
    subroutine check_shapes_3d(data1,data2)
        implicit none
        real,dimension(:,:,:),intent(in)::data1,data2
        integer :: i
        do i=1,3
            if (size(data1,i).ne.size(data2,i)) then
                write(*,*) "Restart file 3D dimensions don't match domain"
                write(*,*) shape(data1)
                write(*,*) shape(data2)
                stop
            endif
        enddo
    end subroutine check_shapes_3d
    
    !>------------------------------------------------------------
    !! Swap the last two dimensions of an array
    !!
    !! Call reshape after finding the nx,ny,nz values
    !!
    !! @param data     3D array to be reshaped
    !!
    !!------------------------------------------------------------
    subroutine swap_y_z_dimensions(data)
        implicit none
        real,dimension(:,:,:),intent(inout),allocatable :: data
        real,dimension(:,:,:), allocatable :: temporary_data
        integer :: nx,ny,nz
        
        nx=size(data,1)
        ny=size(data,2)
        nz=size(data,3)
        allocate(temporary_data(nx,nz,ny))
        temporary_data = reshape(data, [nx,nz,ny], order=[1,3,2])

        deallocate(data)
        allocate(data(nx,nz,ny))
        data=temporary_data
        
    end subroutine swap_y_z_dimensions
    
    !>------------------------------------------------------------
    !!  Load restart file
    !!
    !!  If we are restarting from a given point, initialize the domain from the given restart file
    !!
    !! @param domain        Model domain
    !! @param restart_file  Name of the file to read restart data from
    !! @param time_step     Time step in restart_file to read. 
    !!
    !!------------------------------------------------------------
    subroutine load_restart_file(domain,restart_file,time_step)
        implicit none
        type(domain_type), intent(inout) :: domain
        character(len=*),intent(in)::restart_file
        integer,optional,intent(in) :: time_step
        real,allocatable,dimension(:,:,:)::inputdata
        real,allocatable,dimension(:,:)::inputdata_2d
        integer :: timeslice
        
        if (present(time_step)) then
            timeslice=time_step
        else
            timeslice=1
        endif
        
        write(*,*) "Reading atmospheric restart data"
        write(*,*) "   timestep:",trim(str(timeslice))," from file:",trim(restart_file)
        call io_read3d(restart_file,"qv",inputdata,timeslice)
        call swap_y_z_dimensions(inputdata)
        call check_shapes_3d(inputdata,domain%qv)
        domain%qv=inputdata
        deallocate(inputdata)
        
        call io_read3d(restart_file,"qc",inputdata,timeslice)
        call swap_y_z_dimensions(inputdata)
        domain%cloud=inputdata
        deallocate(inputdata)

        call io_read3d(restart_file,"qr",inputdata,timeslice)
        call swap_y_z_dimensions(inputdata)
        domain%qrain=inputdata
        deallocate(inputdata)
        call io_read3d(restart_file,"qi",inputdata,timeslice)
        call swap_y_z_dimensions(inputdata)
        domain%ice=inputdata
        deallocate(inputdata)
        call io_read3d(restart_file,"qs",inputdata,timeslice)
        call swap_y_z_dimensions(inputdata)
        domain%qsnow=inputdata
        deallocate(inputdata)
        if (io_variable_is_present(restart_file,"qg")) then
            call io_read3d(restart_file,"qg",inputdata,timeslice)
            call swap_y_z_dimensions(inputdata)
            domain%qgrau=inputdata
            deallocate(inputdata)
        endif
        if (io_variable_is_present(restart_file,"nr")) then
            call io_read3d(restart_file,"nr",inputdata,timeslice)
            call swap_y_z_dimensions(inputdata)
            domain%nrain=inputdata
            deallocate(inputdata)
        endif
        if (io_variable_is_present(restart_file,"ni")) then
            call io_read3d(restart_file,"ni",inputdata,timeslice)
            call swap_y_z_dimensions(inputdata)
            domain%nice=inputdata
            deallocate(inputdata)
        endif
        call io_read3d(restart_file,"p",inputdata,timeslice)
        call swap_y_z_dimensions(inputdata)
        domain%p=inputdata
        deallocate(inputdata)
        call io_read3d(restart_file,"th",inputdata,timeslice)
        call swap_y_z_dimensions(inputdata)
        domain%th=inputdata
        deallocate(inputdata)
        call io_read3d(restart_file,"rho",inputdata,timeslice)
        call swap_y_z_dimensions(inputdata)
        domain%rho=inputdata
        deallocate(inputdata)
        
        call io_read2d(restart_file,"rain",inputdata_2d,timeslice)
        domain%rain=inputdata_2d
        deallocate(inputdata_2d)
        call io_read2d(restart_file,"snow",inputdata_2d,timeslice)
        domain%snow=inputdata_2d
        deallocate(inputdata_2d)
        if (io_variable_is_present(restart_file,"graupel")) then
            call io_read2d(restart_file,"graupel",inputdata_2d,timeslice)
            domain%graupel=inputdata_2d
            deallocate(inputdata_2d)
        endif
        if (io_variable_is_present(restart_file,"crain")) then
            call io_read2d(restart_file,"crain",inputdata_2d,timeslice)
            domain%crain=inputdata_2d
            deallocate(inputdata_2d)
        endif
        
        if (io_variable_is_present(restart_file,"soil_t")) then
            write(*,*) "Reading land surface restart data"
            call io_read3d(restart_file,"soil_t",inputdata,timeslice)
            call swap_y_z_dimensions(inputdata)
            call check_shapes_3d(inputdata,domain%soil_t)
            domain%soil_t=inputdata
            deallocate(inputdata)
            call io_read3d(restart_file,"soil_w",inputdata,timeslice)
            call swap_y_z_dimensions(inputdata)
            domain%soil_vwc=inputdata
            deallocate(inputdata)
        
            call io_read2d(restart_file,"ts",inputdata_2d,timeslice)
            domain%skin_t=inputdata_2d
            deallocate(inputdata_2d)
            call io_read2d(restart_file,"hfgs",inputdata_2d,timeslice)
            domain%ground_heat=inputdata_2d
            deallocate(inputdata_2d)
            call io_read2d(restart_file,"snw",inputdata_2d,timeslice)
            domain%snow_swe=inputdata_2d
            deallocate(inputdata_2d)
            call io_read2d(restart_file,"canwat",inputdata_2d,timeslice)
            domain%canopy_water=inputdata_2d
            deallocate(inputdata_2d)
        endif
        
        if (io_variable_is_present(restart_file,"rsds")) then
            call io_read2d(restart_file,"rsds",inputdata_2d,timeslice)
            domain%swdown=inputdata_2d
            deallocate(inputdata_2d)
            
            call io_read2d(restart_file,"rlds",inputdata_2d,timeslice)
            domain%lwdown=inputdata_2d
            deallocate(inputdata_2d)
        endif
        
        
    end subroutine load_restart_file
    
    !>------------------------------------------------------------
    !! Calculate the ZNU and ZNW variables
    !!
    !! @param domain    Model domain structure
    !!
    !!------------------------------------------------------------
    subroutine init_znu(domain)
        implicit none
        type(domain_type), intent(inout) :: domain
        integer :: n_levels,i,xpt,ypt
        real    :: ptop,psfc
        
        n_levels=size(domain%p,2)
        
        ! one grid point into the domain gets a non-boundary point
        xpt=2
        ypt=2
        ptop=domain%p(xpt,n_levels,ypt)-(domain%p(xpt,n_levels-1,ypt)-domain%p(xpt,n_levels,ypt))/2.0 !NOT CORRECT
        psfc=domain%p(xpt,1,ypt)+(domain%p(xpt,1,ypt)-domain%p(xpt,2,ypt))/2.0 !NOT CORRECT
        ptop=max(ptop,1.0)
        allocate(domain%znu(n_levels))
        allocate(domain%znw(n_levels))
        do i=1,n_levels
            domain%znu(i)=(domain%p(xpt,i,ypt)-ptop)/(psfc-ptop)
            if (i>1) then
                domain%znw(i)=((domain%p(xpt,i,ypt)+domain%p(xpt,i-1,ypt))/2-ptop)/(psfc-ptop)
            else
                domain%znw(i)=1
            endif
        enddo
    end subroutine init_znu

    
    !>------------------------------------------------------------
    !!  Initialize the boundary conditions (read inital conditions, etc.)
    !!
    !!  This is one of the most important subroutines, and is nearly identical to bc_update
    !!  This reads the first forcing time step in, calling relevant subroutines as dictated by
    !!  the options specified (e.g. mean_winds, external_winds)
    !!
    !! @param domain    Model Domain
    !! @param bc        Model forcing boundary conditions data structure
    !! @param options   Model Options structure 
    !!
    !!------------------------------------------------------------
    subroutine bc_init(domain,bc,options)
        implicit none
        type(domain_type),intent(inout)::domain
        type(bc_type),intent(inout)::bc
        type(options_type),intent(in)::options
        integer,dimension(io_maxDims)::dims !note, io_maxDims is included from io_routines.
        real,dimension(:,:,:),allocatable::inputdata
        logical :: boundary_value
        integer::nx,ny,nz,i
        real::domainsize
        ! MODULE variables : curstep, curfile, nfiles, steps_in_file, file_list
        
!       in case we are using a restart file we have some trickery to do here to find the proper file to be reading from
!       and set the current time step appropriately... should probably be moved to a subroutine. 
        curfile=1
        if (options%restart) then
            curstep=options%restart_step
        else
            curstep=bc_find_step(options)+1
        endif
        nfiles=size(options%boundary_files)
        allocate(file_list(nfiles))
        file_list=options%boundary_files
        call io_getdims(file_list(curfile),options%pvar, dims)
        if (dims(1)==3) then
            steps_in_file=1
        else
            steps_in_file=dims(dims(1)+1) !dims(1) = ndims
        endif
        
        if (.not.options%ideal) then
            do while (curstep>steps_in_file)
                curfile=curfile+1
                if (curfile>nfiles) then
                    stop "Ran out of files to process!"
                endif
                curstep=curstep-steps_in_file !instead of setting=1, this way we can set an arbitrary starting point multiple files in
                call io_getdims(file_list(curfile),options%pvar, dims)
                if (dims(1)==3) then
                    steps_in_file=1
                else
                    steps_in_file=dims(dims(1)+1) !dims(1) = ndims; dims(ndims+1)=ntimesteps
                endif
            enddo
            
            smoothing_window = min(max(int(options%smooth_wind_distance/domain%dx),1),size(domain%lat,1)/5)
            if (options%debug) write(*,*) "  Smoothing winds over ",trim(str(smoothing_window))," grid cells"
        endif
!       load the restart file
        if (options%restart) then
            call load_restart_file(domain,options%restart_file,options%restart_step_in_file)
            write(*,*) "Reading winds"
            write(*,*) "  timestep:",trim(str(curstep)),"  from file:",trim(file_list(curfile))
            if (options%external_winds) then
                call ext_winds_init(domain,bc,options)
            elseif (options%lt_options%remove_lowres_linear) then
                ! remove the low-res linear wind perturbation field 
                call remove_linear_winds(domain,bc,options,file_list(curfile),curstep)
                call smooth_wind(domain%u,smoothing_window,3)
                call smooth_wind(domain%v,smoothing_window,3)
            elseif (options%mean_winds) then
                call mean_winds(domain,file_list(curfile),curstep,options)
            else
                call read_var(domain%u, file_list(curfile),options%uvar,  &
                                bc%u_geo%geolut,bc%u_geo%vert_lut,curstep,boundary_value, &
                                options)
                call read_var(domain%v, file_list(curfile),options%vvar,  &
                                bc%v_geo%geolut,bc%v_geo%vert_lut,curstep,boundary_value, &
                                options)
                call smooth_wind(domain%u,smoothing_window,3)
                call smooth_wind(domain%v,smoothing_window,3)
            endif
            domain%pii=(domain%p/100000.0)**(Rd/cp)
            domain%rho=domain%p/(Rd*domain%th*domain%pii) ! kg/m^3
            ! note, ideally update_winds needs an initializer instead of being called when we are doing a restart
            ! alternatively we could pass an init_only=True into update_winds
            ! so that the linear wind solution isn't calculated for a restart step
            ! however, we ned to call update_winds to allow all variables to be initialized properly
            call update_winds(domain,options)
            ! we can't just call balance_uvw because the linear winds may need to be initialized
            !call balance_uvw(domain,options)
            
            call output_init(domain,options)
            call write_domain(domain,options,-1)
        else
!           else load data from the first Boundary conditions file
            boundary_value=.False.
            nx=size(domain%p,1)
            ny=size(domain%p,3)
            if (options%external_winds) then
                call ext_winds_init(domain,bc,options)
                ! call smooth_wind(domain%u,1,3)
                ! call smooth_wind(domain%v,1,3)
            elseif (options%lt_options%remove_lowres_linear) then
                ! remove the low-res linear wind perturbation field 
                call remove_linear_winds(domain,bc,options,file_list(curfile),curstep)
                call smooth_wind(domain%u,smoothing_window,3)
                call smooth_wind(domain%v,smoothing_window,3)
            elseif (options%mean_winds) then
                call mean_winds(domain,file_list(curfile),curstep,options)
            else
                call read_var(domain%u, file_list(curfile),options%uvar,  &
                                bc%u_geo%geolut,bc%u_geo%vert_lut,curstep,boundary_value, &
                                options)
                call read_var(domain%v, file_list(curfile),options%vvar,  &
                                bc%v_geo%geolut,bc%v_geo%vert_lut,curstep,boundary_value, &
                                options)
                call smooth_wind(domain%u,smoothing_window,3)
                call smooth_wind(domain%v,smoothing_window,3)
            endif
            call read_var(domain%p,    file_list(curfile),   options%pvar,   &
                            bc%geolut, bc%vert_lut, curstep, boundary_value, &
                            options,   bc%lowres_z,domain%z, interp_vertical=.True.)
            call read_var(domain%th,   file_list(curfile),   options%tvar,   &
                            bc%geolut, bc%vert_lut, curstep, boundary_value, &
                            options)
            call read_var(domain%qv,   file_list(curfile),   options%qvvar,  &
                            bc%geolut, bc%vert_lut, curstep, boundary_value, &
                            options)
            if (trim(options%qcvar)/="") then
                call read_var(domain%cloud,file_list(curfile),   options%qcvar,  &
                                bc%geolut, bc%vert_lut, curstep, boundary_value, &
                                options)
            else
                if (options%debug) print*, "No Cloud water content specified"
                domain%cloud=0
            endif
            if (trim(options%qivar)/="") then
                call read_var(domain%ice,  file_list(curfile),   options%qivar,  &
                                bc%geolut, bc%vert_lut, curstep, boundary_value, &
                                options)
            else
                if (options%debug) print*, "No Cloud ice specified"
                domain%ice=0
            endif
            
            if (options%physics%landsurface==kLSM_BASIC) then
                call read_2dvar(domain%sensible_heat,file_list(curfile),options%shvar,  bc%geolut,curstep)
                call read_2dvar(domain%latent_heat,  file_list(curfile),options%lhvar,  bc%geolut,curstep)
                
                if (options%physics%boundarylayer==kPBL_BASIC) then
                    if (trim(options%pblhvar)/="") then
                        call read_2dvar(domain%pbl_height,   file_list(curfile),options%pblhvar,bc%geolut,curstep)
                    endif
                endif
                ! NOTE, this is a kludge to prevent the model from sucking more moisture out of the lower model layer than exists
                where(domain%latent_heat<0) domain%latent_heat=0
            endif
            
            if (options%physics%radiation==kRA_BASIC) then
                if (trim(options%swdown_var)/="") then
                    call read_2dvar(domain%swdown,  file_list(curfile),options%swdown_var,  bc%geolut,curstep)
                endif
                if (trim(options%lwdown_var)/="") then
                    call read_2dvar(domain%lwdown,  file_list(curfile),options%lwdown_var,  bc%geolut,curstep)
                endif
            endif
            if (trim(options%sst_var)/="") then
                call read_2dvar(domain%sst,  file_list(curfile),options%sst_var,  bc%geolut,curstep)
            endif
            
            call update_pressure(domain%p,bc%lowres_z,domain%z)
            
            nz=size(domain%th,2)
            domainsize=size(domain%th,1)*size(domain%th,3)
            if (options%mean_fields) then
                do i=1,nz
                    domain%th(:,i,:)=sum(domain%th(:,i,:))/domainsize
                    domain%qv(:,i,:)=sum(domain%qv(:,i,:))/domainsize
                    domain%cloud(:,i,:)=sum(domain%cloud(:,i,:))/domainsize
                    domain%ice(:,i,:)=sum(domain%ice(:,i,:))/domainsize
                enddo
            endif

            domain%pii=(domain%p/100000.0)**(Rd/cp)
            domain%rho=domain%p/(Rd*domain%th*domain%pii) ! kg/m^3
            call update_winds(domain,options)
        endif
        
        ! calculate znu and znw from domain pressure variable now that we have it
        call init_znu(domain)
        
    end subroutine bc_init


    !>------------------------------------------------------------
    !!  Same as update_dxdt but only for the edges of the domains
    !!
    !!  This is used for fields that are calculated/updated internally
    !!  by the model physics (e.g. temperature and moisture)
    !!  In the output dxdt variable, the first dimension is the z axis, 
    !!  The second dimension is either X or Y (which ever is specified)
    !!  And the third dimension specifies the boundary it applies to
    !!  1=left, 2=right, 3=bottom, 4=top
    !!
    !! @param dx_dt Change in variable X between time periods d1 and d2
    !! @param d1    Value of field X at time period 1
    !! @param d2    Value of field X at time period 2
    !! @retval dx_dt This field is updated along the boundaries to be (d1-d2)
    !!
    !!------------------------------------------------------------
    subroutine update_edges(dx_dt,d1,d2)
        implicit none
        real,dimension(:,:,:), intent(inout) :: dx_dt
        real,dimension(:,:,:), intent(in) :: d1,d2
        integer :: nx,nz,ny,i

        nx=size(d1,1)
        nz=size(d1,2)
        ny=size(d1,3)
        do i=1,nz
            dx_dt(i,:ny,1)=d1(1,i,:) -d2(1,i,:)
            dx_dt(i,:ny,2)=d1(nx,i,:)-d2(nx,i,:)
            dx_dt(i,:nx,3)=d1(:,i,1) -d2(:,i,1)
            dx_dt(i,:nx,4)=d1(:,i,ny)-d2(:,i,ny)
        enddo
    end subroutine update_edges
    
    
    !>------------------------------------------------------------
    !!  Calculate changes between the current boundary conditions and the time step boundary conditions
    !!
    !!  These changes are used to linearly shift all fields between the two times. 
    !!
    !! @param bc        Model boundary forcing conditions data structure (future conditions and changes)
    !! @param domain    Model domain data structure (current conditions) 
    !!
    !!------------------------------------------------------------
    subroutine update_dxdt(bc,domain)
        implicit none
        type(bc_type), intent(inout) :: bc
        type(domain_type), intent(in) :: domain
        
        ! the only "standard" variable that is always computed is p
        bc%dp_dt=bc%next_domain%p-domain%p
        
        ! NOTE these are only used if lsm option = 1, a bunch of wasted zeros otherwise
        bc%dsh_dt  =bc%next_domain%sensible_heat-domain%sensible_heat
        bc%dlh_dt  =bc%next_domain%latent_heat-domain%latent_heat
        ! only if lsm=1 and PBL option = 1
        bc%dpblh_dt=bc%next_domain%pbl_height-domain%pbl_height
        ! only if rad=1 and lsm\=1
        bc%dsw_dt  =bc%next_domain%swdown-domain%swdown
        bc%dlw_dt  =bc%next_domain%lwdown-domain%lwdown
        ! only if water=1
        bc%dsst_dt  =bc%next_domain%sst-domain%sst
        
        ! these fields are only updated along the edges of the domain
        call update_edges(bc%dth_dt,bc%next_domain%th,domain%th)
        call update_edges(bc%dqv_dt,bc%next_domain%qv,domain%qv)
        call update_edges(bc%dqc_dt,bc%next_domain%cloud,domain%cloud)
    end subroutine update_dxdt

    !>------------------------------------------------------------
    !!  Calculate changes in winds between the current and future conditions
    !!
    !!  These changes are used to linearly shift all fields between the two times. 
    !!
    !! @param bc        Model boundary forcing conditions data structure (future conditions and changes)
    !! @param domain    Model domain data structure (current conditions) 
    !!
    !!------------------------------------------------------------
    subroutine update_dwinddt(bc,domain)
        implicit none
        type(bc_type), intent(inout) :: bc
        type(domain_type), intent(in) :: domain
        
        bc%du_dt=bc%next_domain%u-domain%u
        bc%dv_dt=bc%next_domain%v-domain%v
    end subroutine update_dwinddt
    
    !>------------------------------------------------------------
    !!  Adjust the pressure field for the vertical shift between the low and high-res domains
    !!
    !!  Ideally this should include temperature... but it isn't entirely clear
    !!  what it would mean to do that, what temperature do you use? Current even though you are adjusting future P?
    !!  Alternatively, could adjust input pressure to SLP with future T then adjust back to elevation with current T?
    !!  Currently if T is supplied, it uses the mean of the high and low-res T to split the difference.
    !!  Equations from : http://www.wmo.int/pages/prog/www/IMOP/meetings/SI/ET-Stand-1/Doc-10_Pressure-red.pdf
    !!  excerpt from CIMO Guide, Part I, Chapter 3 (Edition 2008, Updated in 2010) equation 3.2
    !!  http://www.meteormetrics.com/correctiontosealevel.htm
    !!
    !! @param pressure  The pressure field to be adjusted
    !! @param z_lo      The 3D vertical coordinate of the input pressures
    !! @param z_hi      The 3D vertical coordinate of the computed/adjusted pressures
    !! @param lowresT   OPTIONAL 3D temperature field of the input pressures
    !! @param lowresT   OPTIONAL 3D temperature field of the computed/adjusted pressures
    !! @retval pressure The pressure field after adjustment
    !!
    !!------------------------------------------------------------
    subroutine update_pressure(pressure,z_lo,z_hi, lowresT, hiresT)
        implicit none
        real,dimension(:,:,:), intent(inout) :: pressure
        real,dimension(:,:,:), intent(in) :: z_lo,z_hi
        real,dimension(:,:,:), intent(in), optional :: lowresT, hiresT
        ! local variables
        real,dimension(:),allocatable::slp !sea level pressure [Pa]
        ! vapor pressure, change in height, change in temperature with height and mean temperature
        real,dimension(:),allocatable:: dz, tmean !, e, dTdz
        integer :: nx,ny,nz,i,j
        nx=size(pressure,1)
        nz=size(pressure,2)
        ny=size(pressure,3)
        
        if (present(lowresT)) then
            !$omp parallel shared(pressure, z_lo,z_hi, lowresT, hiresT) &
            !$omp private(i,j, dz, tmean) firstprivate(nx,ny,nz)  !! private(e, dTdz)
            allocate(dz(nx))
            allocate(tmean(nx))
            ! allocate(e(nx))
            ! allocate(dTdz(nx))
            !$omp do 
            do j=1,ny
                ! is an additional loop over z more cache friendly? 
                do i=1,nz
                    ! vapor pressure
!                     e = qv(:,:,j) * pressure(:,:,j) / (0.62197+qv(:,:,j))
                    ! change in elevation (note reverse direction from "expected" because the formula is an SLP reduction)
                    dz   = (z_lo(:,i,j) - z_hi(:,i,j))
                    ! lapse rate (not sure if this should be positive or negative)
                    ! dTdz = (loresT(:,:,j) - hiresT(:,:,j)) / dz
                    ! mean temperature between levels
                    if (present(hiresT)) then
                        tmean= (hiresT(:,i,j) + lowresT(:,i,j)) / 2
                    else
                        tmean= lowresT(:,i,j)
                    endif                
                    ! slp= ps*np.exp(((g/R)*Hp) / (ts - a*Hp/2.0 + e*Ch))
                    pressure(:,i,j) = pressure(:,i,j) * exp( ((gravity/Rd) * dz) / tmean )   !&
!                                         (tmean + (e * 0.12) ) )
                    ! alternative formulation M=0.029, R=8.314?
                    ! p= p0*(t0/(t0+dtdz*z))**((g*M)/(R*dtdz))
                    ! do i=1,nz
                    !     pressure(:,i,j) = pressure(:,i,j)*(t0/(tmean(:,i)+dTdz(:,i)*z))**((g*M)/(R*dtdz))
                    ! enddo
                enddo
            enddo
            !$omp end do
            deallocate(dz, tmean)
            ! deallocate(e, dTdz)
            !$omp end parallel
        else
            ! this is pretty foolish to convert to sea level pressure and back... should be done in one step
            ! just need to test that the relationship can work that way h=(zhi-zlo)
            ! this doesn't get used much though (only from bc_init) so it doesnt seem worth the time...
            !$omp parallel shared(pressure, z_lo,z_hi) &
            !$omp private(slp,i,j) firstprivate(nx,ny,nz)
            allocate(slp(nx))
            !$omp do 
            do j=1,ny
                do i=1,nz
                    ! slp = pressure(:,i,j) / (1 - 2.25577E-5 * z_lo(:,i,j))**5.25588
                    pressure(:,i,j) = pressure(:,i,j) * (1 - 2.25577e-5 * (z_hi(:,i,j)-z_lo(:,i,j)))**5.25588
                enddo
            enddo
            !$omp end do
            deallocate(slp)
            !$omp end parallel
        endif
    end subroutine update_pressure
    
!>------------------------------------------------------------
!!  Update the external wind field
!!
!!  Read U and V and rotate into the domain 3D grid, store in bc%next_domain
!!
!! @param bc        Model forcing boundary conditions data structure
!! @param options   Model Options structure 
!!
!!------------------------------------------------------------
    subroutine update_ext_winds(bc,options)
        implicit none
        type(bc_type),intent(inout)::bc
        type(options_type),intent(in)::options
        integer,dimension(io_maxDims)::dims !note, io_maxDims is included from io_routines.
        logical :: use_boundary,use_interior
        ! MODULE variables : ext_winds_ curstep, curfile, nfiles, steps_in_file, file_list
        
        ext_winds_curstep=ext_winds_curstep+1
        if (ext_winds_curstep>ext_winds_steps_in_file) then
            ext_winds_curfile=ext_winds_curfile+1
            ext_winds_curstep=1
            call io_getdims(ext_winds_file_list(ext_winds_curfile),options%uvar, dims)
            if (dims(1)==3) then
                ext_winds_steps_in_file=1
            else
                ext_winds_steps_in_file=dims(dims(1)+1) !dims(1) = ndims
            endif
        endif
        if (ext_winds_curfile>ext_winds_nfiles) then
            stop "Ran out of files to process!"
        endif
        
        use_interior=.False.
        use_boundary=.True.
        call read_var(bc%next_domain%u,    ext_winds_file_list(ext_winds_curfile),options%uvar, &
                      bc%ext_winds%u_geo%geolut,bc%u_geo%vert_lut,ext_winds_curstep,use_interior,options)
        call read_var(bc%next_domain%v,    ext_winds_file_list(ext_winds_curfile),options%vvar, &
                      bc%ext_winds%v_geo%geolut,bc%v_geo%vert_lut,ext_winds_curstep,use_interior,options)
        call rotate_ext_wind_field(bc%next_domain,bc%ext_winds)
    
    end subroutine update_ext_winds
    
    !>------------------------------------------------------------
    !!  Read in the next timestep of input data and apply to the dXdt grids as appropriate. 
    !!  
    !!  Nearly identical to bc_init. Primary differences are that: 
    !!  Restart fields are not checked for and read, and data are stored in bc%next_domain instead of domain
    !!  In addition, the time_varying_z condition is handled here with an additional interpolation. 
    !!  Finally, current domain conditions are copied into next_domain (after computing dxdt values)
    !!  for use in the linear_wind calculations (i.e. linear wind field for the next time step is calculated
    !!  based on atmospheric stability conditions for the current time step). 
    !!
    !! @param domain    Model Domain
    !! @param bc        Model forcing boundary conditions data structure
    !! @param options   Model Options structure 
    !!
    !!------------------------------------------------------------
    subroutine bc_update(domain,bc,options)
        implicit none
        type(domain_type),intent(inout)::domain
        type(bc_type),intent(inout)::bc
        type(options_type),intent(in)::options
        integer,dimension(io_maxDims)::dims !note, io_maxDims is included from io_routines.
        type(bc_type) :: newbc ! just used for updating z coordinate
        real, allocatable, dimension(:,:,:) :: zbase ! may be needed to temporarily store PHB data
        logical :: use_boundary,use_interior
        integer::i,nz,nx,ny
        ! MODULE variables : curstep, curfile, nfiles, steps_in_file, file_list
        
        if (.not.options%ideal) then
            curstep=curstep+1
            do while (curstep>steps_in_file)
                curfile=curfile+1
                if (curfile>nfiles) then
                    stop "Ran out of files to process!"
                endif
                curstep=curstep-steps_in_file !instead of setting=1, this way we can set an arbitrary starting point multiple files in
                call io_getdims(file_list(curfile),options%pvar, dims)
                if (dims(1)==3) then
                    steps_in_file=1
                else
                    steps_in_file=dims(dims(1)+1) !dims(1) = ndims; dims(ndims+1)=ntimesteps
                endif
            enddo
        endif
        use_interior=.False. ! this is passed to the use_boundary flag in geo_interp
        use_boundary=.True.
        
        if (options%time_varying_z) then
            ! read in the updated vertical coordinate
            if (allocated(newbc%z)) deallocate(newbc%z)
            call io_read3d(file_list(curfile), options%zvar, newbc%z, curstep)
            nx=size(newbc%z,1)
            ny=size(newbc%z,2)
            nz=size(newbc%z,3)
            if (trim(options%zbvar)/="") then
                call io_read3d(file_list(curfile),options%zbvar, zbase, curstep)
                newbc%z=newbc%z+zbase
                deallocate(zbase)
            endif
            if (options%z_is_geopotential) then
                newbc%z=newbc%z / gravity
            endif
            if (options%z_is_on_interface) then
                newbc%z(:,:,1:nz-1)=(newbc%z(:,:,1:nz-1) + newbc%z(:,:,2:nz))/2
            endif
            
            ! now simply generate a look up table to convert the current z coordinate to the original z coordinate
            call vLUT_forcing(bc,newbc)
            ! set a maximum on z so we don't try to interpolate data above mass grid
            where(newbc%vert_lut%z==nz) newbc%vert_lut%z=nz-1
            nz=size(newbc%z,3)
            ! generate a new high-res z dataset as well (for pressure interpolations)
            ! call geo_interp(bc%lowres_z, reshape(newbc%z,[nx,nz,ny],order=[1,3,2]), bc%geolut,use_interior)
            
        endif
        
        ! first read in and handle winds
        if (options%external_winds) then
            call update_ext_winds(bc,options)
!           call smooth_wind(bc%next_domain%u,1,3)
!           call smooth_wind(bc%next_domain%v,1,3)
        elseif (options%lt_options%remove_lowres_linear) then
            call remove_linear_winds(bc%next_domain,bc,options,file_list(curfile),curstep)
            call smooth_wind(bc%next_domain%u,smoothing_window,3)
            call smooth_wind(bc%next_domain%v,smoothing_window,3)
        elseif (options%mean_winds) then
            call mean_winds(bc%next_domain,file_list(curfile),curstep,options)
        else
            ! general case, just read in u and v data
            call read_var(bc%next_domain%u,  file_list(curfile), options%uvar, &
                          bc%u_geo%geolut,   bc%u_geo%vert_lut,  curstep,      &
                          use_interior, options, time_varying_zlut=newbc%vert_lut)
            call read_var(bc%next_domain%v,  file_list(curfile), options%vvar, &
                          bc%v_geo%geolut,   bc%v_geo%vert_lut,  curstep,      &
                          use_interior, options, time_varying_zlut=newbc%vert_lut)
            call smooth_wind(bc%next_domain%u,smoothing_window,3)
            call smooth_wind(bc%next_domain%v,smoothing_window,3)
        endif
        
        ! now read in remaining variables
        ! for pressure do not apply vertical interpolation on IO, we will adjust it more accurately
        call read_var(bc%next_domain%p,       file_list(curfile), options%pvar,   &
                      bc%geolut, bc%vert_lut, curstep, use_interior,              &
                      options, time_varying_zlut=newbc%vert_lut, interp_vertical=.True.)
        ! for pressure adjustment, we need temperature on the original model grid, 
        ! so read it without vertical interpolation (and for interior points too)
        call read_var(bc%next_domain%th,      file_list(curfile), options%tvar,   &
                      bc%geolut, bc%vert_lut, curstep, use_interior,              &
                      options, time_varying_zlut=newbc%vert_lut, interp_vertical=.True. )
        ! for pressure update we need real temperature, not potential t to compute an exner function
        ! bc%next_domain%pii=(bc%next_domain%p/100000.0)**(Rd/cp)
        ! now update pressure using the high res T field
        call update_pressure(bc%next_domain%p,bc%lowres_z,domain%z)!,  & 
                            !  lowresT = bc%next_domain%th * bc%next_domain%pii, &
                            !  hiresT  = domain%th * domain%pii)
        
        ! not necessary as long as we are interpolating above
        call read_var(bc%next_domain%th,      file_list(curfile), options%tvar,   &
                      bc%geolut, bc%vert_lut, curstep, use_boundary,              &
                      options, time_varying_zlut=newbc%vert_lut)
                      
        call read_var(bc%next_domain%qv,      file_list(curfile), options%qvvar,  &
                      bc%geolut, bc%vert_lut, curstep, use_boundary,              &
                      options, time_varying_zlut=newbc%vert_lut)
                      
        if (trim(options%qcvar)/="") then
            call read_var(bc%next_domain%cloud,   file_list(curfile), options%qcvar,  &
                          bc%geolut, bc%vert_lut, curstep, use_boundary,              &
                          options, time_varying_zlut=newbc%vert_lut)
        endif
        if (trim(options%qivar)/="") then
            call read_var(bc%next_domain%ice,     file_list(curfile), options%qivar,  &
                          bc%geolut, bc%vert_lut, curstep, use_boundary,              &
                          options, time_varying_zlut=newbc%vert_lut)
        endif

        ! finally, if we need to read in land surface forcing read in those 2d variables as well. 
        if (options%physics%landsurface==kLSM_BASIC) then
            call read_2dvar(bc%next_domain%sensible_heat,file_list(curfile),options%shvar,  bc%geolut,curstep)
            call read_2dvar(bc%next_domain%latent_heat,  file_list(curfile),options%lhvar,  bc%geolut,curstep)
            ! note this is nested in the landsurface=LSM_BASIC condition, because that is the only time it makes sense. 
            if (options%physics%boundarylayer==kPBL_BASIC) then
                if (trim(options%pblhvar)/="") then
                    call read_2dvar(bc%next_domain%pbl_height,   file_list(curfile),options%pblhvar,bc%geolut,curstep)
                endif
            endif
            ! NOTE, this is a kludge to prevent the model from sucking more moisture out of the lower model layer than exists
            where(domain%latent_heat<0) domain%latent_heat=0
        endif
        
        if (options%physics%radiation==kRA_BASIC) then
            if (trim(options%swdown_var)/="") then
                call read_2dvar(bc%next_domain%swdown,  file_list(curfile),options%swdown_var,  bc%geolut,curstep)
            endif
            if (trim(options%lwdown_var)/="") then
                call read_2dvar(bc%next_domain%lwdown,  file_list(curfile),options%lwdown_var,  bc%geolut,curstep)
            endif
        endif
        
        if (trim(options%sst_var)/="") then
            call read_2dvar(bc%next_domain%sst,  file_list(curfile),options%sst_var,  bc%geolut,curstep)
        endif
        
        ! if we want to supply mean forcing fields on the boundaries, compute those here. 
        if (options%mean_fields) then
            call mean_boundaries(bc%next_domain%th)
            call mean_boundaries(bc%next_domain%qv)
            call mean_boundaries(bc%next_domain%cloud)
            call mean_boundaries(bc%next_domain%ice)
        endif
        
        
        ! update scalar dXdt tendency fields first so we can then overwrite them with 
        ! the current model state
        call update_dxdt(bc,domain)
        
        ! we need the internal values of these fields to be in sync with the high res model
        ! for the linear wind calculations... albeit these are for time t, and winds are for time t+1
        bc%next_domain%qv=domain%qv
        bc%next_domain%th=domain%th
        bc%next_domain%cloud=domain%cloud + domain%ice + domain%qrain + domain%qsnow
        
        ! these are required by update_winds for most options
        bc%next_domain%pii=(bc%next_domain%p/100000.0)**(Rd/cp)
        bc%next_domain%rho=bc%next_domain%p/(Rd*domain%th*bc%next_domain%pii) ! kg/m^3
        
        
        call update_winds(bc%next_domain,options)
        ! copy it to the primary domain for output purposes (could also be used for convection or blocking parameterizations?)
        if (options%physics%windtype==kWIND_LINEAR) then
            domain%nsquared=bc%next_domain%nsquared
        endif
        
        ! then updated with wind dXdt fields after updating them
        call update_dwinddt(bc,domain)
    end subroutine bc_update
end module boundary_conditions
