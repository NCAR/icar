!>----------------------------------------------------------
!!  Handle reading and writing Linear wind Look up tables to/from disk
!! 
!!  Look up tables are written from disk in a documented netCDF file (e.g. with all LUT options)
!!  When a LUT is read from disk, the options stored in that file are checked against the options
!!  specified for the current run. If they are compatible (x, y, and z shape) the LUT will be read
!!  If the options for the current run are not the same as those in the file, a warning is printed
!!  but and an error will be thrown, the model will continue, it will just use the options
!!  stored in the file so the LUT does not need to be regenerated. 
!! 
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module linear_theory_lut_disk_io
    
    use data_structures
    use netcdf
    use io_routines,        only: file_exists, io_read, check
    implicit none
    
    private
    public :: read_lut
    public :: write_lut

contains
!>----------------------------------------------------------
!!  Write linear wind look up tables to disk along with the parameters for that table
!! 
!!  @detail
!!  Takes a filename, the LUTs, and the options structure that defines the LUT parameters
!!  A file is created with appropriate dimension names and attributes, along with variables
!!  for each of the relevant. 
!!
!!  @sideeffect
!!  If the parameters do not match the parameters in the LUT file, the parameters structure is updated to match
!! 
!!  @param filename Name of LUT file to write
!!  @param uLUT     Array to hold u wind look up table
!!  @param vLUT     Array to hold v wind look up table
!!  @param options  Linear wind options to compare with the LUT file and update if necessary
!!  @retval error   0 on success other values indicate and error
!!
!!----------------------------------------------------------
function write_lut(filename, uLUT, vLUT, options) result(error)
    implicit none
    character(len=*), intent(in) :: filename
    real, dimension(:,:,:,:,:,:), intent(in) :: uLUT, vLUT
    type(lt_options_type), intent(in) :: options
    integer :: error
    error=0
    
    if (file_exists(filename)) then
        write(*,*) "WARNING: Linear Theory look-up-table file already exists."
        error = 1
        return
    endif
    ! allocate(rev_u_LUT(n_spd_values,n_dir_values,n_nsq_values,nxu,nz,ny))
    call write_var(filename,"uLUT",uLUT, dims=["nspd","ndir","nnsq","nxu","nz","ny"]))
    call write_var(filename,"vLUT",vLUT, dims=["nspd","ndir","nnsq","nx","nz","nyv"])
    
    call add_attribute(filename, "dirmax", options%dirmax)
    call add_attribute(filename, "dirmin", options%dirmin)
    call add_attribute(filename, "spdmax", options%spdmax)
    call add_attribute(filename, "spdmin", options%spdmin)
    call add_attribute(filename, "nsqmax", options%nsqmax)
    call add_attribute(filename, "nsqmin", options%nsqmin)
    call add_attribute(filename, "n_dir_values", options%n_dir_values)
    call add_attribute(filename, "n_nsq_values", options%n_nsq_values)
    call add_attribute(filename, "n_spd_values", options%n_spd_values)
    
end function write_lut

!>----------------------------------------------------------
!!  Write linear wind look up tables to disk along with the parameters for that table
!! 
!!  @detail
!!  Takes a filename, the allocatable LUT arrays, and the options structure that defines the LUT parameters
!!  A file is read and checked for appropriate dimensions, along with parameters. 
!!
!!  @sideeffect
!!  If the parameters do not match the parameters in the LUT file, the parameters structure is updated to match
!! 
!!  @param filename Name of LUT file to write
!!  @param uLUT     Array to hold u wind look up table
!!  @param vLUT     Array to hold v wind look up table
!!  @param dims     1D array containtin dimension sizes of the model domain
!!  @param options  Linear wind options to compare with the LUT file and update if necessary
!!  @retval error   0 if there were no errors, 
!!                  1 if parameters did not match, 
!!                  2 if dimensions don't match,
!!                  3 if file doesn't exist
!!
!!----------------------------------------------------------
function read_LUT(filename, uLUT, vLUT, dims, options) result(error)
    implicit none
    character(len=*), intent(in) :: filename
    real, allocatable, dimension(:,:,:,:,:,:), intent(inout) :: uLUT, vLUT
    integer, dimension(3,2), intent(in) :: dims
    type(lt_options_type), intent(inout) :: options
    integer :: error
    error=0
    
    if (.not.file_exists(filename)) then
        error = 3
        return
    endif
    if (.not.dims_match(filename, dims)) then
        error = 2
        return
    endif
    
    ! this is the slow part where we actually read in all the data
    call io_read(filename,"uLUT",uLUT)
    call io_read(filename,"vLUT",vLUT)

    ! next read in the parameter values and test them against the input parameters
    ! if the file and input parameters do not match, adjust the input parameter to match
    ! and set the error flag to indicate parameters did not match and print a warning
    call io_read_attribute(filename, "dirmax", dirmax)
    options%dirmax = check_attribute(dirmax, options%dirmax, "dirmax", error)
    
    call io_read_attribute(filename, "dirmin", dirmin)
    options%dirmin = check_attribute(dirmin, options%dirmin, "dirmin", error)
    
    call io_read_attribute(filename, "spdmax", spdmax)
    options%spdmax = check_attribute(spdmax, options%spdmax, "spdmax", error)
    
    call io_read_attribute(filename, "spdmin", spdmin)
    options%spdmin = check_attribute(spdmin, options%spdmin, "spdmin", error)
    
    call io_read_attribute(filename, "nsqmax", nsqmax)
    options%nsqmax = check_attribute(nsqmax, options%nsqmax, "nsqmax", error)
    
    call io_read_attribute(filename, "nsqmin", nsqmin)
    options%nsqmin = check_attribute(nsqmin, options%nsqmin, "nsqmin", error)
    
    call io_read_attribute(filename, "n_dir_values", n_dir_values)
    options%n_dir_values = check_attribute(n_dir_values, options%n_dir_values, "n_dir_values", error)
    
    call io_read_attribute(filename, "n_nsq_values", n_nsq_values)
    options%n_nsq_values = check_attribute(n_nsq_values, options%n_nsq_values, "n_nsq_values", error)
    
    call io_read_attribute(filename, "n_spd_values", n_spd_values)
    options%n_spd_values = check_attribute(n_spd_values, options%n_spd_values, "n_spd_values", error)
    
    ! return error
end function read_LUT


function write_var(filename,varname,data_out, dimnames, open_new_file) result(error)
        implicit none
        ! We are writing 6D data, a nx, nz, ny, na, nb, nc grid. 
        integer, parameter :: ndims = 6
        
        ! This is the name of the file and variable we will write. 
        character(len=*), intent(in) :: filename, varname
        real,intent(in) :: data_out(:,:,:,:,:,:)
        character(len=*), optional, dimension(ndims), intent(in) :: dimnames
        logical, optional, intent(in) :: open_new_file
        integer :: error
        
        ! This will be the netCDF ID for the file and data variable.
        integer :: ncid, varid,dimids(ndims)
        character(len=*), dimension(ndims) :: dims
        integer :: dim_length
        
        if (present(dimnames)) then
            dims = dimnames
        else
            dims = ["x","y","z","a","b","c"]
        endif
        
        ncid=-1
        ! Open the file. NF90_NOCLOBBER tells netCDF we want append to existing files
        if (present(open_new_file)) then
            if (open_new_file) then
                call check( nf90_create(filename, NF90_CLOBBER, ncid), filename)
            endif
        endif
        if (ncid/=-1) then
            call check( nf90_create(filename, NF90_NOCLOBBER, ncid), "write_var:creating file:"//filename )
        endif
        ! define the dimensions or test that the existing dim matches the size it should
        do i=1,ndims
            dim_length = size(data_out, i)
            error = error + (2**i) * setup_dim(ncid, dims(i), dimids(i), dim_length)
        end do
        
        varid = -1
        error = error + 2**(ndims+1) * setup_var(ncid, varname, varid, dimids, dims)
        
        ! End define mode. This tells netCDF we are done defining metadata.
        call check( nf90_enddef(ncid) )
        
        !write the actual data to the file
        if (varid /= -1) then
            call check( nf90_put_var(ncid, varid, data_out), trim(filename)//":"//trim(varname))
        endif
        
        ! Close the file, freeing all resources.
        call check( nf90_close(ncid), filename)
        
        ! return error
end function write_var

subroutine add_attribute(filename, att_name, att_value, varname)
    implicit none
    character(len=*), intent(in)           :: filename
    character(len=*), intent(in)           :: att_name
    real,             intent(in)           :: att_value
    character(len=*), intent(in), optional :: varname
    
    integer :: ncid
    integer :: varid
    
    
    
end subroutine add_attribute

function setup_var(variables) result(error)
    implicit none
    
    call check(nf90_inq_varid(ncid, "time", varid),                 trim(filename)//" : time")
    call check(nf90_inquire_variable(ncid, varid, dimids = dims),   trim(filename)//" : time dims")
    call check(nf90_inquire_dimension(ncid, dims(1), len = ntimes), trim(filename)//" : inq time dim")
end function setup_var

function setup_dim(ncid, dimname, dimid, n) result(error)
    implicit none
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: dimname
    integer, intent(inout) :: dimid
    integer, intent(in), optional :: n
    integer :: error
    integer :: ncerror, dim_length
    
    error = 0
    
    ncerror = NF90_INQ_DIMID(ncid, dimname, dimid)
    if (ncerror/=NF90_NOERR) then
        call check( nf90_def_dim(ncid, dimname, n, dimid), "setup_dim:Creating dim:"//dimname )
    endif
    
    if (present(n)) then
        call check( nf90_inquire_dimension(ncid, dimid, len = dim_length), "setup_dim:Reading dim:"//dimname)
        
        if (dim_length/=n) then
            error = 1
        endif
    endif
    
    ! return error
end function setup_dim

function check_attribute(test_value, default_value, value_name, error) result(output_value)
    implicit none
    real, intent(in) :: test_value, default_value
    character(len=*), intent(in) :: value_name
    integer, intent(inout) :: error
    
    real :: output_value
    
    if (default_value /= test_value) then
        error = 1
        write(*,*) "WARNING parameter option: "//trim(value_name),options%dirmax, " did not match file value: ",dirmax
        options%dirmax=dirmax
    endif
end function check_attribute
