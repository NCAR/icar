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
    use string,         only: str
    use netcdf
    use io_routines,    only: file_exists, io_read, check, io_read_attribute, io_add_attribute, io_maxDims, io_getdims
    implicit none
    
    private
    public :: read_lut
    public :: write_lut
    
    interface check_attribute
        module procedure check_real_attribute, check_integer_attribute
    end interface check_attribute

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
            if (.not.options%overwrite_lt_lut) then
                write(*,*) "WARNING: Existing file will not be overwritten."
                return
            endif
            write(*,*) "WARNING: Existing file *will* be overwritten."
        endif
        ! allocate(rev_u_LUT(n_spd_values,n_dir_values,n_nsq_values,nxu,nz,ny))
        error = write_var(filename,"uLUT",uLUT, dimnames=[character(len=4) :: "nspd","ndir","nnsq","nxu","nz","ny"], open_new_file=.True.)
        if (error/=0) then
            write(*,*) "Error writing LUT file:"//trim(filename)
            return
        endif
        error = write_var(filename,"vLUT",vLUT, dimnames=[character(len=4) :: "nspd","ndir","nnsq","nx","nz","nyv"])
        if (error/=0) then
            write(*,*) "Error writing LUT file:"//trim(filename)
            return
        endif
        
        call io_add_attribute(filename, "dirmax", options%dirmax)
        call io_add_attribute(filename, "spdmax", options%spdmax)
        call io_add_attribute(filename, "spdmin", options%spdmin)
        call io_add_attribute(filename, "dirmin", options%dirmin)
        call io_add_attribute(filename, "nsqmax", options%nsqmax)
        call io_add_attribute(filename, "nsqmin", options%nsqmin)
        call io_add_attribute(filename, "n_dir_values", options%n_dir_values)
        call io_add_attribute(filename, "n_nsq_values", options%n_nsq_values)
        call io_add_attribute(filename, "n_spd_values", options%n_spd_values)
        
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
        type(lt_options_type), intent(in) :: options
        ! type(lt_options_type), intent(inout) :: options
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
        error = check_attribute(filename, options%dirmax, "dirmax")
        error = check_attribute(filename, options%dirmin, "dirmin")
        error = check_attribute(filename, options%spdmax, "spdmax")
        error = check_attribute(filename, options%spdmin, "spdmin")
        error = check_attribute(filename, options%nsqmax, "nsqmax")
        error = check_attribute(filename, options%nsqmin, "nsqmin")
        error = check_attribute(filename, options%n_dir_values, "n_dir_values")
        error = check_attribute(filename, options%n_nsq_values, "n_nsq_values")
        error = check_attribute(filename, options%n_spd_values, "n_spd_values")
        
        ! return error
    end function read_LUT

    function dims_match(filename,dimensions)
        implicit none
        character(len=*), intent(in) :: filename
        integer, dimension(3,2), intent(in) :: dimensions
        integer, dimension(io_maxDims) :: dims
        logical :: dims_match
        integer :: i
        
        dims_match=.True.
        call io_getdims(filename,"uLUT",dims)
        do i=1,3
            if (dimensions(i,1)/=dims(i+1)) then
                dims_match=.False.
                return
            endif
        enddo

        call io_getdims(filename,"vLUT",dims)
        do i=1,3
            if (dimensions(i,2)/=dims(i+1)) then
                dims_match=.False.
                return
            endif
        enddo

    end function dims_match
        
        

    function write_var(filename,varname,data_out, dimnames, open_new_file) result(error)
            implicit none
            ! We are writing 6D data, a nx, nz, ny, na, nb, nc grid. 
            integer, parameter :: ndims = 6
            
            ! This is the name of the file and variable we will write. 
            character(len=*), intent(in) :: filename, varname
            real,intent(in) :: data_out(:,:,:,:,:,:)
            character(len=*), optional, dimension(ndims), intent(in) :: dimnames
            logical, optional, intent(in) :: open_new_file
            integer :: error, i
            
            ! This will be the netCDF ID for the file and data variable.
            integer :: ncid, varid,dimids(ndims)
            character(len=MAXVARLENGTH), dimension(ndims) :: dims
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
            if (error==0) then
                error = error + 2**(ndims+1) * setup_var(ncid, varname, varid, dimids, shape(data_out))
            endif
            
            ! End define mode. This tells netCDF we are done defining metadata.
            call check( nf90_enddef(ncid) )
            
            ! write the actual data to the file
            if (varid /= -1) then
                call check( nf90_put_var(ncid, varid, data_out), trim(filename)//":"//trim(varname))
            endif
            
            ! Close the file, freeing all resources.
            call check( nf90_close(ncid), filename)
            
            ! return error
    end function write_var


    function setup_var(ncid, varname, varid, dimids, dims) result(error)
        implicit none
        
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: varname
        integer, intent(out) :: varid
        integer, intent(in), dimension(:) :: dimids
        integer, intent(in), dimension(:) :: dims
        integer :: error
        
        integer, dimension(:), allocatable :: var_dimids
        integer :: i, dimlen, ndims, var_exist_error
        
        ndims = size(dimids)
        
        var_exist_error = nf90_inq_varid(ncid, varname, varid)
        ! if the variable doesn't exist, create it
        if (var_exist_error/=NF90_NOERR) then
            call check( nf90_def_var(ncid, varname, NF90_REAL, dimids, varid), "Creating var:"//trim(varname))
        endif
        
        ! ! find the number of dimensions
        ! call check(nf90_inquire_variable(ncid, varid, ndims = ndims),varname)
        allocate(var_dimids(ndims))
        call check(nf90_inquire_variable(ncid, varid, dimids = var_dimids),   trim(varname)//" dims")
        do i=1,ndims
            call check(nf90_inquire_dimension(ncid, var_dimids(i), len = dimlen), "inq "//trim(varname)//" dim:"//trim(str(i)))
            if (dimlen/=dims(i)) then
                error=3
            endif
        end do
        deallocate(var_dimids)
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

    function check_real_attribute(filename, default_value, value_name) result(error)
        implicit none
        character(len=*), intent(in) :: filename
        real, intent(in) :: default_value
        ! real, intent(inout) :: default_value
        character(len=*), intent(in) :: value_name
        integer :: error
        real :: test_value
        
        ! default return value assumes no error
        error = 0
        
        call io_read_attribute(filename, value_name, test_value)
        if (abs(default_value - test_value)>1e-20) then
            error = 1
            write(*,*) "WARNING parameter option: "//trim(value_name)//"  "//trim(str(default_value)), &
                       " did not match file value: ",trim(str(test_value))
            ! default_value = test_value
        endif
        ! returns any error
    end function check_real_attribute
    
    function check_integer_attribute(filename, default_value, value_name) result(error)
        implicit none
        character(len=*), intent(in) :: filename
        integer, intent(in) :: default_value
        ! integer, intent(inout) :: default_value
        character(len=*), intent(in) :: value_name
        integer :: error
        integer :: test_value
        
        ! default return value assumes no error
        error = 0
        
        call io_read_attribute(filename, value_name, test_value)
        if (abs(default_value - test_value)>1e-20) then
            error = 1
            write(*,*) "WARNING parameter option: "//trim(value_name)//"  "//trim(str(default_value)), &
                       " did not match file value: ",trim(str(test_value))
            ! default_value = test_value
        endif
        ! returns any error
    end function check_integer_attribute

end module linear_theory_lut_disk_io
