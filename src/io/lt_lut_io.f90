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
    use io_routines,    only: file_exists, io_read, check, &
                              io_read_attribute, io_add_attribute, &
                              io_maxDims, io_getdims
    use options_types,  only: lt_options_type
    implicit none

    private
    public :: read_lut
    public :: write_lut

    character(len=10), parameter :: lt_lut_version="1.1"

    interface write_var
        module procedure write_var_1d,write_var_6d
    end interface write_var

    interface check_attribute
        module procedure check_attribute_r,check_attribute_i,check_attribute_c
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
    function write_lut(filename, uLUT, vLUT, dz, options) result(error)
        implicit none
        character(len=*), intent(in) :: filename
        real, dimension(:,:,:,:,:,:), intent(in) :: uLUT, vLUT
        real, dimension(:),           intent(in) :: dz
        type(lt_options_type),        intent(in) :: options
        integer :: error
        error=0

        if (file_exists(filename)) then
            if (this_image()==1) write(*,*) "WARNING: Linear Theory look-up-table file already exists."
            error = 1
            if (.not.options%overwrite_lt_lut) then
                if (this_image()==1) write(*,*) "WARNING: Existing file will not be overwritten."
                return
            endif
            if (this_image()==1) write(*,*) "WARNING: Existing file *will* be overwritten."
        endif
        ! allocate(rev_u_LUT(n_spd_values,n_dir_values,n_nsq_values,nxu,nz,ny))
        error = write_var(filename,"uLUT",uLUT, dimnames=[character(len=4) :: "nspd","ndir","nnsq","nxu","nz","ny"], open_new_file=.True.)
        if (error/=0) then
            if (this_image()==1) write(*,*) "Error writing uLUT to file:"//trim(filename)//" Error code = "//trim(str(error))
            return
        endif
        error = write_var(filename,"vLUT",vLUT, dimnames=[character(len=4) :: "nspd","ndir","nnsq","nx","nz","nyv"])
        if (error/=0) then
            if (this_image()==1) write(*,*) "Error writing vLUT to file:"//trim(filename)//" Error code = "//trim(str(error))
            return
        endif

        error = write_var(filename,"dz",dz, dimnames=[character(len=4) :: "nz"])
        if (error/=0) then
            if (this_image()==1) write(*,*) "Error writing dz to LUT file:"//trim(filename)//" Error code = "//trim(str(error))
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
        call io_add_attribute(filename, "minimum_layer_size", options%minimum_layer_size)

        call io_add_attribute(filename, "lt_LUT_version", trim(lt_lut_version))

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
    !!  ... or at least they would if options was an inout variable... for now no update is performed an error is
    !!  returned so that the requested LUT can be computed instead.
    !!
    !!  @param filename Name of LUT file to write
    !!  @param uLUT     Array to hold u wind look up table
    !!  @param vLUT     Array to hold v wind look up table
    !!  @param dims     1D array containtin dimension sizes of the model domain
    !!  @param options  Linear wind options to compare with the LUT file and update if necessary
    !!  @retval error   0 if there were no errors,
    !!                  1 if parameters did not match,
    !!                  2 if dimensions don't match,
    !!                  3 if file doesn't exist,
    !!                  4 if the version in the LUT file doesn't match that of the code
    !!
    !!----------------------------------------------------------
    function read_LUT(filename, uLUT, vLUT, dz, dims, options) result(error)
        implicit none
        character(len=*), intent(in) :: filename
        real, allocatable, dimension(:,:,:,:,:,:), intent(inout) :: uLUT, vLUT
        real, dimension(:), intent(in) :: dz
        integer, dimension(3,2), intent(in) :: dims
        type(lt_options_type), intent(in) :: options
        ! type(lt_options_type), intent(inout) :: options
        integer :: error
        error=0

        if (.not.file_exists(filename)) then
            error = 3
            return
        endif
        error = check_attribute(filename,lt_lut_version, "lt_LUT_version")
        if (error/=0) then
            error=4
            if (this_image()==1) write(*,*) "WARNING: LUT not read, file lt_LUT_version does not match code"
            return
        endif
        if (.not.dims_match(filename, dims)) then
            error = 2
            if (this_image()==1) write(*,*) "WARNING: LUT not read, LUT dims and specified dims do not match"
            return
        endif

        ! read in the parameter values and test them against the input parameters
        ! if the file and input parameters do not match, adjust the input parameter to match
        error = error + check_attribute(filename, options%dirmax, "dirmax")
        error = error + check_attribute(filename, options%dirmin, "dirmin")
        error = error + check_attribute(filename, options%spdmax, "spdmax")
        error = error + check_attribute(filename, options%spdmin, "spdmin")
        error = error + check_attribute(filename, options%nsqmax, "nsqmax")
        error = error + check_attribute(filename, options%nsqmin, "nsqmin")
        error = error + check_attribute(filename, options%n_dir_values, "n_dir_values")
        error = error + check_attribute(filename, options%n_nsq_values, "n_nsq_values")
        error = error + check_attribute(filename, options%n_spd_values, "n_spd_values")
        error = error + check_attribute(filename, options%minimum_layer_size, "minimum_layer_size")

        ! check that the model levels are actually in the same vertical positions too.
        error = error + check_dz(filename,dz)

        if (error/=0) return

        ! this is the slow part where we actually read in all the data
        call io_read(filename,"uLUT",uLUT)
        call io_read(filename,"vLUT",vLUT)

        ! return no error
    end function read_LUT

    !>----------------------------------------------------------
    !!  Check that the dimensions of a LUT file and current expected LUT match
    !!
    !!  @detail
    !!  Takes a filename and the dimensions of the LUT the model is expecting
    !!  Checks the expected dimensions against the dimension lengths for both
    !!  the vLUT and the uLUT in the specified file.
    !!
    !!  @param filename     Name of LUT file to check
    !!  @param dimensions   Array to hold u and v dimensions array is (3 x 2)
    !!  @retval dims_match  True if they match, false if not
    !!
    !!----------------------------------------------------------
    function dims_match(filename,dimensions)
        implicit none
        character(len=*), intent(in) :: filename
        integer, dimension(3,2), intent(in) :: dimensions
        integer, dimension(io_maxDims) :: dims
        logical :: dims_match
        integer :: i

        ! default value = true, set to false if they don't match
        dims_match=.True.
        ! read the dimensions for the uLUT first
        call io_getdims(filename,"uLUT",dims)
        ! check all three spatial dimensions (nx, nz, ny)
        do i=1,3
            ! note dims(i+4) because dims(1)=ndims, dims(2,3,4)=nspd,ndir,nnsq
            if (dimensions(i,1)/=dims(i+4)) then
                dims_match=.False.
                return
            endif
        enddo

        ! next perform the same check for the vLUT
        call io_getdims(filename,"vLUT",dims)
        do i=1,3
            if (dimensions(i,2)/=dims(i+4)) then
                dims_match=.False.
                return
            endif
        enddo
        ! if we get here we are returning True
    end function dims_match

    !>----------------------------------------------------------
    !!  Check that the dz level thicknesses of a LUT file and current domain match
    !!
    !!  @detail
    !!  Takes a filename and the dz thicknesses of the current model run.
    !!  Checks the expected dz against the LUT file dz.
    !!
    !!  @param filename     Name of LUT file to check
    !!  @param model_dz     Array to hold current model thicknesses (1D real)
    !!  @retval error       0 if they match, 1 otherwise
    !!
    !!----------------------------------------------------------
    function check_dz(filename, model_dz) result(error)
        implicit none
        character(len=*), intent(in) :: filename ! LUT netcdf filename
        real, intent(in), dimension(:) :: model_dz ! dz being used in the model to verify matches LUT file
        integer :: error ! return value

        ! local variables
        real, allocatable, dimension(:) :: lut_dz ! variable to hold the values read from the LUT file
        integer :: i, nz

        ! default assumes no error, code will be set if an error occurs
        error=0

        ! read in the dz variable from the lut file
        call io_read(filename,"dz",lut_dz)
        nz = size(lut_dz)

        ! I don't think this should be possible because we already checked dims_match for the file
        if (nz/=size(model_dz)) then
            error=1
            return
        endif

        ! loop through all levels verifying that the LUT file matches the current model values
        do i=1,nz
            if (model_dz(i)/=lut_dz(i)) then
                error=1
                return
            endif
        enddo
        ! returning 0 if we get to here
    end function check_dz


    !>----------------------------------------------------------
    !!  Write a 6D variable to a netCDF file
    !!
    !!  @detail
    !!  Takes a filename variable name, data, dimension names and optionally a
    !!  flag to clobber existing files. If the specified dimensions don't exist
    !!  they are created based on the size if the data_out variable. Then the
    !!  netcdf variable is created with those dimensions and the data are written.
    !!
    !!  @param filename     Name of LUT file to open or create
    !!  @param varname      Name of variable to create
    !!  @param data_out     Array with the data to be output to the file (6D real)
    !!  @param dimnames     Array of dimension names (6 element 1D character)
    !!  @param open_new_file Flag to clobber any existing file if true
    !!  @retval error       0 if the file was successfully written
    !!
    !!----------------------------------------------------------
    function write_var_6d(filename,varname,data_out, dimnames, open_new_file) result(error)
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

            error=0
            ncid=-1
            varid=-1
            ! Open the file. NF90_NOCLOBBER tells netCDF we want append to existing files
            !   NF90_NETCDF4 tells netCDF that this may be a really LARGE file and it needs to support >4GB/variable
            if (present(open_new_file)) then
                if (open_new_file) then
                    call check( nf90_create(filename, or(NF90_CLOBBER,NF90_NETCDF4), ncid), filename)
                endif
            endif
            if (ncid==-1) then
                ! open an existing netcdf file
                call check( nf90_open(filename, NF90_WRITE, ncid), "write_var:creating file:"//filename )
                ! put file back into define mode so new variables and dimensions can be created.
                call check( nf90_redef(ncid) )
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
            call check( nf90_enddef(ncid), "enddef file:"//trim(filename) )
            ! write the actual data to the file
            if (error == 0) then
                call check( nf90_put_var(ncid, varid, data_out), trim(filename)//":"//trim(varname))
            endif

            ! Close the file, freeing all resources.
            call check( nf90_close(ncid), filename)

            ! return error
    end function write_var_6d

    !>----------------------------------------------------------
    !!  Same as write_var_6d but for a 1D variable
    !!
    !!  @detail
    !!  Takes a filename variable name, data, dimension names and optionally a
    !!  flag to clobber existing files. If the specified dimensions don't exist
    !!  they are created based on the size if the data_out variable. Then the
    !!  netcdf variable is created with those dimensions and the data are written.
    !!
    !!  @param filename     Name of LUT file to open or create
    !!  @param varname      Name of variable to create
    !!  @param data_out     Array with the data to be output to the file (1D real)
    !!  @param dimnames     Array of dimension names (1 element 1D character)
    !!  @param open_new_file Flag to clobber any existing file if true
    !!  @retval error       0 if the file was successfully written
    !!
    !!----------------------------------------------------------
    function write_var_1d(filename,varname,data_out, dimnames, open_new_file) result(error)
            implicit none
            ! We are writing 1D data expected to be nz grid.
            integer, parameter :: ndims = 1

            ! This is the name of the file and variable we will write.
            character(len=*), intent(in) :: filename, varname
            real,intent(in) :: data_out(:)
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
                dims = ["z"]
            endif

            error=0
            ncid=-1
            varid=-1
            ! Open the file. NF90_NOCLOBBER tells netCDF we want append to existing files
            ! NF90_NETCDF4 is necessary to support variables >4GB in size
            if (present(open_new_file)) then
                if (open_new_file) then
                    call check( nf90_create(filename, or(NF90_CLOBBER,NF90_NETCDF4), ncid), filename)
                endif
            endif
            if (ncid==-1) then
                ! open an existing netcdf file
                call check( nf90_open(filename, NF90_WRITE, ncid), "write_var:creating file:"//filename )
                ! put file back into define mode so new variables and dimensions can be created.
                call check( nf90_redef(ncid) )
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
            if (error == 0) then
                call check( nf90_put_var(ncid, varid, data_out), trim(filename)//":"//trim(varname))
            endif

            ! Close the file, freeing all resources.
            call check( nf90_close(ncid), filename)

            ! return error
    end function write_var_1d


    !>----------------------------------------------------------
    !!  Setup a new variable to be written
    !!
    !!  @detail
    !!  Takes an open netCDF file id, variable name, dimension ids, and dimension lengths.
    !!  If the variable doesn't exist it is created with the matching dimension ids.
    !!  These dimension lengths are then checked against the dimension lengths supplied.
    !!  This check is performed in case a variable with this name but different dims already
    !!  existed in the file.
    !!  Supplied netCDF if must be a file that is in define mode
    !!
    !!  @param      ncid    integer open netcdf file id
    !!  @param      varname name of the variable to be created
    !!  @param[out] varid   id of the existing or newly created variable
    !!  @param      dimids  array of dimension ids to be used when creating the variable
    !!  @param      dims    array of dimension lengths to be checked for existing variables
    !!  @retval     error   0 if successful 3 otherwise
    !!
    !!----------------------------------------------------------
    function setup_var(ncid, varname, varid, dimids, dims) result(error)
        implicit none
        integer,                intent(in) :: ncid
        character(len=*),       intent(in) :: varname
        integer,                intent(out):: varid
        integer, dimension(:),  intent(in) :: dimids
        integer, dimension(:),  intent(in) :: dims
        integer :: error

        integer, dimension(:), allocatable :: var_dimids
        integer :: i, dimlen, ndims, var_exist_error

        ! default return value=0
        error=0
        ! number of dimensions
        ndims = size(dimids)
        allocate(var_dimids(ndims))

        var_exist_error = nf90_inq_varid(ncid, varname, varid)
        ! if the variable doesn't exist, create it
        if (var_exist_error/=NF90_NOERR) then
            call check( nf90_def_var(ncid, varname, NF90_REAL, dimids, varid), "Creating var:"//trim(varname))
        endif

        call check(nf90_inquire_variable(ncid, varid, dimids = var_dimids),   trim(varname)//" dims")
        do i=1,ndims
            call check(nf90_inquire_dimension(ncid, var_dimids(i), len = dimlen), "inq "//trim(varname)//" dim:"//trim(str(i)))
            if (dimlen/=dims(i)) then
                error=3
            endif
        end do
        deallocate(var_dimids)
    end function setup_var

    !>----------------------------------------------------------
    !!  Setup a new dimension in a given netcdf file
    !!
    !!  Takes an open netcdf file id and a dimension name.
    !!  returns the dimension id for the given dimension name and
    !!  if given a dimension length (n) it is checked to make sure
    !!  any existing dimension with this name has the write length.
    !!
    !!  @param      ncid    integer an open netcdf file id
    !!  @param      dimname the name of the dimension to be checked
    !!  @param[out] dimid   the existing or newly created dimension id.
    !!  @param      n       The length of the dimension
    !!  @retval     error   0 if successful, 1 otherwise
    !!
    !!----------------------------------------------------------
    function setup_dim(ncid, dimname, dimid, n) result(error)
        implicit none
        integer,          intent(in)    :: ncid
        character(len=*), intent(in)    :: dimname
        integer,          intent(inout) :: dimid
        integer,          intent(in)    :: n
        integer :: error
        integer :: ncerror, dim_length

        error = 0

        ! first check for an existing dimension with the given name
        ncerror = NF90_INQ_DIMID(ncid, dimname, dimid)
        ! if the dimension does not exist, create it here.
        if (ncerror/=NF90_NOERR) then
            call check( nf90_def_dim(ncid, dimname, n, dimid), "setup_dim:Creating dim:"//dimname )
        else
            ! else, check to make sure the existing dimension has the correct length
            call check( nf90_inquire_dimension(ncid, dimid, len = dim_length), "setup_dim:Reading dim:"//dimname)

            if (dim_length/=n) then
                error = 1
                if (this_image()==1) write(*,*) "Dim setup error: ", trim(dimname), dimid, dim_length,"!=",n
            endif
        endif
        ! return error
    end function setup_dim

    !>----------------------------------------------------------
    !!  Check that a named real attribute matches the value in a file
    !!
    !!  Takes a filename, attribute name, and expected attribute value
    !!  Reads the attribute from the file, and compares the results to
    !!  the expected value.  If the results are within 1e-20 of each other.
    !!  they are assumed to match, and 0 is returned, else 1 is returned
    !!
    !!  @param  filename        name of the netCDF file
    !!  @param  default_value   expected value of the attribute
    !!  @param  value_name      name of the attribute to read
    !!  @retval error           0 if expected matches read, 1 otherwise
    !!
    !!----------------------------------------------------------
    function check_attribute_r(filename, default_value, value_name) result(error)
        implicit none
        character(len=*), intent(in) :: filename
        real, intent(in) :: default_value
        ! real, intent(inout) :: default_value
        character(len=*), intent(in) :: value_name
        integer :: error
        real :: test_value

        ! default return value assumes no error
        error = 0

        call io_read_attribute(filename, value_name, test_value, error=error)
        if (error/=NF90_NOERR) then
            error=1
            if (this_image()==1) write(*,*) "WARNING: Error reading attribute: ",trim(value_name)," from: ",trim(filename)
        else
            if (default_value/=test_value) then
                error = 1
                if (this_image()==1) write(*,*) "WARNING parameter option: "//trim(value_name)//"  "//trim(str(default_value)), &
                           " did not match file value: ",trim(str(test_value))
            endif
        endif
        ! returns any error
    end function check_attribute_r

    !>----------------------------------------------------------
    !!  Check that a named integer attribute matches the value in a file
    !!
    !!  Takes a filename, attribute name, and expected attribute value
    !!  Reads the attribute from the file, and compares the results to
    !!  the expected value.  If the results match the function is successful
    !!  and 0 is returned, else 1 is returned
    !!
    !!  @param  filename        name of the netCDF file
    !!  @param  default_value   expected value of the attribute
    !!  @param  value_name      name of the attribute to read
    !!  @retval error           0 if expected matches read, 1 otherwise
    !!
    !!----------------------------------------------------------
    function check_attribute_i(filename, default_value, value_name) result(error)
        implicit none
        character(len=*), intent(in) :: filename
        integer, intent(in) :: default_value
        ! integer, intent(inout) :: default_value
        character(len=*), intent(in) :: value_name
        integer :: error
        integer :: test_value

        ! default return value assumes no error
        error = 0

        call io_read_attribute(filename, value_name, test_value, error=error)
        if (error/=NF90_NOERR) then
            error=1
            if (this_image()==1) write(*,*) "WARNING: Error reading attribute: ",trim(value_name)," from: ",trim(filename)
        else
            if (default_value/=test_value) then
                error = 1
                if (this_image()==1) write(*,*) "WARNING parameter option: "//trim(value_name)//"  "//trim(str(default_value)), &
                           " did not match file value: ",trim(str(test_value))
            endif
        endif
        ! returns any error
    end function check_attribute_i

    !>----------------------------------------------------------
    !!  Check that a named character attribute matches the value in a file
    !!
    !!  Takes a filename, attribute name, and expected attribute value
    !!  Reads the attribute from the file, and compares the results to
    !!  the expected value.  If the results match the function is successful
    !!  and 0 is returned, else 1 is returned
    !!
    !!  @param  filename        name of the netCDF file
    !!  @param  default_value   expected value of the attribute
    !!  @param  value_name      name of the attribute to read
    !!  @retval error           0 if expected matches read, 1 otherwise
    !!
    !!----------------------------------------------------------
    function check_attribute_c(filename, default_value, value_name) result(error)
        implicit none
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: default_value
        character(len=*), intent(in) :: value_name

        integer :: error
        character(len=MAXVARLENGTH)  :: test_value

        ! default return value assumes no error
        error = 0

        call io_read_attribute(filename, value_name, test_value, error=error)
        if (error/=NF90_NOERR) then
            error=1
            if (this_image()==1) write(*,*) "WARNING: Error reading attribute: ",trim(value_name)," from: ",trim(filename)
        else
            if (default_value/=test_value) then
                error = 1
                if (this_image()==1) write(*,*) "WARNING parameter option: "//trim(value_name)//"  "//trim(default_value), &
                           " did not match file value: ",trim(test_value)
            endif
        endif
        ! returns any error
    end function check_attribute_c

end module linear_theory_lut_disk_io
