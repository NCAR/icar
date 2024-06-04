!>----------------------------------------------------------
!!  Define the interface for the output object
!!
!!  Output objects store all of the data and references to data necessary to write
!!  an output file.  This includes primarily internal netcdf related IDs.
!!  Output objects also store an array of variables to output.
!!  These variables maintain pointers to the data to be output as well as
!!  Metadata (e.g. dimension names, units, other attributes)
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module output_interface
  use netcdf

  use icar_constants
  use variable_interface, only : variable_t
  use domain_interface,   only : domain_t
  use meta_data_interface,only : meta_data_t
  use options_interface,  only : options_t
  use time_object,        only : Time_type, THREESIXTY, GREGORIAN, NOCALENDAR, NOLEAP

  implicit none

  private
  public :: output_t

  !>----------------------------------------------------------
  !! Output type definition
  !!
  !!----------------------------------------------------------
  type, extends(meta_data_t) :: output_t

      ! all components are private and should only be modified through procedures
      private

      ! Store the variables to be written
      ! Note n_variables may be smaller then size(variables) so that it doesn't
      ! have to keep reallocating variables whenever something is added or removed
      integer, public :: n_variables = 0
      type(variable_t), public, allocatable :: variables(:)
      ! time variable , publicis stored outside of the variable list... probably need to think about that some
      type(variable_t) :: time

      ! store status of the object
      logical :: is_initialized = .false.
      logical :: creating = .false.

      ! The filename of the netcdf file to write
      character(len=kMAX_FILE_LENGTH) :: filename

      ! the netcdf ID for an open file
      integer :: ncfile_id

      ! number of dimensions in the file
      integer :: n_dims = 0
      ! list of netcdf dimension IDs
      integer :: dim_ids(kMAX_DIMENSIONS)
      ! name of the dimensions in the file
      character(len=kMAX_DIM_LENGTH) :: dimensions(kMAX_DIMENSIONS)
      character(len=kMAX_NAME_LENGTH) :: time_units
      ! global attribute restarted_from
      ! values: 2000-01-01_00:00:00 or "Not Restarted"
      character(len=20) :: restarted_from

  contains

      procedure, public  :: init
      procedure, public  :: add_to_output
      procedure, public  :: add_variables
      procedure, public  :: save_file
      procedure, public  :: set_domain
      procedure, public  :: set_restart_variable

      procedure, private :: init_variables
      procedure, private :: increase_var_capacity
  end type

  interface

      !>----------------------------------------------------------
      !! Initialize the object: set domain, add variables, set restart variable
      !!
      !!----------------------------------------------------------
      module subroutine init(this, domain, options, file_date_format)
          implicit none
          class(output_t),   intent(inout)  :: this
          type(domain_t),    intent(in)     :: domain
          type(options_t),   intent(in)     :: options
          character(len=49), intent(in)     :: file_date_format
      end subroutine

      !>----------------------------------------------------------
      !! Initialize the object (e.g. allocate the variables array)
      !!
      !!----------------------------------------------------------
      module subroutine init_variables(this)
          implicit none
          class(output_t),   intent(inout)  :: this
      end subroutine

      !>----------------------------------------------------------
      !! Increase the size of the variables array if necessary
      !!
      !!----------------------------------------------------------
      module subroutine increase_var_capacity(this)
          implicit none
          class(output_t),   intent(inout)  :: this
      end subroutine

      !>----------------------------------------------------------
      !! Set the domain data structure to be used when writing
      !!
      !!----------------------------------------------------------
      module subroutine set_domain(this, domain)
          implicit none
          class(output_t),  intent(inout)  :: this
          type(domain_t),   intent(in)     :: domain
      end subroutine

      !>----------------------------------------------------------
      !! Set the restart variable to be used when writing global
      !! attributes
      !!----------------------------------------------------------
      module subroutine set_restart_variable(this, restart_time)
          implicit none
          class(output_t),  intent(inout)  :: this
          character(len=25), intent(in) :: restart_time
      end subroutine

      !>----------------------------------------------------------
      !! Add a variable to the list of output variables
      !!
      !!----------------------------------------------------------
      module subroutine add_to_output(this, variable)
          implicit none
          class(output_t),   intent(inout)  :: this
          type(variable_t),  intent(in)     :: variable
      end subroutine

      !>----------------------------------------------------------
      !! Add multiple variables to the list of output variables given their integer constant index
      !!
      !!----------------------------------------------------------
      module subroutine add_variables(this, var_list, domain)
          implicit none
          class(output_t), intent(inout)  :: this
          integer,         intent(in)     :: var_list(:)
          type(domain_t),  intent(in)     :: domain
      end subroutine

      !>----------------------------------------------------------
      !! Save a new timestep (time) to the file "filename"
      !!
      !!----------------------------------------------------------
      module subroutine save_file(this, filename, current_step, time)
          implicit none
          class(output_t),  intent(inout) :: this
          character(len=*), intent(in)    :: filename
          integer,          intent(in)    :: current_step
          type(Time_type),  intent(in)    :: time
      end subroutine

  end interface
end module
