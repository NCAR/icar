submodule(meta_data_interface) meta_data_implementation
  implicit none

contains

    ! add a given attribute name/value pair to the metadata object
    module subroutine add_attribute(this, input_name, input_value)
        class(meta_data_t), intent(inout) :: this
        character(len=*),   intent(in)    :: input_name
        character(len=*),   intent(in)    :: input_value

        if (.not.allocated(this%attribute_values))       call allocate_attributes(this)
        if (this%n_attrs == size(this%attribute_values)) call increase_holding_capacity(this)

        this%n_attrs = this%n_attrs+1
        this%attribute_names( this%n_attrs) = input_name
        this%attribute_values(this%n_attrs) = input_value

    end subroutine

    ! allocate attribute arrays
    subroutine allocate_attributes(this)
        implicit none
        class(meta_data_t), intent(inout) :: this

        ! note, this is just to initialize the arrays. They will grow dynamically if more than
        ! 16 attributes are added to a metadata object using increase_holding_capacity
        if (.not.allocated(this%attribute_names))   allocate(this%attribute_names(16))
        if (.not.allocated(this%attribute_values))  allocate(this%attribute_values(16))

    end subroutine allocate_attributes

    ! increase the size of the attribute arrays
    subroutine increase_holding_capacity(this)
        implicit none
        class(meta_data_t),   intent(inout)  :: this
        character(len=kMAX_ATTR_LENGTH), allocatable :: attribute_names(:)
        character(len=kMAX_ATTR_LENGTH), allocatable :: attribute_values(:)

        ! assert allocated(this%attribute_names)
        ! assert allocated(this%attribute_values)
        allocate(attribute_names,  source=this%attribute_names)
        allocate(attribute_values, source=this%attribute_values)

        deallocate(this%attribute_names)
        deallocate(this%attribute_values)

        allocate(this%attribute_names(size(attribute_names)*2))
        this%attribute_names(:size(attribute_names)) = attribute_names

        allocate(this%attribute_values(size(attribute_values)*2))
        this%attribute_values(:size(attribute_values)) = attribute_values

    end subroutine

end submodule
