submodule(meta_data_interface) meta_data_implementation
  implicit none

contains

    ! add a given attribute name/value pair to the metadata object
    module subroutine add_attribute(this, input_name, input_value)
        class(meta_data_t), intent(inout) :: this
        character(len=*),   intent(in)    :: input_name
        character(len=*),   intent(in)    :: input_value

        if (.not.allocated(this%attributes))       call allocate_attributes(this)
        if (this%n_attrs == size(this%attributes)) call increase_holding_capacity(this)

        this%n_attrs = this%n_attrs+1
        this%attributes( this%n_attrs)%name = input_name
        this%attributes(this%n_attrs)%value = input_value

    end subroutine

    ! allocate attribute arrays
    subroutine allocate_attributes(this)
        implicit none
        class(meta_data_t), intent(inout) :: this

        ! note, this is just to initialize the arrays. They will grow dynamically if more than
        ! 16 attributes are added to a metadata object using increase_holding_capacity
        if (.not.allocated(this%attributes))   allocate(this%attributes(16))

    end subroutine allocate_attributes

    ! increase the size of the attribute arrays
    subroutine increase_holding_capacity(this)
        implicit none
        class(meta_data_t),   intent(inout)  :: this
        type(attribute_t), allocatable :: attributes(:)

        ! assert allocated(this%attributes)
        allocate(attributes,  source=this%attributes)

        deallocate(this%attributes)

        allocate(this%attributes(size(attributes)*2))
        this%attributes(:size(attributes)) = attributes

    end subroutine

end submodule
