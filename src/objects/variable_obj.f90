submodule(variable_interface) variable_implementation
  implicit none


contains

    module function new_variable(name)
        character(len=kMAX_NAME_LENGTH), intent(in) :: name
        type(variable_t) :: new_variable

        new_variable%name = name

        new_variable%var_id = -9999

    end function new_variable

end submodule
