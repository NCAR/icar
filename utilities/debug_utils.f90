module debug_module
    use data_structures, only   : domain_type
    use string,          only   : str

    implicit none
contains

    subroutine domain_check(domain, error_msg)
        implicit none
        type(domain_type),  intent(inout)   :: domain
        character(len=*),   intent(in)      :: error_msg
        
        call check_var(domain%th,       name="th",      msg=error_msg, less_than    =100.0)
        call check_var(domain%th,       name="th",      msg=error_msg, greater_than =600.0)
        call check_var(domain%qv,       name="qv",      msg=error_msg, less_than    =0.0, fix=.True.)
        call check_var(domain%cloud,    name="cloud",   msg=error_msg, less_than    =0.0, fix=.True.)
        call check_var(domain%ice,      name="ice",     msg=error_msg, less_than    =0.0, fix=.True.)
        call check_var(domain%nice,     name="nice",    msg=error_msg, less_than    =0.0, fix=.True.)
        call check_var(domain%qsnow,    name="qsnow",   msg=error_msg, less_than    =0.0, fix=.True.)
        call check_var(domain%qrain,    name="qrain",   msg=error_msg, less_than    =0.0, fix=.True.)
        call check_var(domain%nrain,    name="nrain",   msg=error_msg, less_than    =0.0, fix=.True.)
        
    end subroutine domain_check
    
    subroutine check_var(var, name, msg, greater_than, less_than, fix)
        implicit none
        real,               intent(inout), dimension(:,:,:) :: var
        character(len=*),   intent(in)                      :: name, msg
        real,               intent(in),    optional         :: greater_than, less_than
        logical,            intent(in),    optional         :: fix
        real :: vmax, vmin
        
        if (present(greater_than)) then
            vmax = maxval(var)
            if (vmax > greater_than) then
                write(*,*) trim(msg)
                write(*,*) trim(name)//" is greater than "//str(greater_than)//" : "//str(vmax)
                
                if (present(fix)) then
                    if (fix) then
                        write(*,*) "Fixing..."
                        where(var > greater_than) var = greater_than
                    endif
                endif
            endif
        endif

        if (present(less_than)) then
            vmin = minval(var)
            if (vmin < less_than) then
                write(*,*) trim(msg)
                write(*,*) trim(name)//" is less than "//str(less_than)//" : "//str(vmin)
                
                if (present(fix)) then
                    if (fix) then
                        write(*,*) "Fixing..."
                        where(var < less_than) var = less_than
                    endif
                endif
            endif
        endif
        
    end subroutine check_var
    
end module debug_module
