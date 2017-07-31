module debug_module
    use data_structures, only   : domain_type
    use string,          only   : str

    implicit none
contains

    subroutine domain_check(domain, error_msg)
        implicit none
        type(domain_type),  intent(inout)   :: domain
        character(len=*),   intent(in)      :: error_msg
        
        call check_var(domain%th,       name="th",      msg=error_msg, less_than    =100.0, fix=.True.)
        call check_var(domain%th,       name="th",      msg=error_msg, greater_than =600.0, fix=.True.)
        call check_var(domain%qv,       name="qv",      msg=error_msg, less_than    =0.0, fix=.True.)
        call check_var(domain%cloud,    name="cloud",   msg=error_msg, less_than    =0.0, fix=.True.)
        call check_var(domain%ice,      name="ice",     msg=error_msg, less_than    =0.0, fix=.True.)
        call check_var(domain%nice,     name="nice",    msg=error_msg, less_than    =0.0, fix=.True.)
        call check_var(domain%qsnow,    name="qsnow",   msg=error_msg, less_than    =0.0, fix=.True.)
        call check_var(domain%nsnow,    name="nsnow",   msg=error_msg, less_than    =0.0, fix=.True.)
        call check_var(domain%qrain,    name="qrain",   msg=error_msg, less_than    =0.0, fix=.True.)
        call check_var(domain%nrain,    name="nrain",   msg=error_msg, less_than    =0.0, fix=.True.)
        call check_var(domain%qgrau,    name="qgrau",   msg=error_msg, less_than    =0.0, fix=.True.)
        call check_var(domain%ngraupel, name="ngraupel",msg=error_msg, less_than    =0.0, fix=.True.)
        
    end subroutine domain_check
    
    subroutine check_var(var, name, msg, greater_than, less_than, fix)
        implicit none
        real,               intent(inout), dimension(:,:,:), allocatable :: var
        character(len=*),   intent(in)                      :: name, msg
        real,               intent(in),    optional         :: greater_than, less_than
        logical,            intent(in),    optional         :: fix
        real :: vmax, vmin
        
        if (.not.allocated(var)) then
            return
        endif
        
        if (present(greater_than)) then
            vmax = maxval(var)
            if (vmax > greater_than) then
                write(*,*) trim(msg)
                write(*,*) trim(name)//" is greater than "//trim(str(greater_than))//" : "//trim(str(vmax))
                
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
                write(*,*) trim(name)//" is less than "//trim(str(less_than))//" : "//trim(str(vmin))
                
                if (present(fix)) then
                    if (fix) then
                        write(*,*) "Fixing..."
                        where(var < less_than) var = less_than
                    endif
                endif
            endif
        endif
        
    end subroutine check_var

    subroutine domain_fix(domain)
        implicit none
        type(domain_type),  intent(inout)   :: domain
        
        call fix_var(domain%th,     less_than    =100.0)
        call fix_var(domain%th,     greater_than =600.0)
        call fix_var(domain%qv,     less_than    =0.0)
        call fix_var(domain%cloud,  less_than    =0.0)
        call fix_var(domain%ice,    less_than    =0.0)
        call fix_var(domain%nice,   less_than    =0.0)
        call fix_var(domain%qsnow,  less_than    =0.0)
        call fix_var(domain%nsnow,  less_than    =0.0)
        call fix_var(domain%qrain,  less_than    =0.0)
        call fix_var(domain%nrain,  less_than    =0.0)
        call fix_var(domain%qgrau,  less_than    =0.0)
        call fix_var(domain%ngraupel,less_than   =0.0)

    end subroutine domain_fix
    
    subroutine fix_var(var, greater_than, less_than)
        implicit none
        real,               intent(inout), dimension(:,:,:), allocatable :: var
        real,               intent(in),    optional         :: greater_than, less_than
        
        if (.not.allocated(var)) then
            return
        endif
        
        if (present(greater_than)) then
            where(var > greater_than) var = greater_than
        endif

        if (present(less_than)) then
            where(var < less_than) var = less_than
        endif
        
    end subroutine fix_var

    
end module debug_module
