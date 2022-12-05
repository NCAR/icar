module debug_module
    use domain_interface, only  : domain_t
    use string,           only  : str
    use ieee_arithmetic

    implicit none
contains

    subroutine domain_check(domain, error_msg, fix)
        implicit none
        type(domain_t),     intent(inout)   :: domain
        character(len=*),   intent(in)      :: error_msg
        logical,            intent(in), optional :: fix
        logical :: fix_data

        fix_data = .False.
        if (present(fix)) fix_data = fix

        call check_var(domain%potential_temperature%data_3d, name="th",      msg=error_msg, less_than    =100.0, fix=fix_data)
        call check_var(domain%potential_temperature%data_3d, name="th",      msg=error_msg, greater_than =600.0, fix=fix_data)
        call check_var(domain%water_vapor%data_3d,           name="qv",      msg=error_msg, less_than    =-1e-10,fix=fix_data)
        call check_var(domain%cloud_water_mass%data_3d,      name="cloud",   msg=error_msg, less_than    =-1e-10,fix=fix_data)
        call check_var(domain%cloud_ice_mass%data_3d,        name="ice",     msg=error_msg, less_than    =-1e-10,fix=fix_data)
        call check_var(domain%cloud_ice_number%data_3d,      name="nice",    msg=error_msg, less_than    =-1e-1, fix=fix_data)
        call check_var(domain%snow_mass%data_3d,             name="qsnow",   msg=error_msg, less_than    =-1e-10,fix=fix_data)
        call check_var(domain%snow_number%data_3d,           name="nsnow",   msg=error_msg, less_than    =-1e-1, fix=fix_data)
        call check_var(domain%rain_mass%data_3d,             name="qrain",   msg=error_msg, less_than    =-1e-10,fix=fix_data)
        call check_var(domain%rain_number%data_3d,           name="nrain",   msg=error_msg, less_than    =-1e-1, fix=fix_data)
        call check_var(domain%graupel_mass%data_3d,          name="qgrau",   msg=error_msg, less_than    =-1e-10,fix=fix_data)
        call check_var(domain%graupel_number%data_3d,        name="ngrau",   msg=error_msg, less_than    =-1e-1, fix=fix_data)
        call check_var(domain%w%data_3d,                     name="w",       msg=error_msg, less_than    =-1e5,  fix=fix_data)
        call check_var(domain%w%data_3d,                     name="w",       msg=error_msg, greater_than =1e5,   fix=fix_data)
        call check_var2d(domain%sensible_heat%data_2d,       name="hfx",     msg=error_msg) ! check for NaN's only.
        call check_var2d(domain%latent_heat%data_2d,         name="lfx",     msg=error_msg)
        call check_var2d(domain%skin_temperature%data_2d,    name="tskin",   msg=error_msg)
        call check_var2d(domain%roughness_z0%data_2d,        name="z0",      msg=error_msg)
        call check_var2d(domain%surface_pressure%data_2d,    name="psfc",    msg=error_msg)
        ! call check_var2d(domain%ustar,                       name="ustar", msg=error_msg)
        call check_var(domain%exner%data_3d,                 name="pii",     msg=error_msg)
        call check_var(domain%pressure_interface%data_3d,    name="pi",      msg=error_msg)
        call check_var(domain%pressure%data_3d,              name="p",       msg=error_msg)

    end subroutine domain_check


    subroutine check_var(var, name, msg, greater_than, less_than, fix)
        implicit none
        real,               intent(inout), dimension(:,:,:), pointer :: var
        character(len=*),   intent(in)                      :: name, msg
        real,               intent(in),    optional         :: greater_than, less_than
        logical,            intent(in),    optional         :: fix
        integer :: n
        real :: vmax, vmin
        logical :: printed

        printed = .False.

        if (.not.associated(var)) then
            return
        endif

        if (any(ieee_is_nan(var))) then
            n = COUNT(ieee_is_nan(var))
            ! ALLOCATE(IsNanIdx(n))
            ! IsNanIdx = PACK( (/(i,i=1,SIZE(var))/), MASK=IsNan(var) )  ! if someone can get this to work it would be nice to have.
            write(*,*) trim(msg)
            write(*,*) trim(name)//" has", n," NaN(s) "
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

                if ((vmin - less_than) < -1e-10) then
                ! we only want to hard stop if there is a significant difference.
                ! Numerical precision can mean that advecting hydrometeors ends up with -1e-30 type values which we can ignore
                block
                    integer :: i,j,k,nx,ny,nz

                    nx = size(var,1)
                    nz = size(var,2)
                    ny = size(var,3)
                    do j=lbound(var,3),ubound(var,3)
                        do k=lbound(var,2),ubound(var,2)
                            do i=lbound(var,1),ubound(var,1)
                                if (var(i,k,j) < less_than) then

                                    if (.not.printed) then
                                        print*, "First Error was in grid cell:", i,k,j, var(i,k,j)
                                        printed = .True.
                                    endif
                                    if (.not.present(fix)) then
                                        error stop
                                    endif
                                endif
                            enddo
                        enddo
                    enddo
                end block
                endif

                if (present(fix)) then
                    if (fix) then
                        write(*,*) "Fixing..."
                        where(var < less_than) var = less_than
                    endif
                endif
            endif
        endif

    end subroutine check_var

    subroutine check_var2d(var, name, msg, greater_than, less_than, fix)
        implicit none
        real,               intent(inout), dimension(:,:), pointer :: var
        character(len=*),   intent(in)                      :: name, msg
        real,               intent(in),    optional         :: greater_than, less_than
        logical,            intent(in),    optional         :: fix
        integer :: n, i, j
        real :: vmax, vmin
        logical :: printed

        printed = .False.

        if (.not.associated(var)) then
            return
        endif

        if (any(ieee_is_nan(var))) then
            n = COUNT(ieee_is_nan(var))
            ! ALLOCATE(IsNanIdx(n))
            ! IsNanIdx = PACK( (/(i,i=1,SIZE(var))/), MASK=IsNan(var) )  ! if someone can get this to work it would be nice to have.
            write(*,*) trim(msg)
            write(*,*) trim(name)//" has", n," NaN(s) "
            if (n < 9) then
                do j = lbound(var,2), ubound(var,2)
                    do i = lbound(var,1), ubound(var,1)
                        if (ieee_is_nan(var(i,j))) print*, "NaN in ",i,j
                    enddo
                enddo
            else
                print*, "Too many NaNs, stopping on image:", this_image()
                error stop
            endif

        endif
    end subroutine check_var2d

    ! subroutine domain_fix(domain)
    !     implicit none
    !     type(domain_t),  intent(inout)   :: domain
    !
    !     call fix_var(domain%th,     less_than    =100.0)
    !     call fix_var(domain%th,     greater_than =600.0)
    !     call fix_var(domain%qv,     less_than    =0.0)
    !     call fix_var(domain%cloud,  less_than    =0.0)
    !     call fix_var(domain%ice,    less_than    =0.0)
    !     call fix_var(domain%nice,   less_than    =0.0)
    !     call fix_var(domain%qsnow,  less_than    =0.0)
    !     call fix_var(domain%nsnow,  less_than    =0.0)
    !     call fix_var(domain%qrain,  less_than    =0.0)
    !     call fix_var(domain%nrain,  less_than    =0.0)
    !     call fix_var(domain%qgrau,  less_than    =0.0)
    !     call fix_var(domain%ngraupel,less_than   =0.0)
    !
    ! end subroutine domain_fix

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
