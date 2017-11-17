!>------------------------------------------------------------
!! Module to manage the ICAR wind field, including calls to linear winds
!! importantly it also rotates the wind field into the ICAR grid and
!! balances the U, V, and W fields for "mass" conservation
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module wind
    use linear_theory_winds, only : linear_perturb
    use mod_blocking,        only : add_blocked_flow
    use data_structures
!   use output, only: write_domain
    implicit none
    private
    public::update_winds, balance_uvw, init_winds
    real, parameter::deg2rad=0.017453293 !2*pi/360
contains

    !>------------------------------------------------------------
    !! Forces u,v, and w fields to balance
    !!       du/dx+dv/dy = dw/dz
    !!
    !! Starts by setting w out of the ground=0 then works through layers
    !!
    !!------------------------------------------------------------
    subroutine balance_uvw(domain,options)
        implicit none
        type(domain_type),intent(inout)::domain
        type(options_type),intent(in)::options
        real, allocatable,dimension(:,:) ::du,dv,divergence,rhou,rhov,rhow
        integer ::nx,ny,nz,i
        nx=size(domain%w,1)
        nz=size(domain%w,2)
        ny=size(domain%w,3)

        !------------------------------------------------------------
        ! These could be module level to prevent lots of allocation/deallocation/reallocations
        ! but these are relatively small, and should be allocated efficiently on the stack
        !------------------------------------------------------------
        if (options%advect_density) then
            allocate(rhou(nx-1,ny-2))
            allocate(rhov(nx-2,ny-1))
            allocate(rhow(nx-2,ny-2))
        endif

        allocate(du(nx-2,ny-2))
        allocate(dv(nx-2,ny-2))
        allocate(divergence(nx-2,ny-2))

        ! If this becomes a bottle neck in the code it could be parallelized over y
        ! loop over domain levels
        do i=1,nz
            !------------------------------------------------------------
            ! If we are incorporating density into the advection equation
            ! then it needs to be incorporated when balancing the wind field
            !
            ! Note that the "else" case below does the same thing without density
            ! and is much easier to understand
            !------------------------------------------------------------
            if (options%advect_density) then
                ! first calculate density on the u and v grid points
                rhou = (domain%rho(1:nx-1,i,2:ny-1) + domain%rho(2:nx,i,2:ny-1))/2
                rhov = (domain%rho(2:nx-1,i,1:ny-1) + domain%rho(2:nx-1,i,2:ny))/2

                ! for most grid points interpolate between the current grid and the one above it
                if (i<nz) then
                    rhow = (domain%rho(2:nx-1,i,2:ny-1) + domain%rho(2:nx-1,i+1,2:ny-1))/2
                else
                    ! for the top grid cell extrapolate upwards based on the current grid and the one below it
                    rhow = 2*domain%rho(2:nx-1,i,2:ny-1) - domain%rho(2:nx-1,i-1,2:ny-1)
                endif

                ! calculate horizontal divergence based on the wind * density field
                domain%vr(2:nx-1,i,2:ny) = rhov * domain%v(2:nx-1,i,2:ny) * domain%dz_inter(1,i,1) * domain%dx
                dv = domain%vr(2:nx-1,i,3:ny) - domain%vr(2:nx-1,i,2:ny-1)

                domain%ur(2:nx,i,2:ny-1) = rhou * domain%u(2:nx,i,2:ny-1) * domain%dz_inter(1,i,1) * domain%dx
                du = domain%ur(3:nx,i,2:ny-1) - domain%ur(2:nx-1,i,2:ny-1)

                divergence = du + dv
                if (i==1) then
                    ! if this is the first model level start from 0 at the ground
                    domain%wr(2:nx-1,i,2:ny-1) = 0 - divergence
                    domain%w (2:nx-1,i,2:ny-1) = 0 - divergence / rhow / (domain%dx**2)
                else
                    ! else calculate w as a change from w at the level below
                    domain%wr(2:nx-1,i,2:ny-1) = domain%wr(2:nx-1,i-1,2:ny-1) - divergence
                    domain%w (2:nx-1,i,2:ny-1) = domain%w (2:nx-1,i-1,2:ny-1) - divergence / rhow / (domain%dx**2)
                endif
            else
                !------------------------------------------------------------
                ! If we are not incorporating density this is simpler
                !------------------------------------------------------------
                ! calculate horizontal divergence
                !   in the North-South direction
                dv = domain%v(2:nx-1,i,3:ny) - domain%v(2:nx-1,i,2:ny-1)
                !   in the East-West direction
                du = domain%u(3:nx,i,2:ny-1) - domain%u(2:nx-1,i,2:ny-1)
                !   in net
                divergence = du + dv

                ! Then calculate w to balance
                if (i==1) then
                    ! if this is the first model level start from 0 at the ground
                    domain%w(2:nx-1,i,2:ny-1) = 0 - divergence
                else
                    ! else calculate w as a change from w at the level below
                    domain%w(2:nx-1,i,2:ny-1) = domain%w(2:nx-1,i-1,2:ny-1) - divergence
                endif

                !------------------------------------------------------------
                ! Now do the same for the convective wind field if needed
                !------------------------------------------------------------
                if (options%physics%convection > 0) then
                    ! calculate horizontal divergence
                    dv = domain%v_cu(2:nx-1,i,3:ny) - domain%v_cu(2:nx-1,i,2:ny-1)
                    du = domain%u_cu(3:nx,i,2:ny-1) - domain%u_cu(2:nx-1,i,2:ny-1)
                    divergence = du + dv
                    ! Then calculate w to balance
                    if (i==1) then
                        ! if this is the first model level start from 0 at the ground
                        domain%w_cu(2:nx-1,i,2:ny-1) = 0 - divergence
                    else
                        ! else calculate w as a change from w at the level below
                        domain%w_cu(2:nx-1,i,2:ny-1) = domain%w_cu(2:nx-1,i-1,2:ny-1) - divergence
                    endif
                endif

            endif
        enddo

    end subroutine balance_uvw

    !>------------------------------------------------------------
    !! Correct for a grid that is locally rotated with respect to EW,NS
    !!
    !! Assumes forcing winds are EW, NS relative, not grid relative.
    !!
    !!------------------------------------------------------------
    subroutine make_winds_grid_relative(domain)
        type(domain_type), intent(inout) :: domain
        real,dimension(:),allocatable :: u,v
        integer::nx,nz,ny,k,j

        nx=size(domain%p,1)
        nz=size(domain%p,2)
        ny=size(domain%p,3)

        allocate(u(nx))
        allocate(v(nx))
        !assumes u and v come in on a staggered Arakawa C-grid with one additional grid point in x/y for u/v respectively
        ! destagger to a centered grid (the mass grid)
        domain%u(:nx,:,:) = (domain%u(:nx,:,:)+domain%u(2:,:,:))/2
        domain%v(:,:,:ny) = (domain%v(:,:,:ny)+domain%v(:,:,2:))/2
        do j=1,ny
            do k=1,nz
                ! rotate wind field to the real grid
                u=domain%u(:nx,k,j)*domain%costheta(:nx,j) + domain%v(:nx,k,j)*domain%sintheta(:nx,j)
                v=domain%v(:nx,k,j)*domain%costheta(:nx,j) + domain%u(:nx,k,j)*domain%sintheta(:nx,j)
                domain%u(:nx,k,j)=u
                domain%v(:nx,k,j)=v
            enddo
        enddo
        deallocate(u,v)

        ! put the fields back onto a staggered grid, having effectively lost two grid cells in the staggered directions
        ! estimate the "lost" grid cells by extrapolating beyond the remaining
        domain%u(2:nx,:,:) = (domain%u(1:nx-1,:,:)+domain%u(2:nx,:,:))/2
        domain%u(1,:,:)    = 2*domain%u(1,:,:)  - domain%u(2,:,:)
        domain%u(nx+1,:,:) = 2*domain%u(nx,:,:) - domain%u(nx-1,:,:)

        domain%v(:,:,2:ny) = (domain%v(:,:,1:ny-1)+domain%v(:,:,2:ny))/2
        domain%v(:,:,1)    = 2*domain%v(:,:,1)  - domain%v(:,:,2)
        domain%v(:,:,ny+1) = 2*domain%v(:,:,ny) - domain%v(:,:,ny-1)

    end subroutine make_winds_grid_relative


    !>------------------------------------------------------------
    !! Apply wind field physics and adjustments
    !!
    !! This will call the linear wind module if necessary, otherwise it just updates for
    !! This should ONLY be called once for each forcing step, otherwise effects will be additive.
    !!
    !!------------------------------------------------------------
    subroutine update_winds(domain,options)
        implicit none
        type(domain_type),intent(inout)::domain
        type(options_type),intent(in)::options
        real,allocatable,dimension(:,:,:)::temparray
        integer::nx,ny,nz,i,j

        ! rotate winds from cardinal directions to grid orientation (e.g. u is grid relative not truly E-W)
        if (.not.options%ideal) then
            call make_winds_grid_relative(domain)
        endif

        ! flow blocking parameterization
        if (options%block_options%block_flow) then
            call add_blocked_flow(domain, options)
        endif

        ! linear winds
        if (options%physics%windtype==kWIND_LINEAR) then
            call linear_perturb(domain,options,options%lt_options%vert_smooth,.False.,options%advect_density)
        endif
        ! else assumes even flow over the mountains

        ! use horizontal divergence (convergence) to calculate vertical convergence (divergence)
        call balance_uvw(domain,options)

    end subroutine update_winds

    !>------------------------------------------------------------
    !! Setup initial fields (i.e. grid relative rotation fields)
    !!
    !!------------------------------------------------------------
    subroutine init_winds(domain,options)
        type(domain_type), intent(inout) :: domain
        type(options_type), intent(in) :: options
        integer:: i,j,nx,ny,starti,endi
        real::dist,dlat,dlon

        call allocate_winds(domain)

        nx=size(domain%lat,1)
        ny=size(domain%lat,2)
        do j=1,ny
            do i=1,nx
                ! in case we are in the first or last grid, reset boundaries
                if (i==1) then
                    starti=i
                    endi=i+1
                elseif (i==nx) then
                    starti=i-1
                    endi=i
                else
                    starti=i-1
                    endi=i+1
                endif

                ! change in latitude
                dlat = domain%lat(endi,j) - domain%lat(starti,j)
                ! change in longitude
                dlon = (domain%lon(endi,j) - domain%lon(starti,j)) * cos(deg2rad*domain%lat(i,j))
                ! distance between two points
                dist = sqrt(dlat**2 + dlon**2)
                ! sin/cos of angles for use in rotating fields later
                domain%sintheta(i,j) = (-1) * dlat / dist
                domain%costheta(i,j) = abs(dlon / dist)
            enddo
        enddo
        if (options%debug) then
            print*, ""
            print*, "Domain Geometry"
            print*, "MAX / MIN SIN(theta) (ideally 0)"
            print*, "   ", maxval(domain%sintheta), minval(domain%sintheta)
            print*, "MAX / MIN COS(theta) (ideally 1)"
            print*, "   ", maxval(domain%costheta), minval(domain%costheta)
            print*, ""
        endif

        ! dzdx/y used to add the effect of terrain following grid on W component of wind field
        domain%dzdx = (domain%terrain(2:nx,:) - domain%terrain(1:nx-1,:)) / domain%dx
        domain%dzdy = (domain%terrain(:,2:ny) - domain%terrain(:,1:ny-1)) / domain%dx

    end subroutine init_winds

    !>------------------------------------------------------------
    !! Allocate memory used in various wind related routines
    !!
    !!------------------------------------------------------------
    subroutine allocate_winds(domain)
        type(domain_type), intent(inout) :: domain
        integer::nx,ny

        nx=size(domain%lat,1)
        ny=size(domain%lat,2)

        if (.not.allocated(domain%sintheta)) then
            allocate(domain%sintheta(nx,ny))
        endif
        if (.not.allocated(domain%costheta)) then
            allocate(domain%costheta(nx,ny))
        endif

        if (.not.allocated(domain%dzdx)) then
            allocate(domain%dzdx(nx-1,ny))
        endif
        if (.not.allocated(domain%dzdy)) then
            allocate(domain%dzdy(nx,ny-1))
        endif

    end subroutine allocate_winds

    !>------------------------------------------------------------
    !! Provides a routine to deallocate memory allocated in allocate_winds
    !!
    !!------------------------------------------------------------
    subroutine finalize_winds(domain)
        type(domain_type), intent(inout) :: domain

        if (allocated(domain%sintheta)) then
            deallocate(domain%sintheta)
        endif
        if (allocated(domain%costheta)) then
            deallocate(domain%costheta)
        endif
        if (allocated(domain%dzdx)) then
            deallocate(domain%dzdx)
        endif
        if (allocated(domain%dzdy)) then
            deallocate(domain%dzdy)
        endif

    end subroutine finalize_winds
end module wind
