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
    use data_structures
!   use output, only: write_domain
    implicit none
    private
    public::update_winds,balance_uvw,init_winds
    real, parameter::deg2rad=0.017453293 !2*pi/360
contains

! Forces u,v, and w fields to balance
!       du/dx+dv/dy = dw/dz
! Starts by setting w out of the ground=0 then works through layers
    subroutine balance_uvw(domain,options)
        implicit none
        type(domain_type),intent(inout)::domain
        type(options_type),intent(in)::options
        real, allocatable,dimension(:,:) ::du,dv,divergence,rhou,rhov,rhow
        integer ::nx,ny,nz,i
        nx=size(domain%w,1)
        nz=size(domain%w,2)
        ny=size(domain%w,3)
        
        ! these could be module level to prevent lots of allocation/deallocation/reallocations
        ! but these are relatively small, and this only gets called once for each forcing input
        if (options%advect_density) then
            allocate(rhou(nx-1,ny-2))
            allocate(rhov(nx-2,ny-1))
            allocate(rhow(nx-2,ny-2))
        endif
        allocate(du(nx-2,ny-2))
        allocate(dv(nx-2,ny-2))
        allocate(divergence(nx-2,ny-2))
        
! If this becomes a bottle neck in the code it could be parallelized over y  
!       loop over domain levels
        do i=1,nz
            ! if we are incorporating density into the advection equation (as we should?)
            if (options%advect_density) then
                ! first calculate density on the u and v grid points
                rhou=(domain%rho(1:nx-1,i,2:ny-1) + domain%rho(2:nx,i,2:ny-1))/2
                rhov=(domain%rho(2:nx-1,i,1:ny-1) + domain%rho(2:nx-1,i,2:ny))/2
                ! for most grid points interpolate between the current grid and the one above it
                if (i<nz) then
                    rhow=(domain%rho(2:nx-1,i,2:ny-1) + domain%rho(2:nx-1,i+1,2:ny-1))/2
                else
                    ! for the top grid cell extrapolate upwards based on the current grid and the one below it
                    rhow=2*domain%rho(2:nx-1,i,2:ny-1)-domain%rho(2:nx-1,i-1,2:ny-1)
                endif
                ! calculate horizontal divergence based on the wind * density field
                domain%vr(2:nx-1,i,2:ny)=rhov*domain%v(2:nx-1,i,2:ny) * domain%dz_inter(1,i,1)*domain%dx
                dv=domain%vr(2:nx-1,i,3:ny) - domain%vr(2:nx-1,i,2:ny-1)
                
                domain%ur(2:nx,i,2:ny-1)=rhou*domain%u(2:nx,i,2:ny-1) * domain%dz_inter(1,i,1)*domain%dx
                du=domain%ur(3:nx,i,2:ny-1) - domain%ur(2:nx-1,i,2:ny-1)
                
                divergence=du+dv
                if (i==1) then
                    ! if this is the first model level start from 0 at the ground
                    domain%wr(2:nx-1,i,2:ny-1)=0-divergence
                    domain%w(2:nx-1,i,2:ny-1)=0-divergence/rhow/(domain%dx**2)
                else
                    ! else calculate w as a change from w at the level below
                    domain%wr(2:nx-1,i,2:ny-1)=domain%wr(2:nx-1,i-1,2:ny-1)-divergence
                    domain%w(2:nx-1,i,2:ny-1)=domain%w(2:nx-1,i-1,2:ny-1)-divergence/rhow/(domain%dx**2)
                endif
            else
                ! calculate horizontal divergence
                dv=domain%v(2:nx-1,i,3:ny) - domain%v(2:nx-1,i,2:ny-1)
                du=domain%u(3:nx,i,2:ny-1) - domain%u(2:nx-1,i,2:ny-1)
                divergence=du+dv
                if (i==1) then
                    ! if this is the first model level start from 0 at the ground
                    domain%w(2:nx-1,i,2:ny-1)=0-divergence
                else
                    ! else calculate w as a change from w at the level below
                    domain%w(2:nx-1,i,2:ny-1)=domain%w(2:nx-1,i-1,2:ny-1)-divergence
                endif
            endif
        enddo
        ! NOTE w is scaled by dx/dz because it is balancing divergence through a dx*dz cell face with a flow through a dx*dx face
        ! could rescale below but for now this is left in the advection code (the only place it matters for now)
        ! this makes it easier to play with varying options for including varying dz and rho in advection. 
        ! domain%w=domain%w/domain%dx * domain%dz_inter
        
        deallocate(du,dv,divergence)
        if (options%advect_density) then
            deallocate(rhou,rhov,rhow)
        endif
    end subroutine balance_uvw
    
!   Correct for a grid that is locally rotated with respect to EW,NS (e.g. at the edges of the domain)
!   Assumes forcing winds are EW,NS relative, not grid relative. 
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


    ! apply wind field physics and adjustments
    ! this will call the linear wind module if necessary, otherwise it just updates for 
    ! this should ONLY be called once for each forcing step, otherwise effects will be additive. 
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
        
        ! linear winds
        if (options%physics%windtype==kWIND_LINEAR) then
            call linear_perturb(domain,options,options%lt_options%vert_smooth,.False.,options%advect_density)
        endif
        ! else assumes even flow over the mountains
        
        ! use horizontal divergence (convergence) to calculate vertical convergence (divergence)
        call balance_uvw(domain,options)
        
    end subroutine update_winds
    
    !setup initial fields (i.e. grid relative rotation fields)
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
                dlat=domain%lat(endi,j)-domain%lat(starti,j)
                ! change in longitude
                dlon=(domain%lon(endi,j)-domain%lon(starti,j))*cos(deg2rad*domain%lat(i,j))
                ! distance between two points
                dist=sqrt(dlat**2+dlon**2)
                ! sin/cos of angles for use in rotating fields later
                domain%sintheta(i,j)=(-1)*dlat/dist
                domain%costheta(i,j)=abs(dlon/dist)
            enddo
        enddo
!         if (options%debug) then
!             print*, "Domain Geometry"
!             print*, "MAX / MIN SIN(theta) (ideally 0)"
!             print*, "   ", maxval(domain%sintheta), minval(domain%sintheta)
!             print*, "MAX / MIN COS(theta) (ideally 1)"
!             print*, "   ", maxval(domain%costheta), minval(domain%costheta)
!             print*, " "
!         endif
!       dzdx/y used effect of terrain following grid on W component of wind field
        domain%dzdx = (domain%terrain(2:nx,:)-domain%terrain(1:nx-1,:)) / domain%dx
        domain%dzdy = (domain%terrain(:,2:ny)-domain%terrain(:,1:ny-1)) / domain%dx
        
    end subroutine init_winds
    
! allocate memory
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

! deallocate memory
    subroutine finalize_winds(domain)
        type(domain_type), intent(inout) :: domain
        
        if (allocated(domain%sintheta)) then
            deallocate(domain%sintheta)
        endif
        if (allocated(domain%costheta)) then
            deallocate(domain%costheta)
        endif
    
    end subroutine finalize_winds
end module wind
