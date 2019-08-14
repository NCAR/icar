!>----------------------------------------------------------
!!  Simple PBL diffusion package for ICAR
!!
!!  Local-K diffusion type PBL as in Louis (1979) as documented in Hong and Pan (1996) = HP96
!!  Hong and Pan used this for their free atmosphere diffusion, but noted differences
!!  used in the "current operational model" notably the asymptotic length scale lambda
!!
!! <pre>
!! HP96 = Hong,S.-Y. and H.-L. Pan (1996) Monthly Weather Review v127 p2322
!!       Nonlocal Boundary Layer Vertical Diffusion in a Medium Range Forecast Model
!!
!! Implemented with K,shear,stability... on half levels
!!  rho on half levels for f=k*rho*dq/dz*dt
!!  rho on full levels for q=q+f/rho
!!   q,U,V on full levels
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module pbl_simple
    use data_structures
    use domain_interface,   only : domain_t
    use options_interface,  only : options_t

    private
    public :: simple_pbl, finalize_simple_pbl, init_simple_pbl

!   NOTE *_m indicates module level variables
!   these variables are declared as module level variables so that they do not need to be allocated
!   deallocated, and re-allocated all the time.

!   gradient in the virtual potential temperature
    real, allocatable, dimension(:,:,:) :: virt_pot_temp_zgradient_m
!   gradient in the richardson number
    real, allocatable, dimension(:,:,:) :: rig_m
!   vertical wind shear (dwind / dz)
    real, allocatable, dimension(:,:,:) :: shear_m
!   atmospheric stability
    real, allocatable, dimension(:,:,:) :: stability_m
!   length scale that asymptotes from karman*z to l0 (250m below)
    real, allocatable, dimension(:,:,:) :: l_m
!   diffusion term for momentum
    real, allocatable, dimension(:,:,:) :: K_m
!   diffusion term for scalars (K/prandtl)
    real, allocatable, dimension(:,:,:) :: Kq_m
!   prandtl number to convert K for momentum to K for scalars
    real, allocatable, dimension(:,:,:) :: prandtl_m
!   input qv field to use in calculateing the qv_pbl_tendency
    real, allocatable, dimension(:,:,:) :: lastqv_m
    integer :: nx,nz,ny !NOTE these are subset from full domain e.g. nx-2,ny-2,nz-1

!   limits on Pr noted in HP96 page 2325 below eqn 13
    real, parameter :: pr_upper_limit = 4.0 !Prandtl number for stability
    real, parameter :: pr_lower_limit = 0.25 !Prandtl number for stability
    real, parameter :: asymp_length_scale = 1/250.0 !m from Hong and Pan (1996)
    ! note, they actually use 30m because they only use this for free-atmosphere mixing
    ! but they note that 250m is used in the operational model for the full PBL mixing
    real, parameter :: N_substeps=10. ! number of substeps to allow (puts a cap on K to match CFL)
    real, parameter :: diffusion_reduction=10.0 ! used to reduce diffusion rates


contains
    subroutine simple_pbl(th, qv, cloud, ice, qrain, qsnow, u, v, pii, rho, z, dz, tend_qv_pbl, terrain, dt)
        real, intent(inout), dimension(:,:,:) :: th            ! potential temperature [K]
        real, intent(inout), dimension(:,:,:) :: qv            ! water vapor mixing ratio [kg/kg]
        real, intent(inout), dimension(:,:,:) :: cloud         ! cloud water mixing ratio [kg/kg]
        real, intent(inout), dimension(:,:,:) :: ice           ! cloud ice mixing ratio [kg/kg]
        real, intent(inout), dimension(:,:,:) :: qrain         ! rain water mixing ratio [kg/kg]
        real, intent(inout), dimension(:,:,:) :: qsnow         ! snow mixing ratio [kg/kg]
        real, intent(in),    dimension(:,:,:) :: u             ! east-west wind speed [m/s]
        real, intent(in),    dimension(:,:,:) :: v             ! north south wind speed [m/s]
        real, intent(in),    dimension(:,:,:) :: pii           ! exner function
        real, intent(in),    dimension(:,:,:) :: rho           ! air density [kg / m^3]
        real, intent(in),    dimension(:,:,:) :: z             ! model level heights [m]
        real, intent(in),    dimension(:,:,:) :: dz            ! model level thickness [m]
        real, intent(in),    dimension(:,:)   :: terrain       ! terrain height above sea level [m]
        real, intent(inout),   dimension(:,:,:) :: tend_qv_pbl   ! output water vapor tendency [kg/kg/s] (for use in other physics)
        real, intent(in) :: dt                                 ! time step [s]

        ! local
        integer :: i,j,k

!       OpenMP parallelization small static chunk size because we typically get a small area that takes most of the time (because of substepping)
        !$omp parallel shared(th, qv, cloud, ice, qrain, qsnow, u, v, pii, rho, z, dz, terrain, tend_qv_pbl)     &
        !$omp shared(l_m, K_m, Kq_m, stability_m, prandtl_m, virt_pot_temp_zgradient_m, rig_m, shear_m, lastqv_m) &
        !$omp firstprivate(nx, nz, ny, dt) private(i, k, j)

        !$omp do schedule(static, 2)
        do j = 1, ny
            lastqv_m(:,:,j+1) = qv(:,:,j+1)
            call calc_shear(u, v, dz, j)
            call calc_virt_pot_temp_zgradient(th, qv, dz, cloud, ice, qrain, qsnow,j)
            call calc_richardson_gradient(th, pii, j)
            call calc_pbl_stability_function(j)

!           from eqn 12 in HP96
            do k=1,nz
                l_m(:,k,j) = 1 / (1/(karman*(z(2:nx+1,k,j+1)-terrain(2:nx+1,j+1))) + asymp_length_scale)
            enddo
!           diffusion for momentum... can I ignore this term and go directly to scalars (to save memory)?
!           from HP96 eqn 11
!           k = l**2 * stability * shear * dt/dz
            K_m(:,:,j) = l_m(:,:,j)**2 * stability_m(:,:,j) * shear_m(:,:,j)
!           diffusion for scalars
            Kq_m(:,:,j) = K_m(:,:,j) / prandtl_m(:,:,j)
!           enforce limits specified in HP96
            do k=1,nz
                do i=1,nx
                    if (Kq_m(i,k,j)>1000) then
                        Kq_m(i,k,j)=1000
                    elseif (Kq_m(i,k,j)<1) then
                        Kq_m(i,k,j)=1
                    endif
                enddo
            enddo
            ! arbitrarily rescale diffusion to cut down on what seems to be excessive mixing
            Kq_m(:,:,j) = Kq_m(:,:,j) / diffusion_reduction
            Kq_m(:,:,j)= Kq_m(:,:,j)* dt / ((dz(2:nx+1,2:,j+1)+dz(2:nx+1,:nz,j+1))/2)

            call pbl_diffusion(qv, th, cloud, ice, qrain, qsnow, rho, dz ,j)

            ! tend_qv_pbl(:,:,j+1) = (qv(:,:,j+1) - lastqv_m(:,:,j+1)) / dt
        enddo
        !$omp end do
        !$omp end parallel

    end subroutine simple_pbl

    subroutine diffuse_variable(q,rhomean,rho_dz,j)
        real, intent(inout),dimension(nx+2,nz+1,ny+2) :: q
        real, intent(in),dimension(nx,nz) :: rhomean,rho_dz
        integer,intent(in) :: j
        real,dimension(nx,nz)::fluxes

!       Eventually this should be made into an implicit solution to avoid substepping
        ! if gradient is downward (qv[z+1]>qv[z]) then flux is negative
        fluxes=Kq_m(:,:,j)*rhomean*(q(2:nx+1,:nz,j+1)-q(2:nx+1,2:,j+1))
        ! first layer assumes no flow through the surface, that comes from the LSM
        q(2:nx+1,1,j+1)    = q(2:nx+1,1,j+1)    - fluxes(:,1) / rho_dz(:,1)
        ! middle layers (no change for top layer assuming flux in = flux out)
        q(2:nx+1,2:nz,j+1) = q(2:nx+1,2:nz,j+1) - (fluxes(:,2:nz)-fluxes(:,:nz-1)) / rho_dz(:,2:nz)

    end subroutine diffuse_variable

    subroutine pbl_diffusion(qv, th, cloud, ice, qrain, qsnow, rho, dz ,j)
        real, intent(inout), dimension(:,:,:) :: th            ! potential temperature [K]
        real, intent(inout), dimension(:,:,:) :: qv            ! water vapor mixing ratio [kg/kg]
        real, intent(inout), dimension(:,:,:) :: cloud         ! cloud water mixing ratio [kg/kg]
        real, intent(inout), dimension(:,:,:) :: ice           ! cloud ice mixing ratio [kg/kg]
        real, intent(inout), dimension(:,:,:) :: qrain         ! rain water mixing ratio [kg/kg]
        real, intent(inout), dimension(:,:,:) :: qsnow         ! snow mixing ratio [kg/kg]
        real, intent(in),    dimension(:,:,:) :: rho           ! air density [kg / m^3]
        real, intent(in),    dimension(:,:,:) :: dz            ! model level thickness [m]
        integer,intent(in) :: j

        ! locals
        integer :: i, nsubsteps, t
        real, dimension(nx,nz) :: fluxes, rhomean, rho_dz

        rhomean = (rho(2:nx+1,1:nz,j+1) + rho(2:nx+1,2:nz+1,j+1)) / 2
        rho_dz  = dz(2:nx+1,2:nz+1,j+1) * rho(2:nx+1,2:nz+1,j+1)

        ! note Kq_m already has dt/dz embedded in it
        ! diffusion fluxes within the PBL
        ! q = q + (k dq/dz)/dz *dt

        !if K >1 we are in violation of the CFL condition and we need to subset (or make implicit...)
        ! for most regions it is < 0.5, for the small regions it isn't just substep for now.
        ! nsubsteps will often be 1, but allow up to N sub-steps in extreme cases
        where((Kq_m(:,:,j))>N_substeps*dz(2:nx+1,:nz,j)) Kq_m(:,:,j)=dz(2:nx+1,:nz,j)*N_substeps
        nsubsteps = ceiling( 2 * maxval(Kq_m(:,:,j) / dz(2:nx+1,:nz,j)))
        Kq_m(:,:,j) = Kq_m(:,:,j) / nsubsteps
        do t = 1, nsubsteps
            ! First water vapor
            call diffuse_variable(qv, rhomean, rho_dz, j)
            ! ditto for potential temperature
            call diffuse_variable(th, rhomean, rho_dz, j)
            ! and cloud water
            call diffuse_variable(cloud, rhomean, rho_dz, j)
            ! and cloud ice
            call diffuse_variable(ice, rhomean, rho_dz, j)
            ! and snow
            call diffuse_variable(qsnow, rhomean, rho_dz, j)
            ! and rain
            call diffuse_variable(qrain, rhomean, rho_dz, j)
            ! don't bother with graupel assuming they are falling fast *enough* not entirely fair...
        enddo
    end subroutine pbl_diffusion

    subroutine calc_shear(u, v, dz, j)
        real, intent(in),    dimension(:,:,:) :: u             ! east-west wind speed [m/s]
        real, intent(in),    dimension(:,:,:) :: v             ! north south wind speed [m/s]
        real, intent(in),    dimension(:,:,:) :: dz            ! model level thickness [m]
        integer,intent(in) :: j
        real,dimension(nx,nz+1)::centered_winds
        integer::last_wind,k

        centered_winds(:,:) = sqrt( ((u(2:nx+1,:,j+1) + u(3:nx+2,:,j+1)) / 2)**2 &
                                   +((v(2:nx+1,:,j)   + v(2:nx+1,:,j+1)) / 2)**2)

        shear_m(:,:,j) = abs(centered_winds(:,2:nz+1) - centered_winds(:,1:nz))  &
                         / ((dz(2:nx+1,1:nz,j+1) + dz(2:nx+1,2:nz+1,j+1)) * 0.5)

        where(shear_m(:,:,j)<1e-5) shear_m(:,:,j)=1e-5
    end subroutine calc_shear

!   calculate the vertical gradient in virtual potential temperature
    subroutine calc_virt_pot_temp_zgradient(th, qv, dz, cloud, ice, qrain, qsnow, j)
        real, intent(inout), dimension(:,:,:) :: th            ! potential temperature [K]
        real, intent(inout), dimension(:,:,:) :: qv            ! water vapor mixing ratio [kg/kg]
        real, intent(inout), dimension(:,:,:) :: cloud         ! cloud water mixing ratio [kg/kg]
        real, intent(inout), dimension(:,:,:) :: ice           ! cloud ice mixing ratio [kg/kg]
        real, intent(inout), dimension(:,:,:) :: qrain         ! rain water mixing ratio [kg/kg]
        real, intent(inout), dimension(:,:,:) :: qsnow         ! snow mixing ratio [kg/kg]
        real, intent(in),    dimension(:,:,:) :: dz            ! model level thickness [m]
        integer,intent(in)::j
        integer::k

!       first calculate the virtual potential temperature
!       vth=th*(1+0.61*qv-(qc+qi+qr+qs))
!       note that domain variables are accessed at 2:nx-1 and j+1 because vth is only over the processing domain (1:n-2)
        virt_pot_temp_zgradient_m(:,:,j) = th(2:nx+1,:,j+1)* &
                (1+0.61*qv(2:nx+1,:,j+1) &
                 -(cloud(2:nx+1,:,j+1)+ice(2:nx+1,:,j+1)+qrain(2:nx+1,:,j+1)+qsnow(2:nx+1,:,j+1)))
        do k=1,nz
            virt_pot_temp_zgradient_m(:,k,j)=(virt_pot_temp_zgradient_m(:,k+1,j)-virt_pot_temp_zgradient_m(:,k,j)) &
                                             / ((dz(2:nx+1,k,j+1)+dz(2:nx+1,k+1,j+1))*0.5)
        enddo

    end subroutine calc_virt_pot_temp_zgradient

    ! calculate the stability function based on HP96
    subroutine calc_pbl_stability_function(j)
        integer,intent(in) :: j
        integer :: i, k
!       HP96 eqn 13
        stability_m(:,:,j) = exp(-8.5*rig_m(:,:,j)) + 0.15 / (rig_m(:,:,j)+3)
!       HP96 eqn 13 continued
        prandtl_m(:,:,j) = 1.5+3.08*rig_m(:,:,j)

!       Impose limits as specified
!       on Pr noted in HP96
        do k=1,size(prandtl_m,2)
            do i=1,size(prandtl_m,1)
                if (prandtl_m(i,k,j)>pr_upper_limit) then
                    prandtl_m(i,k,j)=pr_upper_limit
                elseif (prandtl_m(i,k,j)<pr_lower_limit) then
                    prandtl_m(i,k,j)=pr_lower_limit
                endif
            enddo
        enddo
!       alternatively... which is faster? the following *might* vectorize better, but memory bandwidth may be the limit anyway
!       where(prandtl_m(:,:,j)>4) prandtl_m(:,:,j)=4
!       where(prandtl_m(:,:,j)<0.25) prandtl_m(:,:,j)=0.25
    end subroutine calc_pbl_stability_function

!   calculate the gradient in the richardson number as specified in HP96
!   rig = Richardson number Gradient
    subroutine calc_richardson_gradient(th, pii, j)
        ! calculate the local gradient richardson number as in eqn. between 11 and 12 in HP96
        real, intent(inout), dimension(:,:,:) :: th            ! potential temperature [K]
        real, intent(in),    dimension(:,:,:) :: pii           ! exner function

        integer,intent(in)::j
        real, dimension(nx,nz+1) :: temperature
        temperature = th(2:nx+1,:,j+1) * pii(2:nx+1,:,j+1)
!       might be slightly better to interpolate theta to half levels, then recalc p and pii at half levels
        temperature(:,:nz)= (temperature(:,:nz)+temperature(:,2:))*0.5
        rig_m(:,:,j) =  gravity/temperature(:,:nz)  &
                       * virt_pot_temp_zgradient_m(:,:nz,j) * 1/(shear_m(:,:,j)**2)
        where(rig_m(:,:,j)<-2.9) rig_m(:,:,j)=-2.9

    end subroutine calc_richardson_gradient


! memory allocation and deallocation
! Can/Should also add parameter definition from options%pbl (or something like that)
    subroutine init_simple_pbl(domain,options)
        type(domain_t), intent(in) :: domain
        type(options_t),intent(in) :: options
        nx = domain%nx - 2
        nz = domain%nz - 1
        ny = domain%ny - 2

        allocate(virt_pot_temp_zgradient_m(nx,nz+1,ny))
        allocate(rig_m(nx,nz,ny))
        allocate(stability_m(nx,nz,ny))
        allocate(shear_m(nx,nz,ny))
        allocate(prandtl_m(nx,nz,ny))
        allocate(K_m(nx,nz,ny))
        allocate(Kq_m(nx,nz,ny))
        allocate(l_m(nx,nz,ny))
        allocate(lastqv_m(nx+2,nz+1,ny+2))
    end subroutine init_simple_pbl

!   deallocate memory if requested
    subroutine finalize_simple_pbl()
        if (allocated(virt_pot_temp_zgradient_m)) then
            deallocate(virt_pot_temp_zgradient_m)
        endif
        if (allocated(rig_m)) then
            deallocate(rig_m)
        endif
        if (allocated(stability_m)) then
            deallocate(stability_m)
        endif
        if (allocated(shear_m)) then
            deallocate(shear_m)
        endif
        if (allocated(prandtl_m)) then
            deallocate(prandtl_m)
        endif
        if (allocated(K_m)) then
            deallocate(K_m)
        endif
        if (allocated(Kq_m)) then
            deallocate(Kq_m)
        endif
        if (allocated(l_m)) then
            deallocate(l_m)
        endif
        if (allocated(lastqv_m)) then
            deallocate(lastqv_m)
        endif
    end subroutine finalize_simple_pbl
end module pbl_simple
