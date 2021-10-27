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

    integer :: ims, ime, jms, jme, kms, kme
    integer :: ids, ide, jds, jde, kds, kde


!   limits on Pr noted in HP96 page 2325 below eqn 13
    real, parameter :: pr_upper_limit = 4.0 !Prandtl number for stability
    real, parameter :: pr_lower_limit = 0.25 !Prandtl number for stability
    real, parameter :: asymp_length_scale = 1/250.0 !m from Hong and Pan (1996)
    ! note, they actually use 30m because they only use this for free-atmosphere mixing
    ! but they note that 250m is used in the operational model for the full PBL mixing
    real, parameter :: N_substeps=10. ! number of substeps to allow (puts a cap on K to match CFL)
    real, parameter :: diffusion_reduction=10.0 ! used to reduce diffusion rates


contains
    subroutine simple_pbl(th, qv, cloud, ice, qrain, qsnow, um, vm, pii, rho, z, dz, terrain, its, ite, jts, jte, kts, kte_in, dt, tend_qv_pbl)
        real,   intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: th            ! potential temperature [K]
        real,   intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: qv            ! water vapor mixing ratio [kg/kg]
        real,   intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: cloud         ! cloud water mixing ratio [kg/kg]
        real,   intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: ice           ! cloud ice mixing ratio [kg/kg]
        real,   intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: qrain         ! rain water mixing ratio [kg/kg]
        real,   intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: qsnow         ! snow mixing ratio [kg/kg]
        real,   intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: um            ! east-west wind speed on mass grid [m/s]
        real,   intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: vm            ! north south wind speed on mass grid [m/s]
        real,   intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: pii           ! exner function
        real,   intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: rho           ! air density [kg / m^3]
        real,   intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: z             ! model level heights [m]
        real,   intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: dz            ! model level thickness [m]
        real,   intent(in),    dimension(ims:ime, jms:jme)          :: terrain       ! terrain height above sea level [m]
        integer,intent(in) :: its, ite, jts, jte, kts, kte_in
        real,   intent(in) :: dt                                 ! time step [s]
        real,   intent(inout), dimension(ims:ime,kms:kme,jms:jme), optional :: tend_qv_pbl   ! output water vapor tendency [kg/kg/s] (for use in other physics)

        ! local
        integer :: i,j,k, kte

        ! don't process the top model level regardless, there shouldn't be any real pbl diffusion occuring there
        kte = min(kme-1, kte_in)
!       OpenMP parallelization small static chunk size because we typically get a small area that takes most of the time (because of substepping)
        !$omp parallel shared(th, qv, cloud, ice, qrain, qsnow, um, vm, pii, rho, z, dz, terrain) & !, tend_qv_pbl)     &
        !$omp shared(l_m, K_m, Kq_m, stability_m, prandtl_m, virt_pot_temp_zgradient_m, rig_m, shear_m, lastqv_m, ims, ime, jms, jme, kms, kme) &
        !$omp firstprivate(its, ite, jts, jte, kts, kte, dt) private(i, k, j)

        !$omp do schedule(static, 2)
        do j = jts, jte
            do k = kts, kte
                lastqv_m(:,k,j) = qv(:,k,j)
                call calc_shear(um, vm, dz, its, ite, k, j)
                call calc_virt_pot_temp_zgradient(th, qv, dz, cloud, ice, qrain, qsnow, its, ite, k, j)
                call calc_richardson_gradient(th, pii, its, ite, k, j)
                call calc_pbl_stability_function(its, ite, k, j)

                ! from eqn 12 in HP96
                l_m(its:ite,k,j) = 1 / (1/(karman*(z(its:ite,k,j) - terrain(its:ite,j))) + asymp_length_scale)

                ! diffusion for momentum... can I ignore this term and go directly to scalars (to save memory)?
                ! from HP96 eqn 11
                ! k = l**2 * stability * shear * dt/dz
                K_m(its:ite,k,j) = l_m(its:ite,k,j)**2 * stability_m(its:ite,k,j) * shear_m(its:ite,k,j)

                ! diffusion for scalars
                Kq_m(its:ite,k,j) = K_m(its:ite,k,j) / prandtl_m(its:ite,k,j)

                ! rescale diffusion to cut down on excessive mixing
                Kq_m(its:ite,k,j) = Kq_m(its:ite, k,j) / diffusion_reduction
                Kq_m(its:ite,k,j) = Kq_m(its:ite, k,j) * dt / ((dz(its:ite,k,j) + dz(its:ite,k+1,j))/2)

                ! enforce limits specified in HP96
                do i=its,ite
                    if (Kq_m(i,k,j)>1000) then
                        Kq_m(i,k,j)=1000
                    elseif (Kq_m(i,k,j)<1) then
                        Kq_m(i,k,j)=1
                    endif
                enddo
            enddo

            call pbl_diffusion(qv, th, cloud, ice, qrain, qsnow, rho, dz, its, ite, kts, kte, j)

            ! tend_qv_pbl(:,:,j+1) = (qv(:,:,j+1) - lastqv_m(:,:,j+1)) / dt
        enddo
        !$omp end do
        !$omp end parallel

    end subroutine simple_pbl

    subroutine diffuse_variable(q, rhomean, rho_dz, its, ite, kts, kte, j)
        real,   intent(inout),  dimension(ims:ime,kms:kme,jms:jme) :: q
        real,   intent(in),     dimension(ims:ime,kms:kme) :: rhomean,rho_dz
        integer,intent(in) :: its, ite, kts, kte, j

        real, dimension(its:ite,kts:kte) :: fluxes
        integer :: k

        do k = kts, kte
            ! Eventually this should be made into an implicit solution to avoid substepping
            ! if gradient is downward (qv[z+1]>qv[z]) then flux is negative
            fluxes(:, k) = Kq_m(its:ite, k, j) * rhomean(its:ite, k) * (q(its:ite, k, j) - q(its:ite, k+1, j))
        enddo
        ! first layer assumes no flow through the surface, that comes from the LSM
        q(its:ite, kts,j) = q(its:ite, kts, j) - fluxes(its:ite,kts) / rho_dz(its:ite, kts)
        ! middle layers (no change for top layer assuming flux in = flux out)
        q(its:ite, kts+1:kte,j) = q(its:ite, kts+1:kte, j) - (fluxes(its:ite,kts+1:kte) - fluxes(its:ite,kts:kte-1)) / rho_dz(its:ite,kts+1:kte)
        ! make sure top fluxes end up in top layer
        q(its:ite, kte+1,j) = q(its:ite, kte+1, j) + fluxes(its:ite,kte) / rho_dz(its:ite, kte)

    end subroutine diffuse_variable

    subroutine pbl_diffusion(qv, th, cloud, ice, qrain, qsnow, rho, dz, its, ite, kts, kte, j)
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: th            ! potential temperature [K]
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: qv            ! water vapor mixing ratio [kg/kg]
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: cloud         ! cloud water mixing ratio [kg/kg]
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: ice           ! cloud ice mixing ratio [kg/kg]
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: qrain         ! rain water mixing ratio [kg/kg]
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: qsnow         ! snow mixing ratio [kg/kg]
        real, intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: rho           ! air density [kg / m^3]
        real, intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: dz            ! model level thickness [m]
        integer,intent(in) :: its, ite, kts, kte, j

        ! locals
        integer :: i, k, nsubsteps, t
        real, dimension(ims:ime, kms:kme) :: fluxes, rhomean, rho_dz

        do k = kts, kte
            rho_dz(its:ite, k)  = dz(its:ite, k,j) * rho(its:ite, k, j)
            rhomean(its:ite, k) = (rho(its:ite, k, j) + rho(its:ite, k+1, j)) / 2
        enddo
        ! note Kq_m already has dt/dz embedded in it
        ! diffusion fluxes within the PBL
        ! q = q + (k dq/dz)/dz *dt

        !if K >1 we are in violation of the CFL condition and we need to subset (or make implicit...)
        ! for most regions it is < 0.5, for the small regions it isn't just substep for now.
        ! nsubsteps will often be 1, but allow up to N sub-steps in extreme cases
        where((Kq_m(its:ite, kts:kte,j)) > N_substeps*dz(its:ite, kts:kte,j))   &
              Kq_m(its:ite, kts:kte, j) = dz(its:ite, kts:kte, j) * N_substeps

        nsubsteps = ceiling( 2 * maxval(Kq_m(its:ite, kts:kte, j) / dz(its:ite, kts:kte, j)))
        Kq_m(its:ite, kts:kte, j) = Kq_m(its:ite, kts:kte, j) / nsubsteps
        do t = 1, nsubsteps
            ! First water vapor
            call diffuse_variable(qv, rhomean, rho_dz, its, ite, kts, kte, j)
            ! ditto for potential temperature
            call diffuse_variable(th, rhomean, rho_dz, its, ite, kts, kte, j)
            ! and cloud water
            call diffuse_variable(cloud, rhomean, rho_dz, its, ite, kts, kte, j)
            ! and cloud ice
            call diffuse_variable(ice, rhomean, rho_dz, its, ite, kts, kte, j)
            ! and snow
            call diffuse_variable(qsnow, rhomean, rho_dz, its, ite, kts, kte, j)
            ! and rain
            call diffuse_variable(qrain, rhomean, rho_dz, its, ite, kts, kte, j)
            ! don't bother with graupel assuming they are falling fast *enough* not entirely fair...
        enddo
    end subroutine pbl_diffusion

    subroutine calc_shear(u_mass, v_mass, dz, its, ite, k, j)
        real, intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: u_mass     ! east-west wind speed [m/s]
        real, intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: v_mass     ! north south wind speed [m/s]
        real, intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: dz         ! model level thickness [m]
        integer,intent(in) :: its, ite, k, j

        shear_m(its:ite, k, j) = sqrt(  (u_mass(its:ite, k+1, j) - u_mass(its:ite, k, j))**2   &
                                      + (v_mass(its:ite, k+1, j) - v_mass(its:ite, k, j))**2)  &
                                  / ((dz(its:ite, k, j) + dz(its:ite, k+1,j)) * 0.5)

        where(shear_m(its:ite, k,j)<1e-5) shear_m(its:ite, k,j) = 1e-5
    end subroutine calc_shear

!   calculate the vertical gradient in virtual potential temperature
    subroutine calc_virt_pot_temp_zgradient(th, qv, dz, cloud, ice, qrain, qsnow, its, ite, k, j)
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: th            ! potential temperature [K]
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: qv            ! water vapor mixing ratio [kg/kg]
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: cloud         ! cloud water mixing ratio [kg/kg]
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: ice           ! cloud ice mixing ratio [kg/kg]
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: qrain         ! rain water mixing ratio [kg/kg]
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: qsnow         ! snow mixing ratio [kg/kg]
        real, intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: dz            ! model level thickness [m]

        real, dimension(ims:ime, 2) :: vth         ! virtual potential temperature
        integer,intent(in) :: its, ite, k, j

        ! first calculate the virtual potential temperature
        ! vth=th*(1+0.61*qv-(qc+qi+qr+qs))
        vth(its:ite, :) = th(its:ite, k:k+1, j)                 &
                          * (1 + 0.61 * qv(its:ite, k:k+1, j)   &
                          - (cloud(its:ite, k:k+1, j) + ice(its:ite, k:k+1, j) + qrain(its:ite, k:k+1, j) + qsnow(its:ite, k:k+1, j)))

        virt_pot_temp_zgradient_m(its:ite, k, j) = (vth(its:ite, 2) - vth(its:ite, 1)) &
                                             / ((dz(its:ite, k, j) + dz(its:ite, k+1, j)) * 0.5)

    end subroutine calc_virt_pot_temp_zgradient

    ! calculate the stability function based on HP96
    subroutine calc_pbl_stability_function(its, ite, k, j)
        integer,intent(in) :: its, ite, k, j
        integer :: i
        ! HP96 eqn 13
        where (rig_m(its:ite, k, j) > 0) stability_m(its:ite, k, j) = exp(-8.5 * rig_m(its:ite, k, j)) + 0.15 / (rig_m(its:ite, k, j)+3)
        where (rig_m(its:ite, k, j) <= 0) stability_m(its:ite, k, j) = 1 / sqrt(1-1.6*rig_m(its:ite, k, j))

        ! HP96 eqn 13 continued
        prandtl_m(its:ite, k, j) = 1.5 + 3.08 * rig_m(its:ite, k, j)

        ! Impose limits as specified
        ! on Pr noted in HP96
        do i = its, ite
            if (prandtl_m(i,k,j) > pr_upper_limit) then
                prandtl_m(i,k,j) = pr_upper_limit
            elseif (prandtl_m(i,k,j) < pr_lower_limit) then
                prandtl_m(i,k,j) = pr_lower_limit
            endif
        enddo
!       prandtl_m(i,j,k) = min( pr_upper_limit, max( pr_lower_limit, prandtl_m(i,j,k)) )
!       alternatively... which is faster? the following *might* vectorize better, but memory bandwidth may be the limit anyway
!       where(prandtl_m(:,:,j)>4) prandtl_m(:,:,j)=4
!       where(prandtl_m(:,:,j)<0.25) prandtl_m(:,:,j)=0.25
    end subroutine calc_pbl_stability_function

!   calculate the gradient in the richardson number as specified in HP96
!   rig = Richardson number Gradient
    subroutine calc_richardson_gradient(th, pii, its, ite, k, j)
        ! calculate the local gradient richardson number as in eqn. between 11 and 12 in HP96
        real,   intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: th            ! potential temperature [K]
        real,   intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: pii           ! exner function
        integer,intent(in) :: its, ite, k, j

        real, dimension(its:ite) :: temperature
        ! might be slightly better to interpolate theta to half levels, then recalc p and pii at half levels
        temperature = (th(its:ite, k, j) * pii(its:ite, k, j) + th(its:ite, k+1, j) * pii(its:ite, k+1, j)) / 2
        rig_m(its:ite, k, j) =  gravity/temperature  &
                       * virt_pot_temp_zgradient_m(its:ite, k, j) * 1 / (shear_m(its:ite, k, j)**2)
        where(rig_m(its:ite, k, j)<-100.0) rig_m(its:ite, k, j)=-100.0

    end subroutine calc_richardson_gradient


! memory allocation and deallocation
! Can/Should also add parameter definition from options%pbl (or something like that)
    subroutine init_simple_pbl(domain,options)
        type(domain_t), intent(in) :: domain
        type(options_t),intent(in) :: options

        ! module level variables assumed constant throughout the run...
        ims = domain%ims
        ime = domain%ime
        jms = domain%jms
        jme = domain%jme
        kms = domain%kms
        kme = domain%kme
        ids = domain%ids
        ide = domain%ide
        jds = domain%jds
        jde = domain%jde
        kds = domain%kds
        kde = domain%kde

        allocate(virt_pot_temp_zgradient_m(ims:ime,kms:kme,jms:jme))
        allocate(rig_m                    (ims:ime,kms:kme,jms:jme))
        allocate(stability_m              (ims:ime,kms:kme,jms:jme))
        allocate(shear_m                  (ims:ime,kms:kme,jms:jme))
        allocate(prandtl_m                (ims:ime,kms:kme,jms:jme))
        allocate(K_m                      (ims:ime,kms:kme,jms:jme))
        allocate(Kq_m                     (ims:ime,kms:kme,jms:jme))
        allocate(l_m                      (ims:ime,kms:kme,jms:jme))
        allocate(lastqv_m                 (ims:ime,kms:kme,jms:jme))
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
