!> ----------------------------------------------------------------------------
!!  Driver to call different convection schemes
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module convection
    use data_structures
    use icar_constants
    use module_cu_tiedtke,  only: tiedtkeinit, CU_TIEDTKE
    ! use module_cu_kf,       only: kfinit, KFCPS

    use options_interface,   only : options_t
    use domain_interface,    only : domain_t

    implicit none
    private
    public :: init_convection, convect, cu_var_request

    logical,allocatable, dimension(:,:) :: CU_ACT_FLAG
    real,   allocatable, dimension(:,:) :: XLAND, RAINCV, PRATEC, NCA
    real,   allocatable, dimension(:,:,:):: W0AVG, w_stochastic

    integer :: ids,ide,jds,jde,kds,kde ! Domain dimensions
    integer :: ims,ime,jms,jme,kms,kme ! Local Memory dimensions
    integer :: its,ite,jts,jte,kts,kte ! Processing Tile dimensions

contains

    subroutine cu_var_request(options)
        implicit none
        type(options_t),intent(inout) :: options

        if (options%physics%convection == kCU_TIEDTKE) then
            call options%alloc_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%temperature,                        &
                         kVARS%cloud_water, kVARS%cloud_ice, kVARS%precipitation, kVARS%convective_precipitation,   &
                         kVARS%exner, kVARS%dz_interface, kVARS%density, kVARS%pressure_interface, kVARS%pressure,  &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u, kVARS%v, kVARS%w, kVARS%land_mask,        &
                         kVARS%tend_qv, kVARS%tend_th, kVARS%tend_qc, kVARS%tend_qi, kVARS%tend_qs, kVARS%tend_qr,  &
                         kVARS%tend_u, kVARS%tend_v, kVARS%tend_qv_pbl, kVARS%tend_qv_adv, kVARS%znu])

             call options%advect_vars([kVARS%potential_temperature, kVARS%water_vapor])

             call options%restart_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%temperature,                        &
                         kVARS%cloud_water, kVARS%cloud_ice, kVARS%precipitation, kVARS%convective_precipitation,   &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u, kVARS%v, kVARS%pressure])
        endif

    end subroutine cu_var_request


    subroutine init_convection(domain,options)
        implicit none
        type(domain_t),  intent(inout) :: domain
        type(options_t), intent(in) :: options

        if (this_image()==1) write(*,*) "Initializing Cumulus Scheme"

        ! module level variables for easy access... need to think about tiling to permit halo processing separately.
        ids = domain%grid%ids
        ide = domain%grid%ide
        jds = domain%grid%jds
        jde = domain%grid%jde
        kds = domain%grid%kds
        kde = domain%grid%kde
        ims = domain%grid%ims
        ime = domain%grid%ime
        jms = domain%grid%jms
        jme = domain%grid%jme
        kms = domain%grid%kms
        kme = domain%grid%kme
        its = domain%grid%its
        ite = domain%grid%ite
        jts = domain%grid%jts
        jte = domain%grid%jte
        kts = domain%grid%kts
        kte = domain%grid%kte

        if (options%physics%convection > 0) then
            allocate(CU_ACT_FLAG(ims:ime,jms:jme), source=.False.)
            CU_ACT_FLAG(its:ite, jts:jte) = .True.

            allocate(RAINCV(ims:ime,jms:jme))
            allocate(PRATEC(ims:ime,jms:jme))

            RAINCV = 0
            PRATEC = 0

            allocate(XLAND(ims:ime,jms:jme))
            XLAND = domain%land_mask
            where(domain%land_mask == 0) XLAND = 2 ! 0 is water if using "LANDMASK" as input
        endif

        if (options%physics%convection == kCU_TIEDTKE) then
            if (this_image()==1) write(*,*) "    Tiedtke Cumulus scheme"

            ! allocate(w_stochastic(ims:ime,kms:kme,jms:jme))
            call tiedtkeinit(domain%tend%th,domain%tend%qv,   &
                             domain%tend%qc,domain%tend%qi,   &
                             domain%tend%u, domain%tend%v,    &
                             .false.,1,1,0,                   &
                             .true.,                          &
                             ids,ide, jds,jde, kds,kde,       &
                             ims,ime, jms,jme, kms,kme,       &
                             its,ite, jts,jte, kts,kte  )

         ! elseif (options%physics%convection==kCU_KAINFR) then
         !     write(*,*) "    Kain-Fritsch Cumulus scheme"
         !     allocate(W0AVG(ids:ide,kds:kde,jds:jde))
         !     W0AVG=0
         !     allocate(NCA(ids:ide,jds:jde))
         !     NCA=0
         !     call kfinit(domain%tend%th,domain%tend%qv,   &
         !                 domain%tend%qc,domain%tend%qr,   &
         !                 domain%tend%qi,domain%tend%qs,   &
         !                 NCA,W0AVG,1,1,                   & !...,P_QI,P_QS
         !                 0,.false.,.true.,                & ! P_FIRST_SCALAR, restart, allowed_to_read
         !                 ids, ide, jds, jde, kds, kde,  &
         !                 ids, ide, jds, jde, kds, kde,    &
         !                 ids, ide, jds, jde, kds, kde-1)

         endif

    end subroutine init_convection

subroutine convect(domain,options,dt_in)
    implicit none
    type(domain_t),  intent(inout) :: domain
    type(options_t), intent(in)    :: options
    real, intent(in) :: dt_in

    integer :: j,itimestep,STEPCU
    real :: internal_dt

    if (options%physics%convection==0) return

    itimestep = 1
    STEPCU = 1

    !$omp parallel private(j) &
    !$omp default(shared)
    !$omp do schedule(static)
    do j=jts,jte
        RAINCV(:,j)     = 0
        domain%tend%qv(:,:,j) = 0
        domain%tend%th(:,:,j) = 0
        domain%tend%qc(:,:,j) = 0
        domain%tend%qi(:,:,j) = 0
        domain%tend%qs(:,:,j) = 0
        domain%tend%qr(:,:,j) = 0
        domain%tend%u(:,:,j)  = 0
        domain%tend%v(:,:,j)  = 0
    enddo
    !$omp end do
    !$omp end parallel

    if (options%physics%convection==kCU_TIEDTKE) then
        ! call random_number(w_stochastic)
        ! block
        !     integer :: i
        !     do i=1,size(w_stochastic,2)
        !         w_stochastic(:,i,:) = domain%sensible_heat/500 + w_stochastic(:,i,:)
        !     enddo
        ! end block
        call CU_TIEDTKE(                                        &
                 dt_in, itimestep, STEPCU                       &
                ,RAINCV, PRATEC                                 &
                ,domain%latent_heat%data_2d/LH_vaporization     &
                ,domain%sensible_heat%data_2d                   &
                ,domain%znu                                     &
                ,domain%u_mass%data_3d                          &
                ,domain%v_mass%data_3d                          &
                ,domain%w_real%data_3d                          &
                ! ,domain%w_real%data_3d+(w_stochastic*20-15)   &
                ,domain%temperature%data_3d                     &
                ,domain%water_vapor%data_3d                     &
                ,domain%cloud_water_mass%data_3d                &
                ,domain%cloud_ice_mass%data_3d                  &
                ,domain%exner%data_3d                           &
                ,domain%density%data_3d                         &
                ,domain%tend%qv_adv                             &
                ,domain%tend%qv_pbl                             &
                ,domain%dz_interface%data_3d                    &
                ,domain%pressure%data_3d                        &
                ,domain%pressure_interface%data_3d              &
                ,XLAND, CU_ACT_FLAG                             &
                ,ids,ide, jds,jde, kds,kde                      &
                ,ims,ime, jms,jme, kms,kme                      &
                ,its,ite, jts,jte, kts,kte-1                    &
                ,domain%tend%th, domain%tend%qv, domain%tend%qc &
                ,domain%tend%qi, domain%tend%u,  domain%tend%v  &
                ,.True.,.True.,.True.,.True.,.True.             &
                )

    ! elseif (options%physics%convection==kCU_KAINFR) then
    !     call KFCPS(                                          &
    !            ids,ide, jds,jde, kds,kde                     & ! domain
    !           ,ids,ide, jds,jde, kds,kde                     & ! memory
    !           ,ids+1,ide-1, jds+1,jde-1, kds,kde             & ! tile
    !           ,dt_in,itimestep,domain%dx,dt_in,.false.       & ! dt KTAU, dx, cu_dt, adapt_step
    !           ,domain%rho                                    & ! rho
    !           ,RAINCV,PRATEC,NCA                             & ! convective rain, convective rain rate,
    !           ,domain%Um,domain%Vm,domain%th,domain%t        &
    !           ,domain%w_real                                 &
    !           ,domain%qv,domain%dz_inter,domain%p,domain%pii &
    !           ,W0AVG                                         & ! "average" vertical velocity (computed internally from w)
    !           ,XLV0,XLV1,XLS0,XLS1,cp,Rd,gravity,EP1         & ! physical "constants"
    !           ,EP2,SVP1,SVP2,SVP3,SVPT0                      & ! physical "constants"
    !           ,STEPCU,CU_ACT_FLAG,.false.                    & ! number of steps between CU calls, boolean grid to act on , warm_rain_only
    !         ! optional arguments
    !         ! ,F_QV    ,F_QC    ,F_QR    ,F_QI    ,F_QS      &
    !           ,.True., .True.,  .True.,  .True.,  .True.     &
    !           ,domain%tend%qv, domain%tend%qc                &
    !           ,domain%tend%qr, domain%tend%qi                &
    !           ,domain%tend%qs, domain%tend%th)
    endif

    ! add domain wide tendency terms
    ! use a separate dt to make it easier to apply on a different dt
    internal_dt = dt_in

    ! $omp parallel private(j) &
    ! $omp default(shared)
    ! $omp do schedule(static)
    do j=jts,jte
        domain%water_vapor%data_3d(:,:,j)           = domain%water_vapor%data_3d(:,:,j)           + domain%tend%qv(:,:,j)*internal_dt
        domain%cloud_water_mass%data_3d(:,:,j)      = domain%cloud_water_mass%data_3d(:,:,j)      + domain%tend%qc(:,:,j)*internal_dt
        domain%potential_temperature%data_3d(:,:,j) = domain%potential_temperature%data_3d(:,:,j) + domain%tend%th(:,:,j)*internal_dt
        domain%cloud_ice_mass%data_3d(:,:,j)        = domain%cloud_ice_mass%data_3d(:,:,j)        + domain%tend%qi(:,:,j)*internal_dt
        ! if (options%physics%convection==kCU_KAINFR) then
        !     domain%qsnow(:,:,j) =domain%qsnow(:,:,j) + domain%tend%Qs(:,:,j)*internal_dt
        !     domain%qrain(:,:,j) =domain%qrain(:,:,j) + domain%tend%Qr(:,:,j)*internal_dt
        ! endif

        domain%accumulated_precipitation%data_2d(:,j)  = domain%accumulated_precipitation%data_2d(:,j) + RAINCV(:,j)
        domain%accumulated_convective_pcp%data_2d(:,j) = domain%accumulated_convective_pcp%data_2d(:,j) + RAINCV(:,j)

        ! if (options%physics%convection==kCU_TIEDTKE) then
        !     domain%u_cu(ids+1:ide,:,j) = 0.999*domain%u_cu(ids+1:ide,:,j) + (domain%tend%u(ids:ide-1,:,j)+domain%tend%u(ids+1:ide,:,j))/2 * internal_dt
        !     if (j>jds) then
        !         domain%v_cu(:,:,j) = 0.999*domain%v_cu(:,:,j) + (domain%tend%v(:,:,j)+domain%tend%v(:,:,j-1))/2 * internal_dt
        !     endif
        ! endif
    enddo
    ! $omp end do
    ! $omp end parallel


end subroutine convect
end module convection
