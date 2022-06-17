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
    use module_cu_nsas,  only: nsasinit, cu_nsas

    use options_interface,   only : options_t
    use domain_interface,    only : domain_t

    implicit none
    private
    public :: init_convection, convect, cu_var_request

    logical,allocatable, dimension(:,:) :: CU_ACT_FLAG
    real,   allocatable, dimension(:,:) :: XLAND, RAINCV, PRATEC, NCA
    real,   allocatable, dimension(:,:,:):: W0AVG, w_stochastic
    real,   allocatable, dimension(:,:)::  lowest_convection_layer, highest_convection_layer, hpbl  ! pbl height I guess? Used in NSAS. 

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

        else if (options%physics%convection == kCU_NSAS) then  ! Not checked thoroughly yet
            call options%alloc_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%temperature,                        &
                         kVARS%cloud_water, kVARS%cloud_ice, kVARS%precipitation, kVARS%convective_precipitation,   &
                         kVARS%exner, kVARS%dz_interface, kVARS%density, kVARS%pressure_interface, kVARS%pressure,  &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u, kVARS%v, kVARS%w, kVARS%land_mask,        &
                         kVARS%tend_qv, kVARS%tend_th, kVARS%tend_qc, kVARS%tend_qi, kVARS%tend_qs, kVARS%tend_qr,  &
                         kVARS%tend_u, kVARS%tend_v])

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
        integer :: i, j

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

            allocate(W0AVG(ims:ime,kms:kme,jms:jme))
            W0AVG = 0

            allocate(XLAND(ims:ime,jms:jme))
            XLAND = domain%land_mask
            where(domain%land_mask == 0) XLAND = 2 ! 0 is water if using "LANDMASK" as input

            allocate(w_stochastic(ims:ime,kms:kme,jms:jme))

        endif

        if (options%physics%convection == kCU_TIEDTKE) then
            if (this_image()==1) write(*,*) "    Tiedtke Cumulus scheme"

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
         !
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


        else if (options%physics%convection==kCU_NSAS) then
            if (this_image()==1) write(*,*) "     NSAS Cumulus scheme"

            allocate(hpbl(ims:ime,jms:jme))  ! ?? PBL height I assume? 
            allocate(lowest_convection_layer(ims:ime,jms:jme))
            allocate(highest_convection_layer(ims:ime,jms:jme))

            ! NSAS expects a PBL height, But since ICAR does not currently compute this, fix at rough esitmate of 1000m ?
            do i=ims,ime
                do j=jms,jme
                    hpbl(i,j)=1000.  ! or take terrain height into account?
                enddo
            enddo

            call nsasinit(  rthcuten=domain%tend%th                             &
                            ,rqvcuten=domain%tend%qv                            &
                            ,rqccuten=domain%tend%qc                            &
                            ,rqicuten=domain%tend%qi                            &
                            ,rucuten=domain%tend%u                              &
                            ,rvcuten=domain%tend%v                              &
                            ,restart=.false.             &  ! options%restart                            &
                            ,p_qc=1                                             & ! copied from Tiedke; no idea as to what these 3 flags represent. 
                            ,p_qi=1                                             &
                            ,p_first_scalar=0                                   &
                            ,allowed_to_read=.true.                             &
                            ,ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde   &
                            ,ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme   &
                            ,its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte &  ! -1 or not? 
                            )


         endif

    end subroutine init_convection

subroutine convect(domain,options,dt_in)
    implicit none
    type(domain_t),  intent(inout) :: domain
    type(options_t), intent(in)    :: options
    real, intent(in) :: dt_in

    integer :: j,itimestep,STEPCU, dx_factor_nsas
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

    ! Stochastic disturbance for  vertical speeds:
    if (options%cu_options%stochastic_cu /= kNO_STOCHASTIC) then
        call random_number(w_stochastic)
        block
            integer :: i
            do i=1,size(w_stochastic,2)
                w_stochastic(:,i,:) = domain%sensible_heat%data_2d/500 + w_stochastic(:,i,:)
            enddo
        end block
        W0AVG = domain%w_real%data_3d+(w_stochastic * options%cu_options%stochastic_cu - options%cu_options%stochastic_cu*0.75) ! e.g. * 20 - 15)
    else
        W0AVG = domain%w_real%data_3d
    endif


    if (options%physics%convection==kCU_TIEDTKE) then

        call CU_TIEDTKE(                                        &
                 dt_in, itimestep, STEPCU                       &
                ,RAINCV, PRATEC                                 &
                ,domain%latent_heat%data_2d/LH_vaporization     &
                ,domain%sensible_heat%data_2d                   &
                ,domain%znu                                     &
                ,domain%u_mass%data_3d                          &
                ,domain%v_mass%data_3d                          &
                ,W0AVG                                          &
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
    elseif (options%physics%convection==kCU_NSAS) then

        ! TO DO:
        !  - calculate height of PBL
        !  - use wrf_constants.f90 ?


        ! set dx_factor_nsas (cu_nsas.f90 line 567):
        if (domain%dx.le.1000) then
            dx_factor_nsas=1  !       if (dx_factor_nsas == 1) then   dx_factor = 250. / delx ! assume 2.5 ms-1 (1km) and 1.125 cms-1 (200km) (cu_nsas.f90 line 567)
        else
            dx_factor_nsas=2
        endif

        call cu_nsas( &
            dt=dt_in               & !-- dt          time step (s)
            ,dx=domain%dx           & !-- dx          grid interval (m)   ! till here
            ,p3di=domain%pressure_interface%data_3d              & !-- p3di        3d pressure (pa) at interface level
            ,p3d=domain%pressure%data_3d                         & !-- p3d         3d pressure (pa)
            ,pi3d=domain%exner%data_3d                           & !-- pi3d        3d exner function (dimensionless)
            ,qc3d=domain%cloud_water_mass%data_3d                & !-- qc3d        cloud water mixing ratio (kg/kg)
            ,qi3d=domain%cloud_ice_mass%data_3d                  & !-- qi3d        cloud ice mixing ratio (kg/kg)
            ,rho3d=domain%density%data_3d                        &
            ,itimestep=itimestep                                 &
            ,stepcu=STEPCU                                       &
            ,hbot=lowest_convection_layer                        & !-- HBOT          index of lowest model layer with convection  intent(out) ???
            ,htop=highest_convection_layer                       & !-- HTOP          index of highest model layer with convection
            ,cu_act_flag=CU_ACT_FLAG                             &
            ,rthcuten=domain%tend%th                             &
            ,rqvcuten=domain%tend%qv                             &
            ,rqccuten=domain%tend%qc                             &
            ,rqicuten=domain%tend%qi                             &
            ,rucuten=domain%tend%u                               &
            ,rvcuten=domain%tend%v                               &
            ,qv3d=domain%water_vapor%data_3d                     & !-- qv3d        3d water vapor mixing ratio (kg/kg)
            ,t3d=domain%temperature%data_3d                      & !-- t3d         temperature (k)
            ,raincv=RAINCV                                       & !-- raincv      cumulus scheme precipitation (mm)
            ,pratec=PRATEC                                       &
            ,xland=XLAND                                         &
            ,dz8w=domain%dz_interface%data_3d                    & !-- dz8w        dz between full levels (m)
            ,w=W0AVG                                             & !-- w           vertical velocity (m/s)
            ,u3d=domain%u_mass%data_3d                           & !-- u3d         3d u-velocity interpolated to theta points (m/s)
            ,v3d=domain%v_mass%data_3d                           & !-- v3d         3d v-velocity interpolated to theta points (m/s)
            ,hpbl=hpbl  &  ! ?? PBL height I assume? But since ICAR does not currently compute this, fix at rough esitmate of 1000m ?
            ,hfx=domain%sensible_heat%data_2d                     & !  HFX  - net upward heat flux at the surface (W/m^2)       
            ,qfx=domain%latent_heat%data_2d/LH_vaporization       & !  QFX  - net upward moisture flux at the surface (kg/m^2/s)
            ,mp_physics=5        & ! - sets ncloud: - integer no_cloud(0),no_ice(1),cloud+ice(2) (see ln 141 cu_nsas.f90)  use options%physics%microphysics?
            ,dx_factor_nsas=dx_factor_nsas                       & !
            ,p_qc=1                                              & ! copied from Tiedke; no idea as to what these 3 flags represent. 
            ,p_qi=1                                              &
            ,p_first_scalar=0                                    &
            ! ,pgcon=""                                            & ! pgcon_use  = 0.55  ! Zhang & Wu (2003,JAS) ! 0.55 is a physically-based value used by GFS
            ,cp=cp                                               & ! cp  = 1012.0    ! J/kg/K   specific heat capacity of moist STP air?
            ,cliq=4190.                                          & ! from wrf_constants.f90
            ,cpv=4.*461.6                                        & ! from wrf_constants.f90
            ,g=gravity                                           &
            ,xlv=2.5E6                                           & ! from wrf_constants.f90
            ,r_d=287.                                            & ! from wrf_constants.f90
            ,r_v=461.6                                           & ! from wrf_constants.f90
            ,ep_1=461.6/287.-1.                                  & ! EP_1=R_v/R_d-1. (wrf_constants.f90)
            ,ep_2=287./461.6                                     & ! EP_2=R_d/R_v
            ,cice=2106.                                           & ! from wrf_constants.f90                                          &
            ,xls=2.85E6                                           & ! from wrf_constants.f90
            ,psat=610.78                                           & ! from wrf_constants.f90
            ,f_qi=.true.                                             & ! optional arg
            ,f_qc=.true.                                             & ! optional arg
            ,ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,  &
            ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,  &
            its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte-1 )
        ! ! ! Definitions (for coupling - temporary):
        ! - - - - - - - - - - - - - - - - -
        ! microphysics scheme --> ncloud 
        !
        !      !! ncloud   - integer no_cloud(0),no_ice(1),cloud+ice(2) 
        !    if (mp_physics .eq. 0) then
        !     ncloud = 0
        !   elseif ( mp_physics .eq. 1 .or. mp_physics .eq. 3 ) then
        !     ncloud = 1
        !   else
        !     ncloud = 2
        !   endif
        ! - - - - - - - - - - - - - - - - -

    endif

    ! add domain wide tendency terms
    ! use a separate dt to make it easier to apply on a different dt
    internal_dt = dt_in

    if (options%physics%convection==kCU_TIEDTKE) then
        ! $omp parallel private(j) &
        ! $omp default(shared)
        ! $omp do schedule(static)
        do j=jts,jte
            if (options%cu_options%tendency_fraction > 0) then
                if (options%cu_options%tend_qv_fraction > 0) domain%water_vapor%data_3d(:,:,j)           = domain%water_vapor%data_3d(:,:,j)           + domain%tend%qv(:,:,j)*internal_dt * options%cu_options%tend_qv_fraction
                if (options%cu_options%tend_qc_fraction > 0) domain%cloud_water_mass%data_3d(:,:,j)      = domain%cloud_water_mass%data_3d(:,:,j)      + domain%tend%qc(:,:,j)*internal_dt * options%cu_options%tend_qc_fraction
                if (options%cu_options%tend_th_fraction > 0) domain%potential_temperature%data_3d(:,:,j) = domain%potential_temperature%data_3d(:,:,j) + domain%tend%th(:,:,j)*internal_dt * options%cu_options%tend_th_fraction
                if (options%cu_options%tend_qi_fraction > 0) domain%cloud_ice_mass%data_3d(:,:,j)        = domain%cloud_ice_mass%data_3d(:,:,j)        + domain%tend%qi(:,:,j)*internal_dt * options%cu_options%tend_qi_fraction
            endif
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
    endif

end subroutine convect
end module convection
