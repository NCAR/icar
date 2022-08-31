!>----------------------------------------------------------
!! This module provides a wrapper to call various PBL models
!! It sets up variables specific to the physics package to be used including both
!!
!! The main entry point to the code is pbl(domain,options,dt)
!!
!! <pre>
!! Call tree graph :
!!  pbl_init->[ external initialization routines]
!!  pbl->[  external PBL routines]
!!  pbl_finalize
!!
!! High level routine descriptions / purpose
!!   pbl_init           - initializes physics package
!!   pbl                - sets up and calls main physics package
!!   pbl_finalize       - permits physics package cleanup (close files, deallocate memory)
!!
!! Inputs: domain, options, dt
!!      domain,options  = as defined in data_structures
!!      dt              = time step (seconds)
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module planetary_boundary_layer
    use data_structures
    use domain_interface,   only : domain_t
    use options_interface,  only : options_t
    use pbl_simple,    only : simple_pbl, finalize_simple_pbl, init_simple_pbl
    use module_bl_ysu, only : ysuinit, ysu
    use mod_atm_utilities, only : calc_Richardson_nr
    use mod_wrf_constants, only : EOMEG
    use icar_constants !, only : karman,stefan_boltzmann
    use mod_pbl_utilities, only : da_sfc_wtq
    use ieee_arithmetic ! for debugging 
    use array_utilities, only : array_offset_x_3d, array_offset_y_3d


    implicit none
    real,allocatable, dimension(:,:)    ::  windspd, Ri, z_atm, zol, hol, hpbl, psim, &
                                            psih, u10d, v10d, t2d, q2d, gz1oz0, CHS, xland_real,regime
    ! integer, allocatable, dimension(:,:) :: kpbl2d
    real, allocatable, dimension(:,:,:) :: tend_u_ugrid, tend_v_vgrid

    private
    public :: pbl_init, pbl, pbl_finalize, pbl_var_request

    integer :: ids, ide, jds, jde, kds, kde,  &
               ims, ime, jms, jme, kms, kme,  &
               its, ite, jts, jte, kts, kte, j, k, i

    logical :: allowed_to_read, restart, flag_qi

contains

    subroutine pbl_var_request(options)
        implicit none
        type(options_t),intent(inout) :: options

        if (options%physics%boundarylayer==kPBL_YSU) then

            call options%alloc_vars( &
                         [kVARS%water_vapor, kVARS%potential_temperature, kVARS%temperature,                &
                         kVARS%exner, kVARS%dz_interface, kVARS%density, kVARS%pressure_interface,          &
                         kVARS%skin_temperature, kVARS%terrain, kVARS%ground_surf_temperature,              &
                         kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m,                  &
                         kVARS%humidity_2m, kVARS%surface_pressure, kVARS%ground_heat_flux,                 &
                         kVARS%znu, kVARS%znw, kVARS%roughness_z0, kVARS%ustar, kVARS%cloud_ice,            &
                         kVARS%tend_th_pbl, kVARS%tend_qc_pbl, kVARS%tend_qi_pbl,  kVARS%temperature_2m,    &
                         kVARS%tend_u, kVARS%tend_v, kVARS%tend_qv_pbl, kVARS%pressure, kVARS%kpbl,         &
                         kVARS%land_mask, kVARS%cloud_water, kVARS%coeff_heat_exchange_3d, kVARS%hpbl ]) !kVARS%tend_qv_adv,kVARS%tend_qv, kVARS%tend_qs, kVARS%tend_qr,, kVARS%u_mass, kVARS%v_mass,
!           kVARS%coeff_momentum_drag, ??
             call options%advect_vars([kVARS%potential_temperature, kVARS%water_vapor, kVARS%cloud_ice, kVARS%cloud_water]) !??

             call options%restart_vars( &
                        [kVARS%water_vapor, kVARS%potential_temperature, kVARS%temperature,                &
                        kVARS%exner, kVARS%dz_interface, kVARS%density, kVARS%pressure_interface,          &
                        kVARS%skin_temperature, kVARS%terrain, kVARS%ground_surf_temperature,              &
                        kVARS%sensible_heat, kVARS%latent_heat, kVARS%u_10m, kVARS%v_10m,                  &
                        kVARS%humidity_2m, kVARS%surface_pressure, kVARS%ground_heat_flux,                 &
                        kVARS%znu, kVARS%znw, kVARS%roughness_z0, kVARS%ustar, kVARS%cloud_ice,            &
                        kVARS%tend_th_pbl, kVARS%tend_qc_pbl, kVARS%tend_qi_pbl,  kVARS%temperature_2m,    &
                        kVARS%tend_u, kVARS%tend_v, kVARS%tend_qv_pbl, kVARS%pressure, kVARS%kpbl,         &
                        kVARS%land_mask, kVARS%cloud_water,kVARS%coeff_heat_exchange_3d, kVARS%hpbl  ]) !kVARS%u_mass, kVARS%v_mass,
        endif
    end subroutine pbl_var_request


    subroutine pbl_init(domain,options)
        implicit none
        type(domain_t),     intent(inout)   :: domain
        type(options_t),    intent(in)      :: options

        ids = domain%ids ; ide = domain%ide ; jds = domain%jds ; jde = domain%jde ; kds = domain%kds ; kde = domain%kde
        ims = domain%ims ; ime = domain%ime ; jms = domain%jms ; jme = domain%jme ; kms = domain%kms ; kme = domain%kme
        its = domain%its ; ite = domain%ite ; jts = domain%jts ; jte = domain%jte ; kts = domain%kts ; kte = domain%kte

        allowed_to_read = .True.
        restart = .False.
        flag_qi = .true.
        if (.not.allocated(domain%tend%qv_pbl)) allocate(domain%tend%qv_pbl(ims:ime,kms:kme,jms:jme))
        domain%tend%qv_pbl=0

        if (this_image()==1) write(*,*) "Initializing PBL Scheme"

        if (options%physics%boundarylayer==kPBL_SIMPLE) then
            if (this_image()==1) write(*,*) "    Simple PBL"
            call init_simple_pbl(domain, options)

        else if (options%physics%boundarylayer==kPBL_YSU) then

            if (this_image()==1) write(*,*) "    YSU PBL"

            ! allocate local vars YSU:
            allocate(windspd(ims:ime, jms:jme))
            allocate(Ri(ims:ime,jms:jme))
            ! Ri = 0
            allocate(z_atm(ims:ime,jms:jme))
            z_atm = domain%z%data_3d(:,kts,:) - domain%terrain%data_2d ! defines the height of the middle of the first model level
            allocate(zol(ims:ime, jms:jme)) ! zol		z/l height over monin-obukhov length - intent(inout) - but appears to not be used really?
            zol = 10
            allocate(hol(ims:ime, jms:jme)) ! hol		pbl height over monin-obukhov length - intent(inout)
            hol = 1000.0
            ! allocate(hpbl(ims:ime, jms:jme))  ! this should go to domain object for convective modules!!
            allocate(psim(ims:ime, jms:jme))
            ! psim= 0.5
            allocate(psih(ims:ime, jms:jme))
            ! psih=0.5
            allocate(u10d(ims:ime, jms:jme))
            allocate(v10d(ims:ime, jms:jme))
            allocate(t2d(ims:ime, jms:jme))
            allocate(q2d(ims:ime, jms:jme))
            allocate(gz1oz0(ims:ime, jms:jme))  !-- gz1oz0      log(z/z0) where z0 is roughness length
            gz1oz0 = log(z_atm / domain%roughness_z0%data_2d)
            ! allocate(kpbl2d(ims:ime, jms:jme)) ! domain%kpbl now
            ! allocate(CHS(ims:ime,jms:jme))
            ! CHS = 0.01
            allocate(xland_real(ims:ime,jms:jme))
            xland_real=real(domain%land_mask)
            allocate(regime(ims:ime,jms:jme))
            allocate(tend_u_ugrid(ims:ime+1, kms:kme, jms:jme)) ! to add the calculated u/v tendencies to the u/v grid
            allocate(tend_v_vgrid(ims:ime, kms:kme, jms:jme+1))

            ! initialize tendencies (this is done in ysu init but only for tiles, not mem (ie its vs ims))
            ! BK: check if this actually matters ???
            if(.not.restart)then
                do j = jms,jme
                do k = kms,kme
                do i = ims,ime
                    domain%tend%u(i,k,j) = 0.
                    domain%tend%v(i,k,j) = 0.
                    domain%tend%th_pbl(i,k,j) = 0.
                    domain%tend%qv_pbl(i,k,j) = 0.
                    domain%tend%qc_pbl(i,k,j) = 0.
                    domain%tend%qi_pbl(i,k,j) = 0.
                enddo
                enddo
                enddo
              endif


            call ysuinit(rublten=domain%tend%u                  &
                        ,rvblten=domain%tend%v                  &
                        ,rthblten=domain%tend%th_pbl            &
                        ,rqvblten=domain%tend%qv_pbl            &
                        ,rqcblten=domain%tend%qc_pbl            &
                        ,rqiblten=domain%tend%qi_pbl            &
                        ,p_qi=1                                 &
                        ,p_first_scalar=1                       &
                        ,restart=restart                        &
                        ,allowed_to_read= allowed_to_read      &
                        ,ids=ids, ide=ide, jds=jds, jde=jde     &
                        ,kds=kds, kde=kde, ims=ims, ime=ime     &
                        ,jms=jms, jme=jme, kms=kms, kme=kme     &
                        ,its=its, ite=ite, jts=jts, jte=jte     &
                        ,kts=kts, kte=kte-1)
        endif
    end subroutine pbl_init

    subroutine pbl(domain, options, dt_in)
        implicit none
        type(domain_t),  intent(inout)  :: domain
        type(options_t), intent(in)     :: options
        real,            intent(in)     :: dt_in  !  =real(dt%seconds())

        if (options%physics%boundarylayer==kPBL_SIMPLE) then
            call simple_pbl(domain% potential_temperature %data_3d,     &
                            domain% water_vapor           %data_3d,     &
                            domain% cloud_water_mass      %data_3d,     &
                            domain% cloud_ice_mass        %data_3d,     &
                            domain% rain_mass             %data_3d,     &
                            domain% snow_mass             %data_3d,     &
                            domain% u_mass                %data_3d,     &
                            domain% v_mass                %data_3d,     &
                            domain% exner                 %data_3d,     &
                            domain% density               %data_3d,     &
                            domain% z                     %data_3d,     &
                            domain% dz_mass               %data_3d,     &
                            domain% terrain               %data_2d,     &
                            its, ite, jts, jte, kts, kte,               &
                            dt_in)
                            ! domain% qv_pbl_tendency     %data_3d)
        endif

        if (options%physics%boundarylayer==kPBL_YSU) then

            ! windspd=sqrt(  domain%u_mass%data_3d(ims:ime, 1, jms:jme)**2 +     &
            !             domain%v_mass%data_3d(ims:ime, 1, jms:jme)**2   )
            windspd = sqrt(domain%u_10m%data_2d**2 + domain%v_10m%data_2d**2) ! as it is done in lsm_driver.
            where(windspd==0) windspd=1e-5

            ! Richardson number
            call calc_Richardson_nr(Ri,domain%temperature%data_3d, domain%skin_temperature%data_2d, z_atm, windspd)

            ! Copied from WRF, to calc psim and psih. ( Not 100% sure this is the way to go)
            call da_sfc_wtq ( psfc=domain%surface_pressure%data_2d           &
                            , tg=domain%ground_surf_temperature%data_2d     & !?
                            , ps=domain%pressure%data_3d(:,1,:)             &
                            , ts=domain%temperature%data_3d(:,1,:)                  &
                            , qs=domain%cloud_water_mass%data_3d(:,1,:)     &
                            , us=domain%u_mass%data_3d(:,1,:)               &
                            , vs=domain%v_mass%data_3d(:,1,:)               &
                            , hs=z_atm                                      & !: height at the lowest half sigma level
                            , roughness=domain%roughness_z0%data_2d &  !s(ims:ime, jms:jme)         &
                            , xland=xland_real                & ! real(domain%land_mask)
                            , dx=domain%dx                                  &
                            , u10=u10d, v10=v10d, t2=t2d, q2=q2d            & ! output only so can be dummies for now
                            , regime=regime                                 &
                            , psim=psim                                     & ! these we want
                            , psih=psih                                     &
                            , has_lsm=.true.                                & !if(options%physics%landsurface>1)
                            , ust_wrf=domain%ustar                          &
                            ! , regime_wrf, qsfc_wrf, znt_wrf, , mol_wrf, hfx, qfx, pblh  & ! optional
                            ,hfx=domain%sensible_heat%data_2d               &
                            ,qfx=domain%latent_heat%data_2d/LH_vaporization &
                            ,ims=ims, ime=ime, jms=jms, jme=jme)


            call ysu(u3d=domain%u_mass%data_3d                           & !-- u3d         3d u-velocity interpolated to theta points (m/s)
                    ,v3d=domain%v_mass%data_3d                           & !-- v3d         3d v-velocity interpolated to theta points (m/s)
                    ,th3d=domain%potential_temperature%data_3d           & 
                    ,t3d=domain%temperature%data_3d                      & 
                    ,qv3d=domain%water_vapor%data_3d                     & 
                    ,qc3d=domain%cloud_water_mass%data_3d                & !-- qc3d        cloud water mixing ratio (kg/kg)
                    ,qi3d=domain%cloud_ice_mass%data_3d                  & !-- qi3d        cloud ice mixing ratio (kg/kg)
                    ,p3d=domain%pressure%data_3d                         & !-- p3d         3d pressure (pa)
                    ,p3di=domain%pressure_interface%data_3d              & !-- p3di        3d pressure (pa) at interface level
                    ,pi3d=domain%exner%data_3d                           & !-- pi3d        3d exner function (dimensionless)
                    ,rublten=domain%tend%u                               & ! i/o
                    ,rvblten=domain%tend%v                  & ! i/o
                    ,rthblten=domain%tend%th_pbl            & ! i/o
                    ,rqvblten=domain%tend%qv_pbl            & ! i/o
                    ,rqcblten=domain%tend%qc_pbl            & ! i/o
                    ,rqiblten=domain%tend%qi_pbl            & ! i/o
                    ,flag_qi=.false.                        & ! not used in ysu code, so can be whatever?
                    ,cp=cp                                  &
                    ,g=gravity                              &
                    ,rovcp=rovcp                            & ! rovcp = Rd/cp
                    ,rd=Rd                                 &  ! J/(kg K) specific gas constant for dry air
                    ,rovg=rovg                              &
                    ,dz8w=domain%dz_interface%data_3d       & !-- dz8w        dz between full levels (m)
                    ,z=domain%z%data_3d                     & !-- z		height above sea level (m)
                    ,xlv=LH_vaporization                    & !-- xlv         latent heat of vaporization (j/kg)
                    ,rv=Rw                                  &  ! J/(kg K) specific gas constant for wet/moist air
                    ,psfc=domain%surface_pressure%data_2d   &
                    ,znu=domain%znu                         & ! znu and znw are only used if mut is provided. 
                    ,znw=domain%znw                         &
                !   ,mut=""  & ! optional - mass in a cell?
                !   ,p_top=""  & !,                                           && optional - only if mut is supplied
                    ,znt=domain%roughness_z0%data_2d       &  ! i/o -- znt		roughness length (m) (input only)
                    ,ust=domain%ustar                       & ! i/o -- ust		u* in similarity theory (m/s)
                    ,zol=zol                                & ! i/o -- zol		z/l height over monin-obukhov length - intent(inout) - but appears to not be used really?
                    ,hol=hol                                & ! i/o -- hol		pbl height over monin-obukhov length - intent(inout) 
                    ,hpbl=domain%hpbl%data_2d               & ! i/o -- hpbl	pbl height (m) - intent(inout)
                    ,psim=psim                              & !-- psim        similarity stability function for momentum - intent(in)
                    ,psih=psih                              & !-- psih        similarity stability function for heat- intent(in)
                    ,xland=real(domain%land_mask)                               &
                    ,hfx=domain%sensible_heat%data_2d                     & !  HFX  - net upward heat flux at the surface (W/m^2)
                    ,qfx=domain%latent_heat%data_2d/LH_vaporization       & !  QFX  - net upward moisture flux at the surface (kg/m^2/s)
                    ,tsk=domain%skin_temperature%data_2d                    &
                    ,gz1oz0=gz1oz0                          & !-- gz1oz0      log(z/z0) where z0 is roughness length
                    ,wspd=windspd                           & ! i/o -- wspd        wind speed at lowest model level (m/s) 
                    ,br=Ri                                  & !-- br          bulk richardson number in surface layer
                    ,dt=dt_in                               & !-- dt		time step (s)
                    ,dtmin=dt_in/60.                        & !-- dtmin	time step (minute)
                    ,kpbl2d=domain%kpbl                          & ! o --     ?? k layer of pbl top?? 
                    ,svp1=SVP1                              & !-- svp1        constant for saturation vapor pressure (kpa)
                    ,svp2=SVP2                              & !-- svp2        constant for saturation vapor pressure (dimensionless)
                    ,svp3=SVP3                              & !-- svp3        constant for saturation vapor pressure (k)
                    ,svpt0=SVPT0                            & !-- svpt0       constant for saturation vapor pressure (k)
                    ,ep1=EP1                                & !-- ep1         constant for virtual temperature (r_v/r_d - 1) (dimensionless)
                    ,ep2=EP2                                & !-- ep2         constant for specific humidity calculation
                    ,karman=karman                          & !-- karman      von karman constant
                    ,eomeg=EOMEG                            & !-- eomeg       angular velocity of earths rotation (rad/s)
                    ,stbolt=stefan_boltzmann                & !-- stbolt      stefan-boltzmann constant (w/m^2/k^4)
                ,exch_h=domain%coeff_heat_exchange_3d%data_3d  & ! i/o -- exch_h ! exchange coefficient for heat, K m/s , but 3d??
                    ,u10=domain%u_10m%data_2d               &
                    ,v10=domain%v_10m%data_2d               &
                    ,ids=ids, ide=ide, jds=jds, jde=jde     &
                    ,kds=kds, kde=kde, ims=ims, ime=ime     &
                    ,jms=jms, jme=jme, kms=kms, kme=kme     &
                    ,its=its, ite=ite, jts=jts, jte=jte     &
                    ,kts=kts, kte=kte-1                     &
                !optional
                    ,regime=regime                          )!  i/o -- regime	flag indicating pbl regime (stable, unstable, etc.) - not used?

                    ! if(this_image()==1 .and. options%parameters%debug) write(*,*) "  pbl height/lev is:", maxval(domain%hpbl%data_2d ),"m/", maxval(domain%kpbl)  ! uncomment if you want to see the pbl height. 

            !> ------------  add tendency terms  ------------
            !
            ! Here the tendency terms that were calculated by the ysu routine are added to the domain-wide fields.
            ! For u and v, we need to re-balance the uvw fields and re-compute dt after we change them. This is done in the
            ! step routine in time_step.f90, after the pbl call.
            !
            !> -----------------------------------------------

            ! Offset u/v tendencies to u and v grid, then add
            call array_offset_x_3d(domain%tend%u , tend_u_ugrid)
            call array_offset_y_3d(domain%tend%v , tend_v_vgrid)

            domain%u%data_3d   =  domain%u%data_3d  +  tend_u_ugrid  * dt_in
            domain%v%data_3d   =  domain%v%data_3d  +  tend_v_vgrid  * dt_in

            ! add mass grid tendencies
            domain%water_vapor%data_3d            =  domain%water_vapor%data_3d            + domain%tend%qv_pbl  * dt_in
            domain%cloud_water_mass%data_3d       =  domain%cloud_water_mass%data_3d       + domain%tend%qc_pbl  * dt_in
            domain%potential_temperature%data_3d  =  domain%potential_temperature%data_3d  + domain%tend%th_pbl  * dt_in
            domain%cloud_ice_mass%data_3d         =  domain%cloud_ice_mass%data_3d         + domain%tend%qi_pbl  * dt_in

            ! Reset tendencies before the next pbl call. (not sure if necessary)
            domain%tend%qv_pbl    = 0
            domain%tend%th_pbl    = 0
            domain%tend%qc_pbl    = 0
            domain%tend%qi_pbl    = 0
            domain%tend%u         = 0
            domain%tend%v         = 0



            ! -------------------- omp loop   - how to deal with offset (v) grid??   ---------------
            ! ! $omp parallel private(j) &
            ! ! $omp default(shared)
            ! ! $omp do schedule(static)
            ! do j=jts,jte ! OMP  loop

                ! domain%u%data_3d(:,:,j)  =  domain%u%data_3d(:,:,j) + tend_u_ugrid(:,:,j) * dt_in
                ! ! domain%v%data_3d(:,:,j)            =  domain%v%data_3d(:,:,j)       + domain%tend%v(:,:,j) * dt_in

                ! domain%water_vapor%data_3d(:,:,j)  =  domain%water_vapor%data_3d(:,:,j)  +  domain%tend%qv_pbl(:,:,j) * dt_in
                ! domain%cloud_water_mass%data_3d(:,:,j)      = domain%cloud_water_mass%data_3d(:,:,j)      + domain%tend%qc_pbl(:,:,j) * dt_in
                ! domain%potential_temperature%data_3d(:,:,j) = domain%potential_temperature%data_3d(:,:,j) + domain%tend%th_pbl(:,:,j) * dt_in
                ! domain%cloud_ice_mass%data_3d(:,:,j)        = domain%cloud_ice_mass%data_3d(:,:,j)        + domain%tend%qi_pbl(:,:,j) * dt_in

                ! ! Reset tendencies before the next pbl call. (necessary?)
                ! domain%tend%qv_pbl(:,:,j)   = 0
                ! domain%tend%th_pbl(:,:,j)   = 0
                ! domain%tend%qc_pbl(:,:,j)   = 0
                ! domain%tend%qi_pbl(:,:,j)   = 0

            ! enddo
            ! ! $omp end do
            ! ! $omp end parallel

        endif ! End YSU call

    end subroutine pbl

    subroutine pbl_finalize(options)
        implicit none
        type(options_t), intent(in) :: options

        if (options%physics%boundarylayer==kPBL_SIMPLE) then
            call finalize_simple_pbl()
        endif

    end subroutine pbl_finalize
end module planetary_boundary_layer
