!>----------------------------------------------------------
!! This module provides a wrapper to call various microphysics models
!! It sets up variables specific to the physics package to be used including
!! history variables not currently stored in the domain level data structure
!!
!! The main entry point to the code is mp(domain,options,dt)
!!
!! <pre>
!! Call tree graph :
!!  mp_init->[ external initialization routines]
!!  mp->[   external microphysics routines]
!!  mp_finish
!!
!! High level routine descriptions / purpose
!!   mp_init            - allocates module data and initializes physics package
!!   mp                 - sets up and calls main physics package
!!   mp_finish          - deallocates module memory, place to do the same for physics
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
module microphysics
    ! use data_structures
    use icar_constants
    use mod_wrf_constants
    use module_mp_thompson_aer,     only: mp_gt_driver_aer, thompson_aer_init
    use module_mp_thompson,         only: mp_gt_driver, thompson_init
    ! use module_mp_morr_two_moment,  only: MORR_TWO_MOMENT_INIT, MP_MORR_TWO_MOMENT
     use module_mp_wsm6,             only: wsm6, wsm6init
    use module_mp_simple,           only: mp_simple_driver, mp_simple_var_request
    use time_object,                only: Time_type
    use options_interface,          only: options_t
    use domain_interface,           only: domain_t

    implicit none

    ! permit the microphysics to update on a longer time step than the advection
    integer :: update_interval
    real*8 :: last_model_time
    ! temporary variables
    real,allocatable,dimension(:,:) :: SR, last_rain, last_snow, this_precip,refl_10cm


    ! microphysics specific flag.  If it returns the current hourly precip (e.g. Morrison), then set this to false.
    ! for use in "distributing" precipitation.
    logical :: precip_delta

    ! amounts of precipitation to "distribute" to surrounding grid cells
    integer, parameter :: npoints=8
    real,    dimension(npoints) :: dist_fraction = [ 0.1,0.15,0.1, 0.15,0.15, 0.1,0.15,0.1]
    integer, dimension(npoints) :: x_list = [ -1,0,1, -1,1, -1,0,1]
    integer, dimension(npoints) :: y_list = [ 1,1,1, 0,0, -1,-1,-1]

    public :: mp, mp_var_request
contains


    !>----------------------------------------------------------
    !! Initialize microphysical routines
    !!
    !! This routine will call the initialization routines for the specified microphysics packages.
    !! It also initializes any module level variables, e.g. update_interval
    !!
    !! @param   options     ICAR model options to specify required initializations
    !!
    !!----------------------------------------------------------
    subroutine mp_init(options)
        implicit none
        type(options_t), intent(inout) :: options

        if (this_image()==1) write(*,*) ""
        if (this_image()==1) write(*,*) "Initializing Microphysics"
        if (options%physics%microphysics    == kMP_THOMPSON) then
            if (this_image()==1) write(*,*) "    Thompson Microphysics"
            call thompson_init(options%mp_options)
            precip_delta=.True.

        elseif (options%physics%microphysics    == kMP_THOMP_AER) then
            if (this_image()==1) write(*,*) "    Thompson Eidhammer Microphysics"
            call thompson_aer_init()
            precip_delta=.True.

        elseif (options%physics%microphysics == kMP_SB04) then
            if (this_image()==1) write(*,*) "    Simple Microphysics"
            precip_delta=.True.
        ! elseif (options%physics%microphysics==kMP_MORRISON) then
        !     write(*,*) "    Morrison Microphysics"
        !     call MORR_TWO_MOMENT_INIT(hail_opt=0)
        !     precip_delta=.False.
        elseif (options%physics%microphysics==kMP_WSM6) then
            if (this_image()==1) write(*,*) "    WSM6 Microphysics"
            call wsm6init(rhoair0,rhowater,rhosnow,cliq,cpv)
        endif

        update_interval = options%mp_options%update_interval
        last_model_time = -999

    end subroutine mp_init


    subroutine mp_thompson_aer_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        ! List the variables that are required to be allocated for the simple microphysics
        call options%alloc_vars( &
                     [kVARS%pressure,    kVARS%potential_temperature,   kVARS%exner,        kVARS%density,      &
                      kVARS%water_vapor, kVARS%cloud_water,             kVARS%rain_in_air,  kVARS%rain_number_concentration, &
                      kVARS%snow_in_air, kVARS%cloud_ice,               kVARS%w,            kVARS%ice_number_concentration,      &
                      kVARS%snowfall,    kVARS%precipitation,           kVARS%graupel,      kVARS%graupel_in_air,     &
                      kVARS%dz, kVARS%re_cloud, kVARS%re_ice, kVARS%re_snow ])

        ! List the variables that are required to be advected for the simple microphysics
        call options%advect_vars( &
                      [kVARS%potential_temperature, kVARS%water_vapor, kVARS%cloud_water,  kVARS%rain_number_concentration, &
                       kVARS%snow_in_air,           kVARS%cloud_ice,   &
                       kVARS%rain_in_air,           kVARS%ice_number_concentration, kVARS%graupel_in_air   ] )

        ! List the variables that are required to be allocated for the simple microphysics
        call options%restart_vars( &
                       [kVARS%pressure,     kVARS%potential_temperature,    kVARS%water_vapor,   &
                        kVARS%cloud_water,  kVARS%rain_in_air,              kVARS%snow_in_air,   &
                        kVARS%precipitation,kVARS%snowfall,                 kVARS%graupel,       &
                        kVARS%dz,           kVARS%snow_in_air,              kVARS%cloud_ice,     &
                        kVARS%rain_number_concentration, kVARS%rain_in_air,                      &
                        kVARS%ice_number_concentration,  kVARS%graupel_in_air,                   &
                        kVARS%re_cloud, kVARS%re_ice, kVARS%re_snow  ] )


    end subroutine

    subroutine mp_wsm6_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        ! List the variables that are required to be allocated for the simple microphysics
        call options%alloc_vars( &
                     [kVARS%pressure,    kVARS%potential_temperature,   kVARS%exner,        kVARS%density,      &
                      kVARS%water_vapor, kVARS%cloud_water,             kVARS%rain_in_air,                      &
                      kVARS%snow_in_air, kVARS%cloud_ice,               kVARS%dz,                               &
                      kVARS%snowfall,    kVARS%precipitation,           kVARS%graupel,      kVARS%graupel_in_air ])

        ! List the variables that are required to be advected for the simple microphysics
        call options%advect_vars( &
                      [kVARS%potential_temperature, kVARS%water_vapor, kVARS%cloud_water,  &
                       kVARS%snow_in_air,           kVARS%cloud_ice,   &
                       kVARS%rain_in_air,           kVARS%graupel_in_air   ] )

        ! List the variables that are required to be allocated for the simple microphysics
        call options%restart_vars( &
                       [kVARS%pressure,     kVARS%potential_temperature,    kVARS%water_vapor,   &
                        kVARS%cloud_water,  kVARS%rain_in_air,              kVARS%snow_in_air,   &
                        kVARS%precipitation,kVARS%snowfall,                 kVARS%graupel,       &
                        kVARS%dz,           kVARS%snow_in_air,              kVARS%cloud_ice,     &
                        kVARS%rain_in_air,  kVARS%graupel_in_air ] )


    end subroutine


    subroutine mp_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        if (options%physics%microphysics    == kMP_THOMPSON) then
            call mp_thompson_aer_var_request(options)

        elseif (options%physics%microphysics    == kMP_THOMP_AER) then
            call mp_thompson_aer_var_request(options)

        elseif (options%physics%microphysics == kMP_SB04) then
            call mp_simple_var_request(options)

        elseif (options%physics%microphysics==kMP_MORRISON) then
            stop "Morrison physics not re-implemented yet"
        elseif (options%physics%microphysics==kMP_WSM6) then
            call mp_wsm6_var_request(options)
        endif

    end subroutine mp_var_request


    subroutine allocate_module_variables(ims,ime,jms,jme,kms,kme)
        implicit none
        integer, intent(in) :: ims,ime,jms,jme,kms,kme

        ! snow rain ratio
        if (.not.allocated(SR)) then
            allocate(SR(ims:ime,jms:jme))
            SR=0
        endif
        ! last snow amount
        if (.not.allocated(last_snow)) then
            allocate(last_snow(ims:ime,jms:jme))
            last_snow=0
        endif
        ! last rain amount
        if (.not.allocated(last_rain)) then
            allocate(last_rain(ims:ime,jms:jme))
            last_rain=0
        endif
        ! temporary precip amount
        if (.not.allocated(this_precip)) then
            allocate(this_precip(ims:ime,jms:jme))
            this_precip=0
        endif
        if (.not.allocated(refl_10cm)) then
            allocate(refl_10cm(ims:ime,jms:jme))
            refl_10cm=0
        endif
        ! if (.not.allocated(domain%tend%qr)) then
        !     allocate(domain%tend%qr(nx,nz,ny))
        !     domain%tend%qr=0
        ! endif
        ! if (.not.allocated(domain%tend%qs)) then
        !     allocate(domain%tend%qs(nx,nz,ny))
        !     domain%tend%qs=0
        ! endif
        ! if (.not.allocated(domain%tend%qi)) then
        !     allocate(domain%tend%qi(nx,nz,ny))
        !     domain%tend%qi=0
        ! endif


    end subroutine

    !>----------------------------------------------------------
    !! Distribute the microphysics precipitation to neighboring grid cells
    !!
    !! Because ICAR can be too aggressive at putting precip on mountain tops, this
    !! routine smooths out the precip by keeping only a fraction of it locally, and
    !! distributing the rest to the neighboring grid cells, weighted by distance.
    !!
    !! @param   [inout]current_precip   accumulated model precip at this time step
    !! @param   [in]last_precip         accumulated model precip prior to microphysics call
    !! @param   [in]local_fraction      fraction of precip to maintain in the local gridcell
    !! @param   [in]precip_delta        Flag to calculate this time steps precip as a delta from previous
    !!
    !!----------------------------------------------------------
    subroutine distribute_precip(current_precip, last_precip, local_fraction, precip_delta)
        implicit none
        real, dimension(:,:), intent(inout) :: current_precip, last_precip
        real, intent(in) :: local_fraction
        logical, intent(in) :: precip_delta
        ! relies on module variable this_precip as a temporary array

        integer :: i,j, nx,ny
        integer :: x,y, point

        nx=size(current_precip,1)
        ny=size(current_precip,2)

        if (precip_delta) then
            do j=2,ny-1
                do i=2,nx-1
                    this_precip(i,j)    = current_precip(i,j)-last_precip(i,j)
                    current_precip(i,j) = last_precip(i,j)
                end do
            end do
        else
            do j=2,ny-1
                do i=2,nx-1
                    this_precip(i,j)    = last_precip(i,j)
                    current_precip(i,j) = current_precip(i,j)-last_precip(i,j)
                end do
            end do
        endif

        do j=2,ny-1
            do i=2,nx-1
                current_precip(i,j) = current_precip(i,j)+this_precip(i,j)*local_fraction
                do point=1,npoints
                    x = i + x_list(point)
                    y = j + y_list(point)
                    current_precip(i,j) = current_precip(i,j) + this_precip(x,y) * (1-local_fraction) * dist_fraction(point)
                end do
            end do
        end do

    end subroutine distribute_precip

    !>----------------------------------------------------------
    !! Bias correct the precipitation
    !!
    !! Apply a pre-computed bias correction multiplier to the currently predicted precipitation
    !! rain_fraction should be 3D where the last dimension corresponds to some fractional distance
    !! through the year.
    !!
    !! @param   [in]current_date        Current date-time in the model simulation (used to compute time in rain_fraction)
    !! @param   [in]rain_fraction       Multiplier to apply to the current precipitation
    !! @param   [inout]current_precip   The current (potentially accumulated) predicted precipitation
    !! @param   [inout]last_precip      The last value of current_precip (in case it is accumulated)
    !! @param   [in]precip_delta        Flag that current_precip is an accumulation
    !!
    !!----------------------------------------------------------
    subroutine apply_rain_fraction(current_date, rain_fraction, current_precip, last_precip, precip_delta)
        implicit none
        type(Time_type),        intent(in)      :: current_date
        real, dimension(:,:,:), intent(in)      :: rain_fraction
        real, dimension(:,:),   intent(inout)   :: current_precip, last_precip
        logical,                intent(in)      :: precip_delta
        ! relies on module variable this_precip as a temporary array

        integer :: i,j, nx,ny
        integer :: x,y, point
        integer :: n_correction_points, correction_step
        real    :: year_fraction

        nx=size(current_precip,1)
        ny=size(current_precip,2)

        n_correction_points = size(rain_fraction,3)
        year_fraction = current_date%year_fraction()
        correction_step = min(  floor(n_correction_points * year_fraction)+1, &
                                n_correction_points)

        if (precip_delta) then
            do j=2,ny-1
                do i=2,nx-1
                    this_precip(i,j)    = current_precip(i,j)-last_precip(i,j)
                    current_precip(i,j) = last_precip(i,j)
                end do
            end do
        else
            ! In this case, last_precip is actually the current precip delta
            do j=2,ny-1
                do i=2,nx-1
                    this_precip(i,j)    = last_precip(i,j)
                    ! We have to remove the currently predicted precip from the accumulated precip
                    current_precip(i,j) = current_precip(i,j)-last_precip(i,j)
                end do
            end do
        endif

        do j=2,ny-1
            do i=2,nx-1
                current_precip(i,j) = current_precip(i,j)+this_precip(i,j)*rain_fraction(i,j,correction_step)
            end do
        end do

    end subroutine apply_rain_fraction


    subroutine process_subdomain(domain, options, dt,       &
                                its,ite, jts,jte, kts,kte,  &
                                ims,ime, jms,jme, kms,kme,  &
                                ids,ide, jds,jde, kds,kde)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real,           intent(in)    :: dt
        integer,        intent(in)    :: its,ite, jts,jte, kts,kte
        integer,        intent(in)    :: ims,ime, jms,jme, kms,kme
        integer,        intent(in)    :: ids,ide, jds,jde, kds,kde


        ! run the thompson microphysics
        if (options%physics%microphysics==kMP_THOMPSON) then
            ! call the thompson microphysics
            call mp_gt_driver(qv = domain%water_vapor%data_3d,                      &
                              th = domain%potential_temperature%data_3d,            &
                              qc = domain%cloud_water_mass%data_3d,                 &
                              qi = domain%cloud_ice_mass%data_3d,                   &
                              ni = domain%cloud_ice_number%data_3d,                 &
                              qr = domain%rain_mass%data_3d,                        &
                              nr = domain%rain_number%data_3d,                      &
                              qs = domain%snow_mass%data_3d,                        &
                              qg = domain%graupel_mass%data_3d,                     &
                              pii= domain%exner%data_3d,                            &
                              p =  domain%pressure%data_3d,                         &
                              dz = domain%dz_mass%data_3d,                          &
                              dt_in = dt,                                           &
                              itimestep = 1,                                        & ! not used in thompson
                              RAINNC = domain%accumulated_precipitation%data_2d,    &
                              RAINNCV = this_precip,                                & ! not used outside thompson (yet)
                              SR = SR,                                              & ! not used outside thompson (yet)
                              SNOWNC = domain%accumulated_snowfall%data_2d,         &
                              GRAUPELNC = domain%graupel%data_2d,       &
                              ids = ids, ide = ide,                   & ! domain dims
                              jds = jds, jde = jde,                   &
                              kds = kds, kde = kde,                   &
                              ims = ims, ime = ime,                   & ! memory dims
                              jms = jms, jme = jme,                   &
                              kms = kms, kme = kme,                   &
                              its = its, ite = ite,                   & ! tile dims
                              jts = jts, jte = jte,                   &
                              kts = kts, kte = kte)

        elseif (options%physics%microphysics==kMP_THOMP_AER) then
            ! call the thompson microphysics
            call mp_gt_driver_aer(qv = domain%water_vapor%data_3d,                      &
                                  th = domain%potential_temperature%data_3d,            &
                                  qc = domain%cloud_water_mass%data_3d,                 &
                                  qi = domain%cloud_ice_mass%data_3d,                   &
                                  ni = domain%cloud_ice_number%data_3d,                 &
                                  qr = domain%rain_mass%data_3d,                        &
                                  nr = domain%rain_number%data_3d,                      &
                                  qs = domain%snow_mass%data_3d,                        &
                                  qg = domain%graupel_mass%data_3d,                     &
                                  pii= domain%exner%data_3d,                            &
                                  p =  domain%pressure%data_3d,                         &
                                  w =  domain%w%data_3d,                                &
                                  dz = domain%dz_mass%data_3d,                          &
                                  dt_in = dt,                                           &
                                  RAINNC = domain%accumulated_precipitation%data_2d,    &
                                  SNOWNC = domain%accumulated_snowfall%data_2d,         &
                                  GRAUPELNC = domain%graupel%data_2d,       &
                                  re_cloud = domain%re_cloud%data_3d,                   &
                                  re_ice   = domain%re_ice%data_3d,                     &
                                  re_snow  = domain%re_snow%data_3d,                    &
                                  has_reqc=1, has_reqi=1, has_reqs=1,                   &
                                  ids = ids, ide = ide,                   & ! domain dims
                                  jds = jds, jde = jde,                   &
                                  kds = kds, kde = kde,                   &
                                  ims = ims, ime = ime,                   & ! memory dims
                                  jms = jms, jme = jme,                   &
                                  kms = kms, kme = kme,                   &
                                  its = its, ite = ite,                   & ! tile dims
                                  jts = jts, jte = jte,                   &
                                  kts = kts, kte = kte)

        elseif (options%physics%microphysics==kMP_SB04) then
            ! call the simple microphysics routine of SB04
            call mp_simple_driver(domain%pressure%data_3d,                  &
                                  domain%potential_temperature%data_3d,     &
                                  domain%exner%data_3d,                     &
                                  domain%density%data_3d,                   &
                                  domain%water_vapor%data_3d,               &
                                  domain%cloud_water_mass%data_3d,          &
                                  domain%rain_mass%data_3d,                 &
                                  domain%snow_mass%data_3d,                 &
                                  domain%accumulated_precipitation%data_2d, &
                                  domain%accumulated_snowfall%data_2d,      &
                                  dt,                                       &
                                  domain%dz_mass%data_3d,                   &
                                  ims = ims, ime = ime,                   & ! memory dims
                                  jms = jms, jme = jme,                   &
                                  kms = kms, kme = kme,                   &
                                  its = its, ite = ite,                   & ! tile dims
                                  jts = jts, jte = jte,                   &
                                  kts = kts, kte = kte)

        ! elseif (options%physics%microphysics==kMP_MORRISON) then
        !     call MP_MORR_TWO_MOMENT(itimestep,                         &
        !                     domain%th, domain%qv, domain%cloud,     &
        !                     domain%qrain, domain%ice, domain%qsnow, &
        !                     domain%qgrau, domain%nice, domain%nsnow,&
        !                     domain%nrain, domain%ngraupel,          &
        !                     domain%rho, domain%pii, domain%p,       &
        !                     mp_dt, domain%dz_inter, domain%w,       &
        !                     domain%rain, last_rain, SR,             &
        !                     domain%snow, last_snow, domain%graupel, &
        !                     this_precip,                            & ! hm added 7/13/13
        !                     refl_10cm, .False., 0,                  & ! GT added for reflectivity calcs
        !                     domain%tend%qr, domain%tend%qs,         &
        !                     domain%tend%qi                          &
        !                    ,IDS,IDE, JDS,JDE, KDS,KDE               & ! domain dims
        !                    ,IDS,IDE, JDS,JDE, KDS,KDE               & ! memory dims
        !                    ,ITS,ITE, JTS,JTE, KTS,KTE               & ! tile   dims            )
        !                )
        elseif (options%physics%microphysics==kMP_WSM6) then
            call wsm6(q = domain%water_vapor%data_3d,                      &
                              th = domain%potential_temperature%data_3d,            &
                              qc = domain%cloud_water_mass%data_3d,                 &
                              qi = domain%cloud_ice_mass%data_3d,                   &
                              qr = domain%rain_mass%data_3d,                        &
                              qs = domain%snow_mass%data_3d,                        &
                              qg = domain%graupel_mass%data_3d,                     &
                              pii= domain%exner%data_3d,                            &
                              p =  domain%pressure%data_3d,                         &
                              delz = domain%dz_mass%data_3d,                          &
                              den = domain%density%data_3d,                   &
                              delt = dt,                                           &
                              g = gravity,                                          &
                              cpd = cp, cpv = cpv, rd = Rd, rv = Rw, t0c = 273.15,          &
                              ep1 = EP1, ep2 = EP2, qmin = epsilon,                                &
                              XLS = XLS, XLV0 = XLV, XLF0 = XLF,                    &
                              den0 = rhoair0, denr = rhowater,                  &
                              cliq = cliq, cice = cice, psat = psat,                                   &
                              rain = domain%accumulated_precipitation%data_2d,    &
                              rainncv = this_precip,                                & ! not used outside thompson (yet)
                              sr = SR,                                              & ! not used outside thompson (yet)
                              snow = domain%accumulated_snowfall%data_2d,         &
                              graupel = domain%graupel%data_2d,       &
                              ids = ids, ide = ide,                   & ! domain dims
                              jds = jds, jde = jde,                   &
                              kds = kds, kde = kde,                   &
                              ims = ims, ime = ime,                   & ! memory dims
                              jms = jms, jme = jme,                   &
                              kms = kms, kme = kme,                   &
                              its = its, ite = ite,                   & ! tile dims
                              jts = jts, jte = jte,                   &
                              kts = kts, kte = kte)

        endif

        ! needs to be converted to work on specified tile or better, maybe moved out of microphysics driver entirely...
        ! if (options%use_bias_correction) then
        !     call apply_rain_fraction(domain%model_time, domain%rain_fraction, domain%rain, last_rain, precip_delta)
        !     call apply_rain_fraction(domain%model_time, domain%rain_fraction, domain%snow, last_snow, precip_delta)
        ! endif
        !
        ! if (options%mp_options%local_precip_fraction<1) then
        !     call distribute_precip(domain%rain, last_rain, options%mp_options%local_precip_fraction, precip_delta)
        !     call distribute_precip(domain%snow, last_snow, options%mp_options%local_precip_fraction, precip_delta)
        ! endif

    end subroutine process_subdomain

    subroutine process_halo(domain, options, dt, halo,  &
                            its,ite, jts,jte, kts,kte,  &
                            ims,ime, jms,jme, kms,kme,  &
                            ids,ide, jds,jde, kds,kde)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real,           intent(in)    :: dt
        integer,        intent(in)    :: halo
        integer,        intent(in)    :: its,ite, jts,jte, kts,kte
        integer,        intent(in)    :: ims,ime, jms,jme, kms,kme
        integer,        intent(in)    :: ids,ide, jds,jde, kds,kde
        integer :: halo_its,halo_ite, halo_jts,halo_jte, halo_kts,halo_kte

        ! process the western halo
        halo_ite = its+halo-1
        call process_subdomain(domain, options, dt, &
                    its,halo_ite, jts,jte, kts,kte, &
                    ims,ime, jms,jme, kms,kme,      &
                    ids,ide, jds,jde, kds,kde)


        ! process the eastern halo
        halo_its = ite-halo+1
        call process_subdomain(domain, options, dt, &
                    halo_its,ite, jts,jte, kts,kte, &
                    ims,ime, jms,jme, kms,kme,      &
                    ids,ide, jds,jde, kds,kde)

        ! for the top and bottom halos, we no longer process the corner elements, so subset i tile
        halo_its = its+halo
        halo_ite = ite-halo

        ! process the southern halo
        halo_jte = jts+halo-1
        call process_subdomain(domain, options, dt,           &
                    halo_its,halo_ite, jts,halo_jte, kts,kte, &
                    ims,ime, jms,jme, kms,kme,                &
                    ids,ide, jds,jde, kds,kde)


        ! process the northern halo
        halo_jts = jte-halo+1
        call process_subdomain(domain, options, dt,           &
                    halo_its,halo_ite, halo_jts,jte, kts,kte, &
                    ims,ime, jms,jme, kms,kme,                &
                    ids,ide, jds,jde, kds,kde)


    end subroutine process_halo


    !>----------------------------------------------------------
    !! Microphysical driver
    !!
    !! This routine handles calling the individual microphysics routine specified, that
    !! includes creating and passing any temporary variables, and checking when to update
    !! the microphysics based on the specified update_interval.
    !!
    !! @param   domain      ICAR model domain structure
    !! @param   options     ICAR model options structure
    !! @param   dt_in       Current driving time step (this is the advection step)
    !!
    !!----------------------------------------------------------
    subroutine mp(domain, options, dt_in, halo, subset)
        implicit none
        type(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        real,           intent(in)    :: dt_in
        integer,        intent(in),   optional :: halo, subset

        real :: mp_dt
        integer ::ids,ide,jds,jde,kds,kde, itimestep=1
        integer ::ims,ime,jms,jme,kms,kme
        integer ::its,ite,jts,jte,kts,kte, nx,ny,nz

        ids = domain%grid%ids;   ims = domain%grid%ims;   its = domain%grid%its
        ide = domain%grid%ide;   ime = domain%grid%ime;   ite = domain%grid%ite
        nx  = domain%grid%nx
        kds = domain%grid%kds;   kms = domain%grid%kms;   kts = domain%grid%kts
        kde = domain%grid%kde;   kme = domain%grid%kme;   kte = domain%grid%kte
        nz  = domain%grid%nz
        jds = domain%grid%jds;   jms = domain%grid%jms;   jts = domain%grid%jts
        jde = domain%grid%jde;   jme = domain%grid%jme;   jte = domain%grid%jte
        ny  = domain%grid%ny

        if (.not.allocated(SR)) call allocate_module_variables(ims, ime, jms, jme, kms, kme)

        ! if this is the first time mp is called, set last time such that mp will update
        if (last_model_time==-999) then
            last_model_time = (domain%model_time%seconds() - max(real(update_interval), dt_in))
        endif


        ! only run the microphysics if the next time step would put it over the update_interval time
        if (((domain%model_time%seconds() + dt_in)-last_model_time)>=update_interval) then

            ! calculate the actual time step for the microphysics
            mp_dt = domain%model_time%seconds()-last_model_time

            ! reset the counter so we know that *this* is the last time we've run the microphysics
            ! NOTE, ONLY reset this when running the inner subset... ideally probably need a separate counter for the halo and subset
            if (.not.present(halo)) then
                last_model_time = domain%model_time%seconds()
            endif

            ! If we are going to distribute the current precip over a few grid cells, we need to keep track of
            ! the last_precip so we know how much fell
            if ((options%mp_options%local_precip_fraction<1).or.(options%parameters%use_bias_correction)) then
                last_rain = domain%accumulated_precipitation%data_2d
                last_snow = domain%accumulated_snowfall%data_2d
            endif

            ! set the current tile to the top layer to process microphysics for
            if (options%mp_options%top_mp_level>0) then
                kte = min(kte, options%mp_options%top_mp_level)
            endif


            if (present(subset)) then
                call process_subdomain(domain, options, mp_dt,                 &
                                       its = its + subset, ite = ite - subset, &
                                       jts = jts + subset, jte = jte - subset, &
                                       kts = kts,          kte = kte,          &
                                       ims = ims, ime = ime,                   & ! memory dims
                                       jms = jms, jme = jme,                   &
                                       kms = kms, kme = kme,                   &
                                       ids = ids, ide = ide,                   & ! domain dims
                                       jds = jds, jde = jde,                   &
                                       kds = kds, kde = kde)
            endif

            if (present(halo)) then
                call process_halo(domain, options, mp_dt, halo, &
                                       its = its, ite = ite,    &
                                       jts = jts, jte = jte,    &
                                       kts = kts, kte = kte,    &
                                       ims = ims, ime = ime,    & ! memory dims
                                       jms = jms, jme = jme,    &
                                       kms = kms, kme = kme,    &
                                       ids = ids, ide = ide,    & ! domain dims
                                       jds = jds, jde = jde,    &
                                       kds = kds, kde = kde)

            endif

            if ((.not.present(halo)).and.(.not.present(subset))) then

                call process_subdomain(domain, options, mp_dt,  &
                                        its = its, ite = ite,    &
                                        jts = jts, jte = jte,    &
                                        kts = kts, kte = kte,    &
                                        ims = ims, ime = ime,    & ! memory dims
                                        jms = jms, jme = jme,    &
                                        kms = kms, kme = kme,    &
                                        ids = ids, ide = ide,    & ! domain dims
                                        jds = jds, jde = jde,    &
                                        kds = kds, kde = kde)

            endif

        endif

    end subroutine mp

    !>----------------------------------------------------------
    !! Finalize microphysical routines
    !!
    !! This routine will call the finalization routines (if any) for the specified microphysics packages.
    !! It also deallocates any module level variables, e.g. SR
    !!
    !! @param   options     ICAR model options to specify required initializations
    !!
    !!----------------------------------------------------------
    subroutine mp_finish(options)
        implicit none
        type(options_t),intent(in)::options

        if (allocated(SR)) then
            deallocate(SR)
        endif
        if (allocated(last_snow)) then
            deallocate(last_snow)
        endif
        if (allocated(last_rain)) then
            deallocate(last_rain)
        endif
        if (allocated(this_precip)) then
            deallocate(this_precip)
        endif
    end subroutine mp_finish
end module microphysics
