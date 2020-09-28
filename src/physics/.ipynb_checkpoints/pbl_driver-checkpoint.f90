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
    ! use module_bl_ysu, only : ysuinit, ysu
    implicit none

    private
    public :: pbl_init, pbl, pbl_finalize

    integer :: ids, ide, jds, jde, kds, kde,  &
               ims, ime, jms, jme, kms, kme,  &
               its, ite, jts, jte, kts, kte

    logical :: allowed_to_read, restart, flag_qi

contains
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
        ! if (.not.allocated(domain%tend%qv_pbl)) allocate(domain%tend%qv_pbl(ims:ime,kms:kme,jms:jme))
        ! domain%tend%qv_pbl=0

        if (this_image()==1) write(*,*) "Initializing PBL Scheme"
        if (options%physics%boundarylayer==kPBL_SIMPLE) then
            if (this_image()==1) write(*,*) "    Simple PBL"
            call init_simple_pbl(domain, options)
        endif
        ! if (options%physics%boundarylayer==kPBL_YSU) then
        !     if (this_image()==1) write(*,*) "    YSU PBL"
        !     if (.not.allocated(domain%tend%th))     allocate(domain%tend%th(ims:ime,kms:kme,jms:jme))
        !     if (.not.allocated(domain%tend%qc))     allocate(domain%tend%qc(ims:ime,kms:kme,jms:jme))
        !     if (.not.allocated(domain%tend%qr))     allocate(domain%tend%qr(ims:ime,kms:kme,jms:jme))
        !     if (.not.allocated(domain%tend%qi))     allocate(domain%tend%qi(ims:ime,kms:kme,jms:jme))
        !     if (.not.allocated(domain%tend%u))      allocate(domain%tend%u(ims:ime,kms:kme,jms:jme))
        !     if (.not.allocated(domain%tend%v))      allocate(domain%tend%v(ims:ime,kms:kme,jms:jme))
        !     call ysuinit(domain%tend%u,domain%tend%v,       &
        !                  domain%tend%th,domain%tend%qv_pbl, &
        !                  domain%tend%qc,domain%tend%qi,1,1, &
        !                  restart, allowed_to_read,          &
        !                  ids, ide, jds, jde, kds, kde,      &
        !                  ims, ime, jms, jme, kms, kme,      &
        !                  its, ite, jts, jte, kts, kte)
        ! endif
    end subroutine pbl_init

    subroutine pbl(domain, options, dt_in)
        implicit none
        type(domain_t),  intent(inout)  :: domain
        type(options_t), intent(in)     :: options
        real,            intent(in)     :: dt_in

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
            stop "YSU PBL not implemented yet"
!             call ysu(domain%Um, domain%Vm,   domain%th, domain%t,               &
!                      domain%qv, domain%cloud,domain%ice,                        &
!                      domain%p,domain%p_inter,domain%pii,                        &
!                      domain%tend%u,domain%tend%v,domain%tend%th,                &
!                      domain%tend%qv_pbl,domain%tend%qc,domain%tend%qi,flag_qi,  &
!                      cp,gravity,rovcp,rd,rovg,                                  &
!                      domain%dz_i, domain%z,    LH_vaporization,rv,domain%psfc,  &
!                      domain%znu,  domain%znw,  domain%mut,domain%p_top,         &
!STARTHERE                      domain%znt,  domain%ustar,zol, hol, hpbl, psim, psih,      &
!                      domain%xland,domain%sensible_heat,domain%latent_heat,      &
!                      domain%tskin,gz1oz0,      wspd, br,                        &
!                      dt,dtmin,kpbl2d,                                           &
!                      svp1,svp2,svp3,svpt0,ep1,ep2,karman,eomeg,stbolt,          &
!                      exch_h,                                                    &
!                      domain%u10,domain%v10,                                     &
!                      ids,ide, jds,jde, kds,kde,                                 &
!                      ims,ime, jms,jme, kms,kme,                                 &
!                      its,ite, jts,jte, kts,kte)
        endif

    end subroutine pbl

    subroutine pbl_finalize(options)
        implicit none
        type(options_t), intent(in) :: options

        if (options%physics%boundarylayer==kPBL_SIMPLE) then
            call finalize_simple_pbl()
        endif

    end subroutine pbl_finalize
end module planetary_boundary_layer
