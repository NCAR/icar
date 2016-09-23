!> ----------------------------------------------------------------------------
!!  Driver to call different convection schemes
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module convection
    use data_structures
    use module_cu_tiedtke,  only: tiedtkeinit, CU_TIEDTKE
    use module_cu_kf,       only: kfinit, KFCPS
    
    implicit none
    private
    public :: init_convection, convect
    logical,allocatable, dimension(:,:) :: CU_ACT_FLAG
    real,   allocatable, dimension(:,:) :: XLAND, RAINCV, PRATEC, NCA
    real,   allocatable, dimension(:,:,:):: W0AVG
    
contains

    subroutine init_convection(domain,options)
        implicit none
        type(domain_type),  intent(inout) :: domain
        type(options_type), intent(in) :: options
        integer :: ids,ide,jds,jde,kds,kde
        
        write(*,*) "Initializing Cumulus Scheme"
        ids=1
        ide=size(domain%qv,1)
        kds=1
        kde=size(domain%qv,2)
        jds=1
        jde=size(domain%qv,3)

        if (options%physics%convection>0) then
            if (.not.allocated(domain%tend%th)) allocate(domain%tend%th(ids:ide,kds:kde,jds:jde))
            domain%tend%th = 0
            if (.not.allocated(domain%tend%qv)) allocate(domain%tend%qv(ids:ide,kds:kde,jds:jde))
            domain%tend%qv = 0
            if (.not.allocated(domain%tend%qc)) allocate(domain%tend%qc(ids:ide,kds:kde,jds:jde))
            domain%tend%qc = 0
            if (.not.allocated(domain%tend%qr)) allocate(domain%tend%qr(ids:ide,kds:kde,jds:jde))
            domain%tend%qr = 0
            if (.not.allocated(domain%tend%qi)) allocate(domain%tend%qi(ids:ide,kds:kde,jds:jde))
            domain%tend%qi = 0
            if (.not.allocated(domain%tend%qs)) allocate(domain%tend%qs(ids:ide,kds:kde,jds:jde))
            domain%tend%qs = 0
            if (.not.allocated(domain%tend%u))  allocate(domain%tend%u(ids:ide,kds:kde,jds:jde))
            domain%tend%u  = 0
            if (.not.allocated(domain%tend%v))  allocate(domain%tend%v(ids:ide,kds:kde,jds:jde))
            domain%tend%v  = 0
        
            allocate(CU_ACT_FLAG(ids:ide,jds:jde))
            CU_ACT_FLAG=.True.
        
            allocate(RAINCV(ids:ide,jds:jde))
            allocate(PRATEC(ids:ide,jds:jde))
        
            RAINCV=0
            PRATEC=0
        endif
        
        if (options%physics%convection==kCU_TIEDTKE) then
            write(*,*) "    Tiedtke Cumulus scheme"
            if (.not.allocated(domain%tend%qv_adv)) then
                allocate(domain%tend%qv_adv(ids:ide,kds:kde,jds:jde))
                domain%tend%qv_adv=0
            endif
            if (.not.allocated(domain%tend%qv_pbl)) then
                allocate(domain%tend%qv_pbl(ids:ide,kds:kde,jds:jde))
                domain%tend%qv_pbl=0
            endif
            
            call tiedtkeinit(domain%tend%th,domain%tend%qv,   &
                             domain%tend%qc,domain%tend%qi,   &
                             domain%tend%u, domain%tend%v,    &
                             .false.,1,1,0,                   &
                             .true.,                          &
                             ids, ide, jds, jde, kds, kde-1,  &
                             ids, ide, jds, jde, kds, kde,    &
                             ids, ide, jds, jde, kds, kde-1)
         elseif (options%physics%convection==kCU_KAINFR) then
             write(*,*) "    Kain-Fritsch Cumulus scheme"
             allocate(W0AVG(ids:ide,kds:kde,jds:jde))
             W0AVG=0
             allocate(NCA(ids:ide,jds:jde))
             NCA=0
             call kfinit(domain%tend%th,domain%tend%qv,   &
                         domain%tend%qc,domain%tend%qr,   &
                         domain%tend%qi,domain%tend%qs,   &
                         NCA,W0AVG,1,1,                   & !...,P_QI,P_QS
                         0,.false.,.true.,                & ! P_FIRST_SCALAR, restart, allowed_to_read
                         ids, ide, jds, jde, kds, kde,  &
                         ids, ide, jds, jde, kds, kde,    &
                         ids, ide, jds, jde, kds, kde-1)
             
         endif
        
    end subroutine init_convection
    
subroutine convect(domain,options,dt_in)
    implicit none
    type(domain_type),  intent(inout) :: domain
    type(options_type), intent(in)    :: options
    real, intent(in) :: dt_in
    
    integer :: ids,ide,jds,jde,kds,kde,j,itimestep,STEPCU
    real :: internal_dt
    
    if (options%physics%convection==0) return

    itimestep=1
    STEPCU=1
    ids=1
    ide=size(domain%qv,1)
    kds=1
    kde=size(domain%qv,2)
    jds=1
    jde=size(domain%qv,3)
    
    if (.not.allocated(XLAND)) then
        ! if domain%landmask is defined before the call to init_convection, then this should be moved. 
        allocate(XLAND(ids:ide,jds:jde))
        where(domain%landmask==1) XLAND=1 ! land
        where(domain%landmask==0) XLAND=2 ! ocean (using LANDMASK as input)
        where(domain%landmask==2) XLAND=2 ! ocean (using XLAND as input)
    endif
    
    
    !$omp parallel private(j) &
    !$omp firstprivate(ids,ide,jds,jde,kds,kde) &
    !$omp default(shared)
    !$omp do schedule(static)
    do j=jds,jde
        RAINCV(:,j)     = 0
        domain%tend%Qv(:,:,j) = 0
        domain%tend%TH(:,:,j) = 0
        domain%tend%Qc(:,:,j) = 0
        domain%tend%Qi(:,:,j) = 0
        domain%tend%Qs(:,:,j) = 0
        domain%tend%Qr(:,:,j) = 0
        domain%tend%u(:,:,j)  = 0
        domain%tend%v(:,:,j)  = 0
    enddo
    !$omp end do
    !$omp end parallel
    
    
    if (options%physics%convection==kCU_TIEDTKE) then


        call CU_TIEDTKE(                                          &
                 dt_in,itimestep,STEPCU                           &
                ,RAINCV,PRATEC,domain%latent_heat/LH_vaporization &
                ,domain%sensible_heat,domain%ZNU                  &
                ,domain%Um,domain%Vm,domain%w_real                &
                ,domain%T,domain%qv                               &
                ,domain%cloud,domain%ice,domain%pii,domain%rho    &
                ,domain%tend%qv_adv,domain%tend%qv_pbl            &
                ,domain%dz_inter,domain%p,domain%p_inter,XLAND,CU_ACT_FLAG &
                ,ids+1,ide-1, jds+1,jde-1, kds,kde-1              & ! domain
                ,ids,ide,     jds,  jde,   kds,kde                & ! memory
                ,ids+1,ide-1, jds+1,jde-1, kds,kde-1              & ! tile
                ,domain%tend%th, domain%tend%qv, domain%tend%qc   &
                ,domain%tend%qi, domain%tend%u,  domain%tend%v    &
                ,.True.,.True.,.True.,.True.,.True.               &
                )
    
    elseif (options%physics%convection==kCU_KAINFR) then
        call KFCPS(                                          &
               ids,ide, jds,jde, kds,kde                     & ! domain
              ,ids,ide, jds,jde, kds,kde                     & ! memory
              ,ids+1,ide-1, jds+1,jde-1, kds,kde             & ! tile
              ,dt_in,itimestep,domain%dx,dt_in,.false.       & ! dt KTAU, dx, cu_dt, adapt_step
              ,domain%rho                                    & ! rho
              ,RAINCV,PRATEC,NCA                             & ! convective rain, convective rain rate, 
              ,domain%Um,domain%Vm,domain%th,domain%t        &
              ,domain%w_real                                 &
              ,domain%qv,domain%dz_inter,domain%p,domain%pii &
              ,W0AVG                                         & ! "average" vertical velocity (computed internally from w)
              ,XLV0,XLV1,XLS0,XLS1,cp,Rd,gravity,EP1         & ! physical "constants"
              ,EP2,SVP1,SVP2,SVP3,SVPT0                      & ! physical "constants"
              ,STEPCU,CU_ACT_FLAG,.false.                    & ! number of steps between CU calls, boolean grid to act on , warm_rain_only
            ! optional arguments
            ! ,F_QV    ,F_QC    ,F_QR    ,F_QI    ,F_QS      &
              ,.True., .True.,  .True.,  .True.,  .True.     &
              ,domain%tend%qv, domain%tend%qc                &
              ,domain%tend%qr, domain%tend%qi                &
              ,domain%tend%qs, domain%tend%th)
    endif
    
    ! add domain wide tendency terms
    ! use a separate dt to make it easier to apply on a different dt
    internal_dt=dt_in
!     $omp firstprivate(ids,ide, jds,jde, dt_in, internal_dt) &
    
    !$omp parallel private(j) &
    !$omp default(shared)
    !$omp do schedule(static)
    do j=jds,jde
        domain%qv(:,:,j)   =domain%qv(:,:,j)   + domain%tend%Qv(:,:,j)*internal_dt
        domain%cloud(:,:,j)=domain%cloud(:,:,j)+ domain%tend%Qc(:,:,j)*internal_dt
        domain%th(:,:,j)   =domain%th(:,:,j)   + domain%tend%TH(:,:,j)*internal_dt
        domain%ice(:,:,j)  =domain%ice(:,:,j)  + domain%tend%Qi(:,:,j)*internal_dt
        if (options%physics%convection==kCU_KAINFR) then
            domain%qsnow(:,:,j) =domain%qsnow(:,:,j) + domain%tend%Qs(:,:,j)*internal_dt
            domain%qrain(:,:,j) =domain%qrain(:,:,j) + domain%tend%Qr(:,:,j)*internal_dt
        endif
                
        domain%rain(:,j)   =domain%rain(:,j)  + RAINCV(:,j)
        domain%crain(:,j)  =domain%crain(:,j) + RAINCV(:,j)
        
        if (options%physics%convection==kCU_TIEDTKE) then
            domain%u(ids+1:ide,:,j) = domain%u(ids+1:ide,:,j) + (domain%tend%u(ids:ide-1,:,j)+domain%tend%u(ids+1:ide,:,j))/2 * internal_dt
            if (j>jds) then
                domain%v(:,:,j) = domain%v(:,:,j) + (domain%tend%v(:,:,j)+domain%tend%v(:,:,j-1))/2 * internal_dt
            endif
        endif
    enddo
    !$omp end do
    !$omp end parallel  
    
    
end subroutine convect  
end module convection
