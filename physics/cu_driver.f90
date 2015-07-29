!> ----------------------------------------------------------------------------
!!
!!  Driver to call different convection schemes
!!
!!  Author: Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module convection
    use data_structures
    use module_cu_tiedtke,  only: tiedtkeinit, CU_TIEDTKE
    use module_cu_kf,       only: kfinit, KFCPS
    
    implicit none
    real,   allocatable, dimension(:,:,:) :: p8,U3D,V3D,T3D, & 
                                             RTHCUTEN,RQVCUTEN,RQCCUTEN,RQICUTEN,    &
                                             RUCUTEN,RVCUTEN, RQRCUTEN, RQSCUTEN, W0AVG
    logical,allocatable, dimension(:,:) :: CU_ACT_FLAG
    real,   allocatable, dimension(:,:) :: XLAND, RAINCV, PRATEC, NCA
    real,   allocatable, dimension(:)   :: ZNU
    
    ! from wrf module_model_constants
    real, parameter :: XLV0         = 3.15E6
    real, parameter :: XLV1         = 2370.
    real, parameter :: XLS0         = 2.905E6
    real, parameter :: XLS1         = 259.532
    real, parameter ::  SVP1=0.6112
    real, parameter ::  SVP2=17.67
    real, parameter ::  SVP3=29.65
    real, parameter ::  SVPT0=273.15
    real, parameter ::  EP1=Rw/Rd-1.
    real, parameter ::  EP2=Rd/Rw
contains
    subroutine znu_init(domain)
        type(domain_type), intent(in) :: domain
        integer :: n_levels,i,xpt,ypt
        real    :: ptop,psfc
        
        n_levels=size(domain%p,2)
        
        xpt=2
        ypt=2
        ptop=domain%p(xpt,n_levels,ypt)-(domain%p(xpt,n_levels-1,ypt)-domain%p(xpt,n_levels,ypt))/2.0 !NOT CORRECT
        psfc=domain%p(xpt,1,ypt)+(domain%p(xpt,1,ypt)-domain%p(xpt,2,ypt))/2.0 !NOT CORRECT
        ptop=max(ptop,10.0)
        allocate(znu(n_levels))
        do i=1,n_levels
            znu(i)=(domain%p(xpt,i,ypt)-ptop)/(psfc-ptop)
        enddo
    end subroutine znu_init

    subroutine init_convection(domain,options)
        implicit none
        type(domain_type),  intent(in) :: domain
        type(options_type), intent(in) :: options
        integer :: ids,ide,jds,jde,kds,kde
        
        write(*,*) "Initializing Cumulus Scheme"
        ids=1
        ide=size(domain%qv,1)
        kds=1
        kde=size(domain%qv,2)
        jds=1
        jde=size(domain%qv,3)

        call znu_init(domain)
        
        allocate(RTHCUTEN(ids:ide,kds:kde,jds:jde))
        RTHCUTEN=0
        allocate(RQVCUTEN(ids:ide,kds:kde,jds:jde))
        RQVCUTEN=0
        allocate(RQCCUTEN(ids:ide,kds:kde,jds:jde))
        RQCCUTEN=0
        allocate(RQRCUTEN(ids:ide,kds:kde,jds:jde))
        RQRCUTEN=0
        allocate(RQICUTEN(ids:ide,kds:kde,jds:jde))
        RQICUTEN=0
        allocate(RQSCUTEN(ids:ide,kds:kde,jds:jde))
        RQSCUTEN=0
        allocate(RUCUTEN(ids:ide,kds:kde,jds:jde))
        RUCUTEN=0
        allocate(RVCUTEN(ids:ide,kds:kde,jds:jde))
        RVCUTEN=0
        
        allocate(p8(ids:ide,kds:kde,jds:jde))
        p8=0
        allocate(U3D(ids:ide,kds:kde,jds:jde))
        U3D=0
        allocate(V3D(ids:ide,kds:kde,jds:jde))
        V3D=0
        allocate(T3D(ids:ide,kds:kde,jds:jde))
        T3D=0
        allocate(CU_ACT_FLAG(ids:ide,jds:jde))
        CU_ACT_FLAG=.True.
        
        allocate(RAINCV(ids:ide,jds:jde))
        allocate(PRATEC(ids:ide,jds:jde))
        
        RAINCV=0
        PRATEC=0
        
        
        if (options%physics%convection==1) then
            write(*,*) "    Tiedtke Cumulus scheme"
            call tiedtkeinit(RTHCUTEN,RQVCUTEN,RQCCUTEN,RQICUTEN,   &
                                 RUCUTEN,RVCUTEN,                   &
                                 .false.,1,1,0,                     &
                                 .true.,                            &
                                 ids, ide, jds, jde, kds, kde-1,      &
                                 ids, ide, jds, jde, kds, kde,      &
                                 ids, ide, jds, jde, kds, kde-1)
         elseif (options%physics%convection==3) then
             write(*,*) "    Kain-Fritsch Cumulus scheme"
             allocate(W0AVG(ids:ide,kds:kde,jds:jde))
             allocate(NCA(ids:ide,jds:jde))
             NCA=-100
             call kfinit(RTHCUTEN,RQVCUTEN,RQCCUTEN,RQRCUTEN, &
                     RQICUTEN,RQSCUTEN,NCA,W0AVG,1,1,         & !...,P_QI,P_QS
                     0,.false.,.true.,                        & ! P_FIRST_SCALAR, restart, allowed_to_read
                     ids, ide, jds, jde, kds, kde-1,          &
                     ids, ide, jds, jde, kds, kde,            &
                     ids, ide, jds, jde, kds, kde-1                 )
             
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
        RQVCUTEN(:,:,j) = 0
        RTHCUTEN(:,:,j) = 0
        RQCCUTEN(:,:,j) = 0
        RQICUTEN(:,:,j) = 0
        RQSCUTEN(:,:,j) = 0
        RQRCUTEN(:,:,j) = 0
        p8(:,:kde-1,j)  = (domain%p(:,kds:kde-1,j) + domain%p(:,kds+1:kde,j))/2
        p8(:,kde,j)     =  domain%p(:,kde,j) + (domain%p(:,kde,j)-p8(:,kde-1,j))
        
        V3D(:,:,j)=(domain%v(:,:,j+1)+domain%v(:,:,j))/2
        U3D(:,:,j)=(domain%u(ids:ide,:,j)+domain%u(ids+1:ide+1,:,j))/2
        T3D(:,:,j)=domain%th(:,:,j)*domain%pii(:,:,j)
    enddo
    !$omp end do
    !$omp end parallel
    
    
    if (options%physics%convection==1) then


        call CU_TIEDTKE(                                          &
                 dt_in,itimestep,STEPCU                           &
                ,RAINCV,PRATEC,domain%latent_heat/LH_vaporization &
                ,domain%sensible_heat,ZNU(kds:kde),U3D,V3D        &
                ,domain%w*domain%dz/domain%dx,T3D,domain%qv       &
                ,domain%cloud,domain%ice,domain%pii,domain%rho    &
                ,domain%qv_adv_tendency,domain%qv_pbl_tendency    &
                ,domain%dz,p8,domain%p,XLAND,CU_ACT_FLAG          &
                ,ids+1,ide-1, jds+1,jde-1, kds,kde-1              & ! domain
                ,ids,ide, jds,jde, kds,kde                        & ! memory
                ,ids+1,ide-1, jds+1,jde-1, kds,kde-1              & ! tile
                ,RTHCUTEN,RQVCUTEN,RQCCUTEN,RQICUTEN              &
                ,RUCUTEN, RVCUTEN                                 &
                ,.True.,.True.,.True.,.True.,.True.               &
                )
    
    elseif (options%physics%convection==3) then
        call KFCPS(                                         &
               ids,ide, jds,jde, kds,kde                     & ! domain
              ,ids,ide, jds,jde, kds,kde                     & ! memory
              ,ids+1,ide-1, jds+1,jde-1, kds,kde             & ! tile
              ,dt_in,itimestep,domain%dx,dt_in,.false.       & ! dt KTAU, dx, cu_dt, adapt_step
              ,domain%rho                                    & ! rho
              ,RAINCV,PRATEC,NCA                             & ! convective rain, convective rain rate, 
              ,U3D,V3D,domain%th,T3D,domain%w*domain%dz/domain%dx      &
              ,domain%qv,domain%dz,domain%p,domain%pii       &
              ,W0AVG                                         & ! "average" vertical velocity (computed internally from w)
              ,XLV0,XLV1,XLS0,XLS1,CP,Rd,gravity,EP1         & ! physical "constants"
              ,EP2,SVP1,SVP2,SVP3,SVPT0                      & ! physical "constants"
              ,STEPCU,CU_ACT_FLAG,.false.                    & ! number of steps between CU calls, boolean grid to act on , warm_rain_only
            ! optional arguments
            ! ,F_QV    ,F_QC    ,F_QR    ,F_QI    ,F_QS      &
              ,.True., .True.,  .True.,  .True.,  .True.     &
              ,RQVCUTEN,RQCCUTEN,RQRCUTEN,RQICUTEN,RQSCUTEN  &
              ,RTHCUTEN                                      &
              )
    endif
    
    ! add domain wide tendency terms
    ! use a separate dt to make it easier to apply on a different dt
    internal_dt=dt_in
    
    !$omp parallel private(j) &
    !$omp firstprivate(ids,ide, jds,jde, dt_in, internal_dt) &
    !$omp default(shared)
    !$omp do schedule(static)
    do j=jds,jde
        domain%qv(:,:,j)   =domain%qv(:,:,j)   + RQVCUTEN(:,:,j)*internal_dt
        domain%cloud(:,:,j)=domain%cloud(:,:,j)+ RQCCUTEN(:,:,j)*internal_dt
        domain%th(:,:,j)   =domain%th(:,:,j)   + RTHCUTEN(:,:,j)*internal_dt
        domain%ice(:,:,j)  =domain%ice(:,:,j)  + RQICUTEN(:,:,j)*internal_dt
        if (options%physics%convection==3) then
            domain%qsnow(:,:,j) =domain%qsnow(:,:,j) + RQSCUTEN(:,:,j)*internal_dt
            domain%qrain(:,:,j) =domain%qrain(:,:,j) + RQRCUTEN(:,:,j)*internal_dt
        endif
                
        domain%rain(:,j)   =domain%rain(:,j)  + RAINCV(:,j)
        domain%crain(:,j)  =domain%crain(:,j) + RAINCV(:,j)
        
        if (options%physics%convection==1) then
            domain%u(ids+1:ide,:,j) = domain%u(ids+1:ide,:,j) + (RUCUTEN(ids:ide-1,:,j)+RUCUTEN(ids+1:ide,:,j))/2 * internal_dt
            if (j>jds) then
                domain%v(:,:,j) = domain%v(:,:,j) + (RVCUTEN(:,:,j)+RVCUTEN(:,:,j-1))/2 * internal_dt
            endif
        endif
    enddo
    !$omp end do
    !$omp end parallel  
    
    
end subroutine convect  
end module convection