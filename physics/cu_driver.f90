!> ----------------------------------------------------------------------------
!!
!!  Driver to call different convection schemes
!!
!!  Author: Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module convection
    use data_structures
    use module_cu_tiedtke
    implicit none
    real,allocatable,dimension(:,:,:)::p8,U3D,V3D,T3D, & 
                                       RTHCUTEN,RQVCUTEN,RQCCUTEN,RQICUTEN,    &
                                       RUCUTEN,RVCUTEN
    logical,allocatable,dimension(:,:)::CU_ACT_FLAG
    real,allocatable,dimension(:,:) ::XLAND, RAINCV, PRATEC
    real,allocatable,dimension(:) :: ZNU ! = (/0.9975, 0.9915, 0.984, 0.975, 0.965, 0.9525, 0.9375, &
                                        !    0.92, 0.9, 0.88, 0.86, 0.835, 0.805, 0.775, 0.745,   &
                                        !    0.71, 0.67, 0.63, 0.59, 0.55, 0.51, 0.47, 0.43, 0.39,&
                                        !    0.355, 0.325, 0.295, 0.27, 0.25, 0.23, 0.21, 0.19,   &
                                        !    0.17, 0.15, 0.13, 0.11, 0.09100001, 0.074, 0.059,    &
                                        !    0.046, 0.035, 0.025, 0.015, 0.005/)
                                        
contains
    subroutine znu_init(domain)
        type(domain_type), intent(in) :: domain
        integer::n_levels,i,xpt,ypt
        real::ptop,psfc
        
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
        type(domain_type), intent(in) :: domain
        type(options_type), intent(in) :: options
        integer::ids,ide,jds,jde,kds,kde
        
        write(*,*) "Initializing Cumulus Scheme"
        ids=1
        ide=size(domain%qv,1)
        kds=1
        kde=size(domain%qv,2)
        jds=1
        jde=size(domain%qv,3)

        call znu_init(domain)
        
        if (options%physics%convection==1) then
            write(*,*) "    Tiedtke Cumulus scheme"
            allocate(RTHCUTEN(ids:ide,kds:kde,jds:jde))
            RTHCUTEN=0
            allocate(RQVCUTEN(ids:ide,kds:kde,jds:jde))
            RQVCUTEN=0
            allocate(RQCCUTEN(ids:ide,kds:kde,jds:jde))
            RQCCUTEN=0
            allocate(RQICUTEN(ids:ide,kds:kde,jds:jde))
            RQICUTEN=0
            allocate(RUCUTEN(ids:ide,kds:kde,jds:jde))
            RUCUTEN=0
            allocate(RVCUTEN(ids:ide,kds:kde,jds:jde))
            RVCUTEN=0
            call tiedtkeinit(RTHCUTEN,RQVCUTEN,RQCCUTEN,RQICUTEN,   &
                                 RUCUTEN,RVCUTEN,                   &
                                 .false.,1,1,0,                     &
                                 .true.,                            &
                                 ids, ide, jds, jde, kds, kde-1,      &
                                 ids, ide, jds, jde, kds, kde,      &
                                 ids, ide, jds, jde, kds, kde-1)
         endif
        
    end subroutine init_convection
    
subroutine convect(domain,options,dt_in)
    implicit none
    type(domain_type), intent(inout) :: domain
    type(options_type), intent(in)   :: options
    real, intent(in) ::dt_in
    integer ::ids,ide,jds,jde,kds,kde,j,itimestep,STEPCU
    
    itimestep=1
    STEPCU=1
    ids=1
    ide=size(domain%qv,1)
    kds=1
    kde=size(domain%qv,2)
    jds=1
    jde=size(domain%qv,3)
    
    if (.not.allocated(p8)) then
!       allocate(pii(ids:ide,kds:kde,jds:jde))
!       pii=1
!       allocate(rho(ids:ide,kds:kde,jds:jde))
!       rho=1
!       allocate(zed1(ids:ide,kds:kde,jds:jde))
!       zed1=0
!       allocate(zed2(ids:ide,kds:kde,jds:jde))
!       zed2=0
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
        
        allocate(XLAND(ids:ide,jds:jde))
        allocate(RAINCV(ids:ide,jds:jde))
        allocate(PRATEC(ids:ide,jds:jde))
        
        where(domain%landmask==1) XLAND=1 ! land
        where(domain%landmask==0) XLAND=2 ! ocean (using LANDMASK as input)
        where(domain%landmask==2) XLAND=2 ! ocean (using XLAND as input)
        RAINCV=0
        PRATEC=0
    endif
    if (options%physics%convection==1) then

        !$omp parallel private(j) &
        !$omp firstprivate(ids,ide,jds,jde,kds,kde) &
        !$omp default(shared)
        !$omp do schedule(static)
        do j=jds,jde
            RAINCV(:,j)=0
            p8(:,:kde-1,j)=(domain%p(:,kds:kde-1,j)+domain%p(:,kds+1:kde,j))/2
            p8(:,kde,j)=p8(:,kde-1,j)+(p8(:,kde-1,j)-p8(:,kde-2,j))
            if (j>jds) then
                V3D(:,:,j)=(domain%v(:,:,j-1)+domain%v(:,:,j))/2
            endif
            U3D(:,:,j)=(domain%u(ids:ide,:,j)+domain%u(ids+1:ide+1,:,j))/2
            T3D(:,:,j)=domain%th(:,:,j)*domain%pii(:,:,j)
!           rho(:,:,j)=domain%p(:,:,j)/(R*T3D(:,:,j))
        enddo
        !$omp end do
        !$omp end parallel  

        call CU_TIEDTKE(                                          &
                 dt_in,itimestep,STEPCU                           &
                ,RAINCV,PRATEC,domain%latent_heat/LH_vaporization,domain%sensible_heat,ZNU(kds:kde) &
                ,U3D,V3D,domain%w,T3D,domain%qv,domain%cloud,domain%ice,domain%pii,domain%rho &
                ,domain%qv_adv_tendency,domain%qv_pbl_tendency    &
                ,domain%dz,p8,domain%p,XLAND,CU_ACT_FLAG          &
                ,ids+1,ide-1, jds+1,jde-1, kds,kde-1              &
                ,ids,ide, jds,jde, kds,kde                        &
                ,ids+1,ide-1, jds+1,jde-1, kds,kde-1              &
                ,RTHCUTEN,RQVCUTEN,RQCCUTEN,RQICUTEN              &
                ,RUCUTEN, RVCUTEN                                 &
                ,.True.,.True.,.True.,.True.,.True.               &
                )
!         write(*,*) "QV range"
!         write(*,*) MAXVAL(RQVCUTEN(2:ide-1,:,2:jde-1)), MINVAL(RQVCUTEN(2:ide-1,:,2:jde-1))
!         write(*,*) "qc range"
!         write(*,*) MAXVAL(RQCCUTEN(2:ide-1,:,2:jde-1)), MINVAL(RQCCUTEN(2:ide-1,:,2:jde-1))
!         write(*,*) "th range"
!         write(*,*) MAXVAL(RTHCUTEN(2:ide-1,:,2:jde-1)), MINVAL(RTHCUTEN(2:ide-1,:,2:jde-1))
!         write(*,*) "qi range"
!         write(*,*) MAXVAL(RQICUTEN(2:ide-1,:,2:jde-1)), MINVAL(RQICUTEN(2:ide-1,:,2:jde-1))
!         write(*,*) "rain range"
!         write(*,*) MAXVAL(RAINCV(2:ide-1,2:jde-1)), MINVAL(RAINCV(2:ide-1,2:jde-1))
        domain%qv=domain%qv+RQVCUTEN*dt_in
        domain%cloud=domain%cloud+RQCCUTEN*dt_in
        domain%th=domain%th+RTHCUTEN*dt_in
        domain%ice=domain%ice+RQICUTEN*dt_in
        domain%rain=domain%rain+RAINCV
        domain%crain=domain%crain+RAINCV
!       write(*,*) MAXVAL(domain%qv(2:ide-1,:,2:jde-1)), MINVAL(domain%qv(2:ide-1,:,2:jde-1))
    endif
    
    
end subroutine convect  
end module convection