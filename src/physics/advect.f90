!> ----------------------------------------------------------------------------
!!  Standard advection scheme with variable order
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!! ----------------------------------------------------------------------------
module adv_std
    use data_structures
    use icar_constants
    use options_interface, only: options_t
    use domain_interface,  only: domain_t
    use adv_fluxcorr,      only: WRF_flux_corr
    implicit none
    private
    real,dimension(:,:,:),allocatable :: U_m, V_m, W_m, rho
    integer :: ims, ime, jms, jme, kms, kme, its, ite, jts, jte, i_s, i_e, j_s, j_e, horder, vorder
    real    :: dx

    ! For use advecting a (convective?) wind field
    ! real,dimension(:,:,:),allocatable :: U_4cu_u, V_4cu_u, W_4cu_u
    ! real,dimension(:,:,:),allocatable :: U_4cu_v, V_4cu_v, W_4cu_v

    public :: adv_std_init, adv_std_var_request, adv_std_advect3d, adv_std_compute_wind

contains

    subroutine adv_std_init(domain,options)
        implicit none
        type(domain_t), intent(in) :: domain
        type(options_t),intent(in) :: options


        ims = domain%ims
        ime = domain%ime
        jms = domain%jms
        jme = domain%jme
        kms = domain%kms
        kme = domain%kme
        its = domain%its
        ite = domain%ite
        jts = domain%jts
        jte = domain%jte
        
        !set order of advection
        horder = options%adv_options%h_order
        vorder = options%adv_options%v_order
        
        !Define bounds of advection computation. If using monotonic flux-limiter, it is necesarry to increase
        !advection bounds by 1. The necesarry extension of the halo is handeled in domain_object
        i_s = its
        i_e = ite
        j_s = jts
        j_e = jte
        
        if (options%adv_options%flux_corr==kFLUXCOR_MONO) then
            i_s = its - 1
            i_e = ite + 1
            j_s = jts - 1
            j_e = jte + 1
        endif
        
        ! if module level arrays are already allocated for some reason, deallocate them first
        if (allocated(U_m)) deallocate(U_m)
        if (allocated(V_m)) deallocate(V_m)
        if (allocated(W_m)) deallocate(W_m)
        !if (allocated(lastqv_m)) deallocate(lastqv_m)

        ! allocate the module level arrays
        allocate(U_m     (ims+1:ime,kms:kme,jms:jme  ))
        allocate(V_m     (ims:ime,  kms:kme,jms+1:jme))
        allocate(W_m     (ims:ime,  kms:kme,jms:jme  ))
        allocate(rho     (ims:ime,  kms:kme,jms:jme  ))
        !allocate(lastqv_m(ims:ime,  kms:kme,jms:jme  ))

        !     if (.not.allocated(U_4cu_u)) then
        !         allocate(U_4cu_u(nx,  nz, ny))
        !         U_4cu_u = 0
        !         allocate(V_4cu_u(nx+1,nz, ny-1))
        !         V_4cu_u = 0
        !         allocate(W_4cu_u(nx+1,nz, ny))
        !         W_4cu_u = 0
        !
        !         allocate(U_4cu_v(nx-1,nz, ny+1))
        !         U_4cu_v = 0
        !         allocate(V_4cu_v(nx,  nz, ny))
        !         V_4cu_v = 0
        !         allocate(W_4cu_v(nx,  nz, ny+1))
        !         W_4cu_v = 0
        !     endif
    end subroutine

    subroutine adv_std_var_request(options)
        implicit none
        type(options_t), intent(inout) :: options

        ! List the variables that are required to be allocated for adv4 advection
        call options%alloc_vars( &
                        [kVARS%u,    kVARS%v,   kVARS%w,     kVARS%dz_interface, kVARS%water_vapor])

        call options%advect_vars([kVARS%water_vapor, kVARS%potential_temperature])

        ! List the variables that are required for restarts with adv4 advection
        call options%restart_vars( &
                        [kVARS%u,    kVARS%v,   kVARS%w,     kVARS%dz_interface, kVARS%water_vapor])

    end subroutine

!     Note this routine has been manually inlined because the compiler didn't seem to optimize it well and made the array copies
!     the routine is left in for documentation.
!     subroutine flux2(l,r,U,nx,nz,ny,f)
!     !     Calculate the donor cell flux function
!     !     l = left gridcell scalar
!     !     r = right gridcell scalar
!     !     U = Courant number (u*dt/dx)
!     !
!     !     If U is positive, return l*U if U is negative return r*U
!     !     By using the mathematical form instead of the logical form,
!     !     we can run on the entire grid simultaneously, and avoid branches
!
!     !   arguments
!         implicit none
!         real, dimension(1:nx,1:nz,1:ny), intent(in) :: l,r,U
!         real, dimension(1:nx,1:nz,1:ny), intent(inout) :: f
!         integer,intent(in) :: ny,nz,nx
!         !   internal parameter
!         integer ::  err,i!,j,Ny,Nz,Nx
!         !   main code
!         f= ((U+ABS(U)) * l + (U-ABS(U)) * r)/2
!
!     end subroutine flux2

    subroutine flux3(q,u,v,w,flux_x,flux_z,flux_y)
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),    intent(in)          :: q
        real, dimension(i_s:i_e,  kms:kme,j_s:j_e), intent(in)     :: w
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e),  intent(in)    :: u
        real, dimension(i_s:i_e,  kms:kme,j_s:j_e+1),intent(in)    :: v
        
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e),intent(inout)      :: flux_x
        real, dimension(i_s:i_e,  kms:kme,j_s:j_e+1),intent(inout)    :: flux_y
        real, dimension(i_s:i_e,  kms:kme+1,j_s:j_e),intent(inout)    :: flux_z
        
        
        if (horder==1) then
            flux_x= ((u + ABS(u)) * q(i_s-1:i_e,:,j_s:j_e) + (u - ABS(u)) * q(i_s:i_e+1,:,j_s:j_e))  / 2
            flux_y= ((v + ABS(v)) * q(i_s:i_e,:,j_s-1:j_e) + (v - ABS(v)) * q(i_s:i_e,:,j_s:j_e+1))  / 2
        else if (horder==3) then
            !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
            flux_x = 7*(q(i_s:i_e+1,:,j_s:j_e)+q(i_s-1:i_e,:,j_s:j_e)) - &
                                (q(i_s+1:i_e+2,:,j_s:j_e)+q(i_s-2:i_e-1,:,j_s:j_e))
            flux_x = (u*flux_x)/12
            !Application of 3rd order diffusive terms
            flux_x = flux_x - (abs(u)/12) * (3*(q(i_s:i_e+1,:,j_s:j_e)-q(i_s-1:i_e,:,j_s:j_e)) - &
                                (q(i_s+1:i_e+2,:,j_s:j_e)-q(i_s-2:i_e-1,:,j_s:j_e)))
                            
            !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
            flux_y = 7*(q(i_s:i_e,:,j_s:j_e+1)+q(i_s:i_e,:,j_s-1:j_e)) - &
                                (q(i_s:i_e,:,j_s+1:j_e+2)+q(i_s:i_e,:,j_s-2:j_e-1))
            flux_y = (v*flux_y)/12
            !Application of 3rd order diffusive terms
            flux_y = flux_y - (abs(v)/12) * (3*(q(i_s:i_e,:,j_s:j_e+1)-q(i_s:i_e,:,j_s-1:j_e)) - &
                            (q(i_s:i_e,:,j_s+1:j_e+2)-q(i_s:i_e,:,j_s-2:j_e-1)))

        else if (horder==5) then
            !Calculation of 6th order fluxes for later application of 5th order diffusive terms
            flux_x = 37*(q(i_s:i_e+1,:,j_s:j_e)+q(i_s-1:i_e,:,j_s:j_e)) - &
                            8*(q(i_s+1:i_e+2,:,j_s:j_e)+q(i_s-2:i_e-1,:,j_s:j_e)) + &
                            (q(i_s+2:i_e+3,:,j_s:j_e)+q(i_s-3:i_e-2,:,j_s:j_e))
            flux_x = (u*flux_x)/60
            !Application of 5th order diffusive terms
            flux_x = flux_x - (abs(u)/60) * (10*(q(i_s:i_e+1,:,j_s:j_e)-q(i_s-1:i_e,:,j_s:j_e)) - &
                            5*(q(i_s+1:i_e+2,:,j_s:j_e)-q(i_s-2:i_e-1,:,j_s:j_e)) + &
                            (q(i_s+2:i_e+3,:,j_s:j_e)-q(i_s-3:i_e-2,:,j_s:j_e)))
                
            !Calculation of 6th order fluxes for later application of 5th order diffusive terms
            flux_y = 37*(q(i_s:i_e,:,j_s:j_e+1)+q(i_s:i_e,:,j_s-1:j_e)) - &
                            8*(q(i_s:i_e,:,j_s+1:j_e+2)+q(i_s:i_e,:,j_s-2:j_e-1)) + &
                            (q(i_s:i_e,:,j_s+2:j_e+3)+q(i_s:i_e,:,j_s-3:j_e-2))
            flux_y = (v*flux_y)/60
            !Application of 5th order diffusive terms
            flux_y = flux_y - (abs(v)/60) * (10*(q(i_s:i_e,:,j_s:j_e+1)-q(i_s:i_e,:,j_s-1:j_e)) - &
                            5*(q(i_s:i_e,:,j_s+1:j_e+2)-q(i_s:i_e,:,j_s-2:j_e-1)) + &
                            (q(i_s:i_e,:,j_s+2:j_e+3)-q(i_s:i_e,:,j_s-3:j_e-2)))
        endif
        
        if (vorder==1) then
            flux_z(:,kms+1:kme,:) = ((w(:,kms:kme-1,:) + ABS(w(:,kms:kme-1,:))) * q(i_s:i_e,kms:kme-1,j_s:j_e) + &
                                     (w(:,kms:kme-1,:) - ABS(w(:,kms:kme-1,:))) * q(i_s:i_e,kms+1:kme,j_s:j_e))  / 2
        else if (vorder==3) then
            !Calculation of 4th order fluxes for later application of 3rd order diffusive terms
            flux_z(:,kms+2:kme-1,:) = 7*(q(i_s:i_e,kms+2:kme-1,j_s:j_e)+q(i_s:i_e,kms+1:kme-2,j_s:j_e)) - &
                                    (q(i_s:i_e,kms+3:kme,j_s:j_e)+q(i_s:i_e,kms:kme-3,j_s:j_e))
            flux_z(:,kms+2:kme-1,:) = (w(:,kms+1:kme-2,:)*flux_z(:,kms+2:kme-1,:))/12
            !Application of 3rd order diffusive terms
            flux_z(:,kms+2:kme-1,:) = flux_z(:,kms+2:kme-1,:) - (abs(w(:,kms+1:kme-2,:))/12) * &
                        (3 * (q(i_s:i_e,kms+2:kme-1,j_s:j_e) - q(i_s:i_e,kms+1:kme-2,j_s:j_e)) - &
                            (q(i_s:i_e,kms+3:kme,j_s:j_e)   - q(i_s:i_e,kms:kme-3,j_s:j_e)))
                            
            !Do simple upwind for the cells who's stencil does not allow higher-order
            flux_z(:,kms+1,:) = ((w(:,kms,:) + ABS(w(:,kms,:))) * q(i_s:i_e,kms,j_s:j_e) + &
                                 (w(:,kms,:) - ABS(w(:,kms,:))) * q(i_s:i_e,kms+1,j_s:j_e))  / 2
            flux_z(:,kme,:) = ((w(:,kme-1,:) + ABS(w(:,kme-1,:))) * q(i_s:i_e,kme-1,j_s:j_e) + &
                               (w(:,kme-1,:) - ABS(w(:,kme-1,:))) * q(i_s:i_e,kme,j_s:j_e))  / 2

        else if (vorder==5) then
            !Calculation of 6th order fluxes for later application of 5th order diffusive terms
            flux_z(:,kms+3:kme-2,:) = 37*(q(i_s:i_e,kms+3:kme-2,j_s:j_e)+q(i_s:i_e,kms+2:kme-3,j_s:j_e)) - &
                                       8*(q(i_s:i_e,kms+4:kme-1,j_s:j_e)+q(i_s:i_e,kms+1:kme-4,j_s:j_e)) + &
                                         (q(i_s:i_e,kms+5:kme,j_s:j_e)+  q(i_s:i_e,kms:kme-5,j_s:j_e))
            flux_z(:,kms+3:kme-2,:) = (w(:,kms+2:kme-3,:)*flux_z(:,kms+3:kme-2,:))/60
            !Application of 5th order diffusive terms
            flux_z(:,kms+3:kme-2,:) = flux_z(:,kms+3:kme-2,:) - (abs(w(:,kms+2:kme-3,:))/60) * &
                        (10 * (q(i_s:i_e,kms+3:kme-2,j_s:j_e) - q(i_s:i_e,kms+2:kme-3,j_s:j_e)) - &
                            5*(q(i_s:i_e,kms+4:kme-1,j_s:j_e) - q(i_s:i_e,kms+1:kme-4,j_s:j_e)) + &
                              (q(i_s:i_e,kms+5:kme,j_s:j_e)   - q(i_s:i_e,kms:kme-5,j_s:j_e)))
                              
            flux_z(:,kms+2,:) = 7*(q(i_s:i_e,kms+2,j_s:j_e)+q(i_s:i_e,kms+1,j_s:j_e)) - &
                                    (q(i_s:i_e,kms+3,j_s:j_e)+q(i_s:i_e,kms,j_s:j_e))
            flux_z(:,kms+2,:) = (w(:,kms+1,:)*flux_z(:,kms+2,:))/12
            !Application of 3rd order diffusive terms
            flux_z(:,kms+2,:) = flux_z(:,kms+2,:) - (abs(w(:,kms+1,:))/12) * &
                        (3 * (q(i_s:i_e,kms+2,j_s:j_e) - q(i_s:i_e,kms+1,j_s:j_e)) - &
                            (q(i_s:i_e,kms+3,j_s:j_e)   - q(i_s:i_e,kms,j_s:j_e)))

            flux_z(:,kme-1,:) = 7*(q(i_s:i_e,kme-1,j_s:j_e)+q(i_s:i_e,kme-2,j_s:j_e)) - &
                                    (q(i_s:i_e,kme,j_s:j_e)+q(i_s:i_e,kme-3,j_s:j_e))
            flux_z(:,kme-1,:) = (w(:,kme-2,:)*flux_z(:,kme-1,:))/12
            !Application of 3rd order diffusive terms
            flux_z(:,kme-1,:) = flux_z(:,kme-1,:) - (abs(w(:,kme-2,:))/12) * &
                        (3 * (q(i_s:i_e,kme-1,j_s:j_e) - q(i_s:i_e,kme-2,j_s:j_e)) - &
                            (q(i_s:i_e,kme,j_s:j_e)   - q(i_s:i_e,kme-3,j_s:j_e)))
                            

            !Do simple upwind for the cells who's stencil does not allow higher-order
            flux_z(:,kms+1,:) = ((w(:,kms,:) + ABS(w(:,kms,:))) * q(i_s:i_e,kms,j_s:j_e) + &
                                 (w(:,kms,:) - ABS(w(:,kms,:))) * q(i_s:i_e,kms+1,j_s:j_e))  / 2
            flux_z(:,kme,:) = ((w(:,kme-1,:) + ABS(w(:,kme-1,:))) * q(i_s:i_e,kme-1,j_s:j_e) + &
                               (w(:,kme-1,:) - ABS(w(:,kme-1,:))) * q(i_s:i_e,kme,j_s:j_e))  / 2

        endif
        
                                                          

        !Handle top and bottom boundaries for z here
        flux_z(:,kms,:) = 0
        flux_z(:,kme+1,:) = q(i_s:i_e,kme,j_s:j_e) * w(:,kme,:)

                                         
    end subroutine flux3

    subroutine adv_std_advect3d(qfluxes,qold,dz,jaco,t_factor_in,flux_corr)

        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(inout)   :: qfluxes
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in)      :: qold
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in)      :: dz
        real, dimension(ims:ime,  kms:kme,jms:jme),  intent(in)      :: jaco
        real, optional,                              intent(in)      :: t_factor_in
        integer, optional,                           intent(in)      :: flux_corr
        ! interal parameters
        real, dimension(i_s:i_e+1,kms:kme,j_s:j_e)   :: flux_x
        real, dimension(i_s:i_e,  kms:kme,j_s:j_e+1) :: flux_y
        real, dimension(i_s:i_e,  kms:kme+1,j_s:j_e) :: flux_z
        real    :: t_factor
        
        ! !$omp parallel shared(qin,q,u,v,w) firstprivate(nx,ny,nz) private(i,f1,f3,f4,f5)
        ! !$omp do schedule(static)
        !do i=jms,jme
        !    q(:,:,i)=qin(:,:,i)
        !enddo
                
        ! !$omp end do
        ! !$omp barrier
        ! !$omp do schedule(static)
            ! by manually inlining the flux2 call we should remove extra array copies that the compiler doesn't remove.
            ! equivalent flux2 calls are left in for reference (commented) to restore recall that f1,f3,f4... arrays should be 3D : n x m x 1
            
        !Initialize t_factor, which is used during RK time stepping to scale the time step
        t_factor = 1.0
        if (present(t_factor_in)) t_factor = t_factor_in
        

        call flux3(qfluxes,U_m(i_s:i_e+1,:,j_s:j_e)*t_factor,&
                            V_m(i_s:i_e,:,j_s:j_e+1)*t_factor,&
                            W_m(i_s:i_e,:,j_s:j_e)*t_factor,flux_x,flux_z,flux_y)
        
        if (present(flux_corr)) then
            if (flux_corr == kFLUXCOR_MONO) then
                !Calculate flux corrections as done by WRF (Wang et al., 2009)
                call WRF_flux_corr(qold,U_m(i_s:i_e+1,:,j_s:j_e),&
                                        V_m(i_s:i_e,:,j_s:j_e+1),&
                                        W_m(i_s:i_e,:,j_s:j_e),flux_x,flux_z,flux_y, &
                                jaco(i_s:i_e,:,j_s:j_e),dz(i_s:i_e,:,j_s:j_e),rho(i_s:i_e,:,j_s:j_e))
            endif
        endif
        qfluxes = qold

        ! perform horizontal advection, from difference terms
        qfluxes(its:ite,:,jts:jte)  = qfluxes(its:ite,:,jts:jte)  - &
                                   ((flux_x(its+1:ite+1,:,jts:jte) - flux_x(its:ite,:,jts:jte)) + &
                                   (flux_y(its:ite,:,jts+1:jte+1) - flux_y(its:ite,:,jts:jte))) &
                                   / (jaco(its:ite,:,jts:jte)*rho(its:ite,:,jts:jte))                      
               ! then vertical (order doesn't matter because fluxes f1-6 are calculated before applying them)
               ! add fluxes to middle layers
        qfluxes(its:ite,:,jts:jte) = qfluxes(its:ite,:,jts:jte)  &
                                   - (flux_z(its:ite,kms+1:kme+1,jts:jte) - flux_z(its:ite,kms:kme,jts:jte)) &
                                   / (dz(its:ite,:,jts:jte)*jaco(its:ite,:,jts:jte)*rho(its:ite,:,jts:jte))

        ! !$omp end do
        ! !$omp end parallel
        
        !qfluxes(its:ite,:,jts:jte)  = qfluxes(its:ite,:,jts:jte)  + &
        !                           20*((qold(its+1:ite+1,:,jts:jte) - 2*qold(its:ite,:,jts:jte) + qold(its-1:ite-1,:,jts:jte)) + &
        !                           (qold(its:ite,:,jts+1:jte+1) - 2*qold(its:ite,:,jts:jte) + qold(its:ite,:,jts-1:jte-1))) &
        !                           / (4*dx)                      

        
    end subroutine adv_std_advect3d

    ! subroutine setup_cu_winds(domain, options, dt)
    !     implicit none
    !     type(domain_type),  intent(in) :: domain
    !     type(options_type), intent(in) :: options
    !     real,               intent(in) :: dt
    !
    !     real    :: dx
    !     integer :: nx,nz,ny
    !
    !     dx = domain%dx
    !     nx = size(domain%dz,1)
    !     nz = size(domain%dz,2)
    !     ny = size(domain%dz,3)
    !
    !
    !     U_4cu_u           =  (dt/dx) * (domain%u(1:nx,:,:)      + domain%u(2:nx+1,:,:)) / 2
    !     V_4cu_u(2:nx,:,:) =  (dt/dx) * (domain%v(1:nx-1,:,2:ny) + domain%v(2:nx,:,2:ny)) / 2
    !     W_4cu_u(2:nx,:,:) =  (dt/dx) * (domain%w(1:nx-1,:,:)    + domain%w(2:nx,:,:)) / 2
    !     call rebalance_cu_winds(U_4cu_u, V_4cu_u, W_4cu_u)
    !
    !     U_4cu_v(:,:,2:ny) =  (dt/dx) * (domain%u(2:nx,:,1:ny-1) + domain%u(2:nx,:,2:ny)) / 2
    !     V_4cu_v           =  (dt/dx) * (domain%v(:,:,1:ny)      + domain%v(:,:,2:ny+1)) / 2
    !     W_4cu_v(:,:,2:ny) =  (dt/dx) * (domain%w(:,:,1:ny-1)    + domain%w(:,:,2:ny)) / 2
    !     call rebalance_cu_winds(U_4cu_v, V_4cu_v, W_4cu_v)
    !
    ! end subroutine setup_cu_winds

    ! subroutine rebalance_cu_winds(u,v,w)
    !     implicit none
    !     ! u, v, w 3D east-west, south-north, and up-down winds repsectively
    !     ! note for this code, u is [nx-1,nz,ny] and v is [nx,nz,ny-1]
    !     real, dimension(:,:,:), intent(inout) :: u, v, w
    !
    !     real, allocatable, dimension(:,:) :: divergence, du, dv
    !     integer :: i,nx,ny,nz
    !
    !     nx = size(w,1)
    !     nz = size(w,2)
    !     ny = size(w,3)
    !
    !     allocate(divergence(nx-2,ny-2))
    !     allocate(du(nx-2,ny-2))
    !     allocate(dv(nx-2,ny-2))
    !
    !     do i=1,nz
    !         ! calculate horizontal divergence
    !         dv = v(2:nx-1,i,2:ny-1) - v(2:nx-1,i,1:ny-2)
    !         du = u(2:nx-1,i,2:ny-1) - u(1:nx-2,i,2:ny-1)
    !         divergence = du + dv
    !         ! Then calculate w to balance
    !         if (i==1) then
    !             ! if this is the first model level start from 0 at the ground
    !             w(2:nx-1,i,2:ny-1) = 0 - divergence
    !         else
    !             ! else calculate w as a change from w at the level below
    !             w(2:nx-1,i,2:ny-1) = w(2:nx-1,i-1,2:ny-1) - divergence
    !         endif
    !     enddo
    !
    ! end subroutine rebalance_cu_winds
    !
    ! subroutine advect_cu_winds(domain, options, dt)
    !     implicit none
    !     type(domain_type),  intent(inout) :: domain
    !     type(options_type), intent(in)    :: options
    !     real,               intent(in)    :: dt
    !
    !     integer :: nx,nz,ny
    !
    !     nx = size(domain%dz,1)
    !     nz = size(domain%dz,2)
    !     ny = size(domain%dz,3)
    !
    !     ! first put the background u,v,w winds on a staggered grid with respect to the u grid
    !     ! then advect the u winds
    !     if (options%advect_density) then
    !         print*, "ERROR: Density advection not enabled when using convective winds"
    !         print*, "   Requires update to wind.f90 balance_uvw and advect.f90 (at least)"
    !         stop
    !     endif
    !
    !     call setup_cu_winds(domain, options, dt)
    !
    !     ! set the top boundary condition for CU winds to 0 to prevent artifacts coming in from the "top"
    !     domain%u_cu(:,nz,:) = 0
    !     domain%v_cu(:,nz,:) = 0
    !
    !     call advect3d(domain%u_cu, U_4cu_u,V_4cu_u,W_4cu_u, domain%rho, domain%dz_inter, nx+1,nz,ny, options)
    !     call advect3d(domain%v_cu, U_4cu_v,V_4cu_v,W_4cu_v, domain%rho, domain%dz_inter, nx,nz,ny+1, options)
    !
    ! end subroutine advect_cu_winds


    subroutine test_divergence(dz)
        implicit none
        real, intent(in) :: dz(ims:ime,kms:kme,jms:jme)

        real, allocatable :: du(:,:), dv(:,:), dw(:,:)
        integer :: i,j,k

        allocate(du(ims+1:ime-1,jms+1:jme-1))
        allocate(dv(ims+1:ime-1,jms+1:jme-1))
        allocate(dw(ims+1:ime-1,jms+1:jme-1))

        do i=ims+1,ime-1
            do j=jms+1,jme-1
                do k=kms,kme
                    du(i,j) = (U_m(i+1,k,j)-U_m(i,k,j))
                    dv(i,j) = (V_m(i,k,j+1)-V_m(i,k,j))
                    if (k==kms) then
                        dw(i,j) = (W_m(i,k,j))/dz(i,k,j)
                    else
                        dw(i,j) = (W_m(i,k,j)-W_m(i,k-1,j))/dz(i,k,j)
                    endif
                    if (abs(du(i,j) + dv(i,j) + dw(i,j)) > 1e-3) then
                        print*, this_image(), i,k,j , abs(du(i,j) + dv(i,j) + dw(i,j))
                        print*, "Winds are not balanced on entry to advect"
                        !error stop
                    endif
                enddo
            enddo
        enddo

    end subroutine test_divergence

    subroutine adv_std_compute_wind(domain, options, dt)
        implicit none

        type(options_t),    intent(in)  :: options
        type(domain_t),  intent(inout) :: domain
        real,intent(in)::dt
        real, allocatable, dimension(:,:,:) :: rho_i

        
        ! if this if the first time we are called, we need to allocate the module level arrays
        ! Could/should be put in an init procedure
        if (.not.allocated(U_m)) then
            allocate(U_m     (ims+1:ime,kms:kme,jms:jme  ))
            allocate(V_m     (ims:ime,  kms:kme,jms+1:jme))
            allocate(W_m     (ims:ime,  kms:kme,jms:jme  ))
            !allocate(lastqv_m(ims:ime,  kms:kme,jms:jme  ))
        endif
        allocate(rho_i     (ims:ime,  kms:kme-1,jms:jme  ))


        ! if (options%physics%convection > 0) then
            ! print*, "Advection of convective winds not enabled in ICAR >=1.5 yet"
            ! stop
            ! U_m = (domain%u_cu(2:nx,:,:) + domain%u(2:nx,:,:)) * (dt/dx)
            ! V_m = (domain%v_cu(:,:,2:ny) + domain%v(:,:,2:ny)) * (dt/dx)
            ! W_m = (domain%w_cu + domain%w)                     * (dt/dx)
            ! call rebalance_cu_winds(U_m,V_m,W_m)
        ! else
             ! Divide only U and V by dx. This minimizes the number of operations per advection step. W cannot be divided by dz,
             ! since non-uniform dz spacing does not allow for the same spacing to be assumed on either side of a k+1/2 interface,
             ! as is required for the adv4 scheme.
            dx = domain%dx
            rho = 1
            if (options%parameters%advect_density) rho = domain%density%data_3d  
            
            rho_i(:,kms:kme-1,:) = ( rho(:,kms:kme-1,:)*domain%advection_dz(:,kms:kme-1,:) + rho(:,kms+1:kme,:)*domain%advection_dz(:,kms+1:kme,:) ) &
                                    / (domain%advection_dz(:,kms:kme-1,:)+domain%advection_dz(:,kms+1:kme,:))
        
            U_m = domain%u%data_3d(ims+1:ime,:,:) * dt * (rho(ims+1:ime,:,:)+rho(ims:ime-1,:,:))*0.5 * &
                    domain%jacobian_u(ims+1:ime,:,:) / domain%dx
            V_m = domain%v%data_3d(:,:,jms+1:jme) * dt * (rho(:,:,jms+1:jme)+rho(:,:,jms:jme-1))*0.5 * &
                    domain%jacobian_v(:,:,jms+1:jme) / domain%dx
                    
            W_m(:,kms:kme-1,:) = domain%w%data_3d(:,kms:kme-1,:) * dt * domain%jacobian_w(:,kms:kme-1,:) * rho_i
            W_m(:,kme,:) = domain%w%data_3d(:,kme,:) * dt * domain%jacobian_w(:,kme,:) * rho(:,kme,:)

    end subroutine adv_std_compute_wind

    subroutine setup_advection_dz(domain, options)
        implicit none
        type(domain_t),  intent(inout) :: domain
        type(options_t), intent(in)    :: options
        integer :: i

        if (.not.allocated(domain%advection_dz)) then
            allocate(domain%advection_dz(ims:ime,kms:kme,jms:jme))
        else
            return
        endif

        do i=kms,kme
            domain%advection_dz(:,i,:) = options%parameters%dz_levels(i)
        enddo

    end subroutine setup_advection_dz


end module adv_std
