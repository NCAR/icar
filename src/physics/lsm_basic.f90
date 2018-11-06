!>----------------------------------------------------------
!! Not really a land surface model, use coarse model fluxes
!!
!! This code was used in very early versions
!! of ICAR and is maintained for now just in case it becomes of interest again. 
!! Takes sensible and latent heat fluxes supplied by coarse model and adds them in
!! this code is generally deprecated. 
!!
!! Also provides simple PBL mixing, and radiative cooling 
!!  (which shouldn't be here.) 
!! 
!! The entry point to the code is lsm_basic(domain,options,dt)
!!
!! <pre>
!! Call tree graph :
!! lsm_basic->
!!  [->],
!!  [->],
!!  [->]
!! 
!! High level routine descriptions / purpose
!! 
!! Driver inputs: domain,options,dt
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module module_lsm_basic
    use data_structures

    implicit none
    integer,allocatable,dimension(:,:)::pblh_i
!     real, parameter :: LH_vaporization=2260000.0 ! J/kg
!     real, parameter :: R=287.058 ! J/(kg K) specific gas constant for air
!     real, parameter :: cp = 1012.0 ! specific heat capacity of moist STP air? J/kg/K
!     real, parameter :: g=9.81 ! gravity m/s^2
    contains

! add the latent and sensible heat fluxes to the temprature and humidity variables in the first model level
    subroutine simple_land_surface_fluxes(th,qv,p,dz,pii, dt, &
                                    sensible_heat, latent_heat, pblh, &
                                    ids,ide,kds,kde,jds,jde,options)
        implicit none
!       input potential temperature (K) and humidity mixing ratio (kg/kg)
        real, dimension(ids:ide,kds:kde,jds:jde),intent(inout) :: th,qv
!       input pressure (Pa), layer thickness (m), and potential temperature conversion function
        real, dimension(ids:ide,kds:kde,jds:jde),intent(in) :: p,dz,pii
!       input land surface fluxes (W/m^2)
        real, dimension(ids:ide,jds:jde),intent(in) :: sensible_heat, latent_heat
        integer, dimension(ids:ide,jds:jde),intent(in) :: pblh
!       input time step (seconds)
        real, intent(in) :: dt
!       array subscripts, shouldn't be necessary in fortran90, but something is getting tripped up
        integer, intent(in) :: ids,ide, jds,jde, kds,kde
        type(options_type),intent(in)::options
!       internal loop counter
        integer :: i,j
!       temporary variables internal to the loop and private to each thread
        real, dimension(ids:ide,kds:kde) :: rho
        real, dimension(ids+1:ide-1) :: temperature,dTemp,lhdQV
        real, dimension(kds:kde)::means,f1,f2,rhomean
        
        
!       not really needed anymore variables for tracking statistics and modifying fluxes
        real :: factor,mxdT,mxdqv,temp,diffusionrate,coolingrate,sbconstant
        
        temp=0.
        mxdT=0.
        mxdqv=0.
        factor=1.0 ! factor to scale land surface fluxes for tuning purposes
        ! diffusionrate*dq/dz must be < 1 (or even 0.5?) for model stability only if dt~>100
        diffusionrate=dt/(1.0*60.0*60.0)*1 ! *10 = rate "constant" this should be based on stability and wind speed?
        sbconstant=5.67e-8
        coolingrate=1.5*(dt/(60.0*60.0*24.0)) *sbconstant / 300.0 !1.5K/day radiative cooling rate (300 = W/m^2 at 270K)

        !       OpenMP commands
        !$omp parallel firstprivate(jds,jde,ids,ide,kds,kde,diffusionrate,coolingrate) &
        !$omp private(i,temperature,dTemp,lhdQV,rho,rhomean,means,f1,f2) &
        !$omp shared(pii,p,th,qv,dz,sensible_heat,latent_heat,pblh)
        !$omp do schedule(static)
        do i=jds+1,jde-1
!           loop over the last (slowest) dimension with the first (fastest) dimension processed internally
!           compute the real temperature from the potential temperature
            temperature=th(ids+1:ide-1,1,i)*pii(ids+1:ide-1,1,i) ! K
!           compute air density
            rho(ids+1:ide-1,:)=p(ids+1:ide-1,:,i)/(Rd*th(ids+1:ide-1,:,i)*pii(ids+1:ide-1,:,i)) ! kg/m^3
!           convert sensible heat flux to a temperature delta term
            ! J/s/m^2 * s / J/(kg*K) => kg*K/m^2 ... /((kg/m^3) * m) => K
            dTemp=(sensible_heat(ids+1:ide-1,i)*dt/cp) / (rho(ids+1:ide-1,1)*dz(ids+1:ide-1,1,i)) / factor! K
!           add temperature delta and convert back to potential temperature
            temperature=temperature+dTemp
            th(ids+1:ide-1,1,i)=temperature/pii(ids+1:ide-1,1,i)
!           convert latent heat flux to a mixing ratio tendancy term
            ! J/s/m^2 * s / J/kg => kg/m^2 ... / (kg/m^3 * m) => kg/kg
            lhdQV=(latent_heat(ids+1:ide-1,i)/LH_vaporization*dt) / (rho(ids+1:ide-1,1)*dz(ids+1:ide-1,1,i)) / factor
!           and add back to the mixing ratio
            qv(ids+1:ide-1,1,i)=qv(ids+1:ide-1,1,i)+lhdQV

            !WARNING: this could mask errors in other places in the model...
!           where(qv(ids+1:ide-1,1,i)<1e-10) qv(ids+1:ide-1,1,i)=1e-10 ! in case lhdQV in negative and abs(lh) > qv

!           stupid diffusion within the PBL I should at least speed up/clean up the math at some point...
            if (options%physics%boundarylayer==kPBL_BASIC) then
                do j=ids+1,ide-1
                    rhomean(kds:kde-1)=(rho(j,kds:kde-1)+rho(j,kds+1:kde))/2
                    rhomean(kde)=rho(j,kde)
                
                    means(kds:kde-1)=(qv(j,kds:kde-1,i)+qv(j,kds+1:kde,i))/2
    !               diffusion fluxes within the PBL
                    f1(1:pblh(j,i)-1)=(qv(j,1:pblh(j,i)-1,i)-means(1:pblh(j,i)-1))*diffusionrate
                    f2(1:pblh(j,i)-1)=(qv(j,2:pblh(j,i),i)-means(1:pblh(j,i)-1))*diffusionrate
    !               diffusion fluxes above the PBL
                    f1(pblh(j,i):kde-1)=(qv(j,pblh(j,i):kde-1,i)-means(pblh(j,i):kde-1))*diffusionrate/5
                    f2(pblh(j,i):kde-1)=(qv(j,pblh(j,i)+1:kde,i)-means(pblh(j,i):kde-1))*diffusionrate/5
    !               add fluxes to grid cells
                    qv(j,kds:kde-1,i)=qv(j,kds:kde-1,i)-f1(kds:kde-1)*rhomean(kds:kde-1)/rho(j,kds:kde-1)
                    qv(j,kds+1:kde,i)=qv(j,kds+1:kde,i)-f2(kds:kde-1)*rhomean(kds+1:kde)/rho(j,kds+1:kde)

    !               same process for potential temperature
                    means(kds:kde-1)=(th(j,kds:kde-1,i)+th(j,kds+1:kde,i))/2
                    f1(1:pblh(j,i)-1)=(th(j,1:pblh(j,i)-1,i)-means(1:pblh(j,i)-1))*diffusionrate
                    f2(1:pblh(j,i)-1)=(th(j,2:pblh(j,i),i)-means(1:pblh(j,i)-1))*diffusionrate
                    f1(pblh(j,i):kde-1)=(th(j,pblh(j,i):kde-1,i)-means(pblh(j,i):kde-1))*diffusionrate/5
                    f2(pblh(j,i):kde-1)=(th(j,pblh(j,i)+1:kde,i)-means(pblh(j,i):kde-1))*diffusionrate/5
                    th(j,kds:kde-1,i)=th(j,kds:kde-1,i)-f1(kds:kde-1)*rhomean(kds:kde-1)/rho(j,kds:kde-1)
                    th(j,kds+1:kde,i)=th(j,kds+1:kde,i)-f2(kds:kde-1)*rhomean(kds+1:kde)/rho(j,kds+1:kde)
                enddo
            endif
            if (options%physics%radiation==kRA_BASIC) then
!           stupid radiative cooling th=th-coolingrate*(T^4)
                th(ids+1:ide-1,:,i)=th(ids+1:ide-1,:,i) &
                            -(((th(ids+1:ide-1,:,i)*pii(ids+1:ide-1,:,i))**4)*coolingrate)
            endif
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine simple_land_surface_fluxes
    
    subroutine calc_pbl_index(z,pblh,pblh_i)
        implicit none
        real, dimension(:,:,:),intent(in)::z
        real, dimension(:,:),intent(in)::pblh
        integer,dimension(:,:),intent(out)::pblh_i
        integer::nx,ny,nz,i,j,k
        
        nx=size(z,1)
        nz=size(z,2)
        ny=size(z,3)
        pblh_i=1
        !$omp parallel firstprivate(nx,ny,nz) &
        !$omp private(i,j,k) &
        !$omp shared(pblh_i,pblh)
        !$omp do schedule(static)
        do j=1,ny
            do k=1,nz-1
                do i=1,nx
                    if (z(i,k,j)<=pblh(i,j)) then
                        pblh_i(i,j)=k
                    endif
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine calc_pbl_index
    
    subroutine lsm_basic(domain,options,dt)
        implicit none
        type(domain_type),intent(inout)::domain
        type(options_type),intent(in)::options
        real,intent(in)::dt
        integer::nx,ny,nz,j
        
        nx=size(domain%p,1)
        nz=size(domain%p,2)
        ny=size(domain%p,3)
        if (.not.allocated(pblh_i)) then
            allocate(pblh_i(nx,ny))
        endif
        call calc_pbl_index(domain%z,domain%pbl_height,pblh_i)
        
        call simple_land_surface_fluxes(domain%th,domain%qv,domain%p,domain%dz,domain%pii,dt, &
                                        domain%sensible_heat,domain%latent_heat,pblh_i, &
                                        1,nx,1,nz,1,ny,options)
        
    end subroutine lsm_basic
end module module_lsm_basic
