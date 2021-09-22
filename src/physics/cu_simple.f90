!>----------------------------------------------------------
!! Simple convection resolving convective parameterization code
!!
!! Because this scheme assume it is being run at convection resolving
!! scales, it does not modify the microphysical or cloud fraction variables
!! Instead the simple cu scheme modifies the wind field to create "resolved" updrafts
!! These in turn should lift air, leading to condensation, heating, and more lifting
!!
!! The initial version of this code does not use traditional convective indices (e.g. CAPE, CIN)
!! prefering to permit the physics resolve these effects outside of the convective scheme
!! This could almost be thought of more as a simple dynamics scheme then as a convective scheme
!!
!! The entry point to the code is cu_simple.
!!
!! High level routine descriptions / purpose
!!   cu_simple           - manages the entire module
!!   bouyancy            - computes a 3D bouyancy field
!!   wind_adjustment     - uses the bouyancy field to compute wind adjustments
!!
!! Driver inputs: pressure,th,pii,rho,qv,qc,qr,qs,rain,snow,dt,dz,nx,ny,nz
!!   pressure   = pressure                      - 3D - input  - Pa     - (nx,nz,ny)
!!   th         = potential temperature         - 3D - in/out - K      - (nx,nz,ny)
!!   pii        = exner function                - 3D - input  - []     - (nx,nz,ny)
!!   rho        = air density                   - 3D - input  - kg/m^3 - (nx,nz,ny)
!!   qv         = specific humidity             - 3D - in/out - kg/kg  - (nx,nz,ny)
!!   qc         = cloud water content           - 3D - in/out - kg/kg  - (nx,nz,ny)
!!   qr         = rain water content            - 3D - in/out - kg/kg  - (nx,nz,ny)
!!   qs         = snow water content            - 3D - in/out - kg/kg  - (nx,nz,ny)
!!   dt         = time step                     - 0D - input  - sec.   - scalar
!!   ims, ime   = start end of x array memory   - 0D - input  - n      - scalar
!!   jms, jme   = start end of y array memory   - 0D - input  - n      - scalar
!!   kms, kme   = start end of z array memory   - 0D - input  - n      - scalar
!!   its, ite   = start end of x tile to process- 0D - input  - n      - scalar
!!   jts, jte   = start end of y tile to process- 0D - input  - n      - scalar
!!   kts, kte   = start end of z tile to process- 0D - input  - n      - scalar
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module cu_simple_mod

    implicit none

    real,    parameter :: wind_effects(2) = [1, -1]
    integer, parameter :: u_x(2) = [[-1, 1], [-1, 1]]
    integer, parameter :: u_y(2) = [[ 0, 0], [ 0, 0]]
    integer, parameter :: v_y(2) = [[-1, 1], [-1, 1]]
    integer, parameter :: v_x(2) = [[ 0, 0], [ 0, 0]]

    real,    parameter :: time_const = 60.0 ! seconds

contains

    subroutine cu_simple(pressure,th,qv,pii,rho,u,v, dt,dz,    &
                         ims, ime, jms, jme, kms, kme, &
                         its, ite, jts, jte, kts, kte)
        implicit none
        real, intent(inout) :: pressure (ims:ime, kms:kme, jms:jme)
        real, intent(inout) :: th       (ims:ime, kms:kme, jms:jme)
        real, intent(inout) :: qv       (ims:ime, kms:kme, jms:jme)
        real, intent(inout) :: pii      (ims:ime, kms:kme, jms:jme)
        real, intent(inout) :: rho      (ims:ime, kms:kme, jms:jme)
        real, intent(inout) :: u        (ims:ime+1,kms:kme,jms:jme)
        real, intent(inout) :: v        (ims:ime, kms:kme, jms:jme+1)
        real, intent(inout) :: z        (ims:ime, kms:kme, jms:jme)
        real, intent(inout) :: rain     (ims:ime, jms:jme)
        real, intent(inout) :: snow     (ims:ime, jms:jme)
        real, intent(in)    :: dt
        integer,intent(in)  :: ims, ime, jms, jme, kms, kme
        integer,intent(in)  :: its, ite, jts, jte, kts, kte

        real    :: bouyancy (ims:ime,kms:kme,jms:jme)
        integer :: i, j, k

        ! $omp parallel default(shared)
        ! $omp private(i, j, k)
        ! $omp do
        do j = jts,jte
            do i = its, ite
                call calc_bouyancy(bouyancy(i,:,j), th(i,:,j), qv(i,:,j), z(i,:,j), kms,kme, kts,kte)

                call adjust_winds(bouyancy, u, v, i, j, dt,     &
                                  ims, ime, jms, jme, kms, kme, &
                                  kts, kte)
                enddo

            enddo
        enddo
        ! $omp end do
        ! $omp end parallel

    end subroutine cu_simple

    subroutine calc_bouyancy(bouyancy, th, qv, z, kms,kme, kts,kte)
        implicit none
        real,   intent(inout) :: bouyancy   (kms:kme)
        real,   intent(in)    :: th         (kms:kme)
        real,   intent(in)    :: qv         (kms:kme)
        real,   intent(in)    :: z          (kms:kme)
        integer,intent(in)    :: kms, kme, kts, kte

        integer :: k, k_top, k_bottom

        print*, "Calculate bouyancy"
        do k = kts, kte
            k_top    = min(k+1, kme)
            k_bottom = max(k-1, kms)
            bouyancy(k) = (th(k_top) - th(k_bottom)) / (z(k_top) - z(k_bottom))
        enddo

    end subroutine calc_bouyancy

end module
