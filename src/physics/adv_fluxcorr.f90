!> ----------------------------------------------------------------------------
!!  A collection of flux correction schemes for advection
!!
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!! ----------------------------------------------------------------------------
module adv_fluxcorr
    use data_structures
    use domain_interface,  only: domain_t

    implicit none
    private
    integer :: ims, ime, jms, jme, kms, kme, its, ite, jts, jte

    public :: WRF_flux_corr, init_fluxcorr

contains

    subroutine init_fluxcorr(domain)
        implicit none
        type(domain_t), intent(in) :: domain


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
    end subroutine init_fluxcorr

    subroutine WRF_flux_corr(q,u,v,w,flux_x,flux_z,flux_y,jaco,dz,rho)
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),    intent(in) :: q
        real, dimension(its-1:ite+1,  kms:kme,jts-1:jte+1),  intent(in) :: w, jaco, dz, rho
        real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+1),  intent(in) :: u
        real, dimension(its-1:ite+1,  kms:kme,jts-1:jte+2),  intent(in) :: v
        
        real, dimension(its-1:ite+2,kms:kme,  jts-1:jte+1),intent(inout)          :: flux_x
        real, dimension(its-1:ite+1,kms:kme,  jts-1:jte+2),intent(inout)          :: flux_y
        real, dimension(its-1:ite+1,kms:kme+1,jts-1:jte+1),intent(inout)    :: flux_z
        
        
        real, dimension(its-1:ite+2,kms:kme,  jts-1:jte+1)         :: upwind_flux_x
        real, dimension(its-1:ite+1,kms:kme,  jts-1:jte+2)         :: upwind_flux_y
        real, dimension(its-1:ite+1,kms:kme+1,jts-1:jte+1)         :: upwind_flux_z
        
        real, dimension(its-1:ite+1,  kms:kme,jts-1:jte+1)   :: flux_in, scale_in, flux_out, scale_out, temp, qmax, qmin

        !Initialize some internal variables
        scale_in = 1.0
        scale_out = 1.0
              
        qmax = q(its-1:ite+1,:,jts-1:jte+1)
        qmin = q(its-1:ite+1,:,jts-1:jte+1)
        
        ! Get upwind fluxes
        call upwind_flux3(q,u,v,w,upwind_flux_x,upwind_flux_z,upwind_flux_y,jaco,dz,rho)

        !Compute concentration if upwind only was used
        temp  = q(its-1:ite+1,:,jts-1:jte+1) - ((upwind_flux_x(its:ite+2,:,:) - upwind_flux_x(its-1:ite+1,:,:)) + &
                                       (upwind_flux_y(:,:,jts:jte+2) - upwind_flux_y(:,:,jts-1:jte+1))) &
                                   / (jaco*rho)                      
        temp = temp - (upwind_flux_z(:,kms+1:kme+1,:) - upwind_flux_z(:,kms:kme,:)) / (dz*jaco*rho)
        
        ! We are now done with the input fluxes, so they can be overwritten to save memory
        ! flux now represents the merged fluxes
        flux_x = flux_x - upwind_flux_x
        flux_y = flux_y - upwind_flux_y
        flux_z = flux_z - upwind_flux_z
        
        ! Next compute max and min possible fluxes
        flux_in = -((min(0.,flux_x(its:ite+2,:,:)) - max(0.,flux_x(its-1:ite+1,:,:))) + &
                        (min(0.,flux_y(:,:,jts:jte+2)) - max(0.,flux_y(:,:,jts-1:jte+1)))) &
                                   / (jaco*rho)                         
        flux_in = flux_in - (min(0.,flux_z(:,kms+1:kme+1,:)) - max(0.,flux_z(:,kms:kme,:))) / (dz*jaco*rho)
        
        flux_out = ((max(0.,flux_x(its:ite+2,:,:)) - min(0.,flux_x(its-1:ite+1,:,:))) + &
                        (max(0.,flux_y(:,:,jts:jte+2)) - min(0.,flux_y(:,:,jts-1:jte+1)))) &
                                   / (jaco*rho)                         
        flux_out = flux_out + (max(0.,flux_z(:,kms+1:kme+1,:)) - min(0.,flux_z(:,kms:kme,:))) / (dz*jaco*rho)
        
                            
        where (u(its-1:ite+1,:,:) > 0)
            qmax(:,:,:) = max(q(its-2:ite,:,jts-1:jte+1),qmax(:,:,:))
            qmin(:,:,:) = min(q(its-2:ite,:,jts-1:jte+1),qmin(:,:,:))
        else where (u(its:ite+2,:,:) < 0)
            qmax(:,:,:) = max(q(its:ite+2,:,jts-1:jte+1),qmax(:,:,:))
            qmin(:,:,:) = min(q(its:ite+2,:,jts-1:jte+1),qmin(:,:,:))
        end where
        
        where (v(:,:,jts-1:jte+1) > 0)
            qmax(:,:,:) = max(q(its-1:ite+1,:,jts-2:jte),qmax(:,:,:))
            qmin(:,:,:) = min(q(its-1:ite+1,:,jts-2:jte),qmin(:,:,:))
        else where (v(:,:,jts:jte+2) < 0)
            qmax(:,:,:) = max(q(its-1:ite+1,:,jts:jte+2),qmax(:,:,:))
            qmin(:,:,:) = min(q(its-1:ite+1,:,jts:jte+2),qmin(:,:,:))
        end where
        
        where (w(:,kms:kme-1,:) > 0)
            qmax(:,kms+1:kme,:) = max(q(its-1:ite+1,kms:kme-1,jts-1:jte+1),qmax(:,kms+1:kme,:))
            qmin(:,kms+1:kme,:) = min(q(its-1:ite+1,kms:kme-1,jts-1:jte+1),qmin(:,kms+1:kme,:))
        else where (w(:,kms:kme-1,:) < 0)
            qmax(:,kms:kme-1,:) = max(q(its-1:ite+1,kms+1:kme,jts-1:jte+1),qmax(:,kms:kme-1,:))
            qmin(:,kms:kme-1,:) = min(q(its-1:ite+1,kms+1:kme,jts-1:jte+1),qmin(:,kms:kme-1,:))
        end where

        
        where(flux_in  > (qmax-temp))  scale_in = max(0.,(qmax-temp)/(flux_in + 1.E-15))
        where(flux_out > (temp-qmin)) scale_out = max(0.,(temp-qmin)/(flux_out+ 1.E-15))
              
        ! We are now done with merged fluxes, so they can be overwritten to save memory
        ! flux now represents the normalized fluxes
        where(flux_x(its:ite+1,:,:) > 0)
            flux_x(its:ite+1,:,:) = min(scale_in(its:ite+1,:,:),scale_out(its-1:ite,:,:))*flux_x(its:ite+1,:,:)
        else where(flux_x(its:ite+1,:,:) < 0)
            flux_x(its:ite+1,:,:) = min(scale_out(its:ite+1,:,:),scale_in(its-1:ite,:,:))*flux_x(its:ite+1,:,:)
        end where

        where(flux_y(:,:,jts:jte+1) > 0)
            flux_y(:,:,jts:jte+1) = min(scale_in(:,:,jts:jte+1),scale_out(:,:,jts-1:jte))*flux_y(:,:,jts:jte+1)
        else where(flux_y(:,:,jts:jte+1) < 0)
            flux_y(:,:,jts:jte+1) = min(scale_out(:,:,jts:jte+1),scale_in(:,:,jts-1:jte))*flux_y(:,:,jts:jte+1)
        end where
        
        where(flux_z(:,kms+1:kme,:) > 0)
            flux_z(:,kms+1:kme,:) = min(scale_in(:,kms+1:kme,:),scale_out(:,kms:kme-1,:))*flux_z(:,kms+1:kme,:)
        else where(flux_z(:,kms+1:kme,:) < 0)
            flux_z(:,kms+1:kme,:) = min(scale_out(:,kms+1:kme,:),scale_in(:,kms:kme-1,:))*flux_z(:,kms+1:kme,:)
        end where

        !Finally, compute the output fluxes
        flux_x = flux_x + upwind_flux_x
        flux_y = flux_y + upwind_flux_y
        flux_z = flux_z + upwind_flux_z
        
    end subroutine WRF_flux_corr

    subroutine upwind_flux3(q,u,v,w,flux_x,flux_z,flux_y,jaco,dz,rho)
        implicit none
        real, dimension(ims:ime,  kms:kme,jms:jme),    intent(in)      :: q
        real, dimension(its-1:ite+2,  kms:kme,jts-1:jte+1),  intent(in)    :: u
        real, dimension(its-1:ite+1,  kms:kme,jts-1:jte+2),  intent(in)    :: v
        real, dimension(its-1:ite+1,  kms:kme,jts-1:jte+1),  intent(in) :: w, jaco, dz, rho

        real, dimension(its-1:ite+2,kms:kme,jts-1:jte+1),intent(inout)          :: flux_x
        real, dimension(its-1:ite+1,kms:kme,jts-1:jte+2),intent(inout)          :: flux_y
        real, dimension(its-1:ite+1,kms:kme+1,jts-1:jte+1),intent(inout)    :: flux_z
        
        real, dimension(ims:ime,  kms:kme,jms:jme) :: dumb_q

        !When using RK3, we may have a time step derived using a CFL constraint larger than 1
        !This means that our upwind advection here may be in violation of the CFL criterion,
        !since this does not take place within the RK3 scheme. Since the possible CFL constraints
        !under RK3 with the available advection orders are all < 2, we can simply do 2 upwind steps

        !Compute upwind fluxes for first step
        flux_x= 0.5*((u + ABS(u)) * q(its-2:ite+1,:,jts-1:jte+1)  + (u - ABS(u)) * q(its-1:ite+2,:,jts-1:jte+1))  / 2

        flux_y= 0.5*((v + ABS(v)) * q(its-1:ite+1,:,jts-2:jte+1) +  (v - ABS(v)) * q(its-1:ite+1,:,jts-1:jte+2))  / 2

        flux_z(:,kms+1:kme,:) = 0.5*((w(:,kms:kme-1,:) + ABS(w(:,kms:kme-1,:))) * q(its-1:ite+1,kms:kme-1,jts-1:jte+1) + &
                                 (w(:,kms:kme-1,:) - ABS(w(:,kms:kme-1,:))) * q(its-1:ite+1,kms+1:kme,jts-1:jte+1))  / 2
                                         
        !Handle top and bottom boundaries for z here
        flux_z(:,kms,:) = 0
        flux_z(:,kme+1,:) = 0.5*q(its-1:ite+1,kme,jts-1:jte+1) * w(:,kme,:)

        !Update intermediate concentration
        dumb_q = q
        dumb_q(its-1:ite+1,:,jts-1:jte+1)  = q(its-1:ite+1,:,jts-1:jte+1) - ((flux_x(its:ite+2,:,:) - flux_x(its-1:ite+1,:,:)) + &
                                       (flux_y(:,:,jts:jte+2) - flux_y(:,:,jts-1:jte+1))) &
                                   / (jaco*rho)                      
        dumb_q(its-1:ite+1,:,jts-1:jte+1) = dumb_q(its-1:ite+1,:,jts-1:jte+1) - (flux_z(:,kms+1:kme+1,:) - flux_z(:,kms:kme,:)) / (dz*jaco*rho)
                        
        !Now compute upwind fluxes after second step
        flux_x= flux_x + 0.5*((u + ABS(u)) * dumb_q(its-2:ite+1,:,jts-1:jte+1)  + (u - ABS(u)) * dumb_q(its-1:ite+2,:,jts-1:jte+1))  / 2

        flux_y= flux_y + 0.5*((v + ABS(v)) * dumb_q(its-1:ite+1,:,jts-2:jte+1) +  (v - ABS(v)) * dumb_q(its-1:ite+1,:,jts-1:jte+2))  / 2

        flux_z(:,kms+1:kme,:) = flux_z(:,kms+1:kme,:) + 0.5*((w(:,kms:kme-1,:) + ABS(w(:,kms:kme-1,:))) * dumb_q(its-1:ite+1,kms:kme-1,jts-1:jte+1) + &
                                 (w(:,kms:kme-1,:) - ABS(w(:,kms:kme-1,:))) * dumb_q(its-1:ite+1,kms+1:kme,jts-1:jte+1))  / 2
                                         
        !Handle top flux again
        flux_z(:,kme+1,:) = flux_z(:,kme+1,:) + 0.5*dumb_q(its-1:ite+1,kme,jts-1:jte+1) * w(:,kme,:)
                                        
    end subroutine upwind_flux3
end module adv_fluxcorr