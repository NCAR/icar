! ------------------------------------------------------------------------------
!  copied from WRF/var/da/da_physics/da_sfc_wtq.inc to support the YSU pbl scheme. 
!  Specifically the calculation of psim, psih, which are made output vars here. 
!
!  Bert Kruyt 2022
!
! ------------------------------------------------------------------------------
module mod_pbl_utilities
    use icar_constants,           only : pi, gravity, Rd, Rw, cp, LH_vaporization,SVPT0
    ! use data_structures
    ! use domain_interface,   only : domain_t
    ! use options_interface,  only : options_t

    implicit none

contains

subroutine da_tp_to_qs( t, p, es, qs)

    !---------------------------------------------------------------------------
    ! Purpose: Convert T/p to saturation specific humidity.
    !
    !  Method: qs = es_alpha * es / ( p - ( 1 - rd_over_rv ) * es ).
    !          use Rogers & Yau (1989) formula: es = a exp( bTc / (T_c + c) )
    !--------------------------------------------------------------------------
 
    implicit none
 
    real, intent(in)  :: t, p
    real, intent(out) :: es, qs
    
    real              :: t_c              ! T in degreesC.
    real, parameter   :: es_alpha = 611.2 ! (= SVP1*1000)
    real, parameter   :: es_beta = 17.67  ! (= SVP2 = 17.67 )
    real, parameter   :: es_gamma = 243.5
    real, parameter    :: rd_over_rv = Rd / Rw! gas_constant / gas_constant_v
    real, parameter    :: rd_over_rv1 = 1.0 - rd_over_rv
    real, parameter    :: t_kelvin =SVPT0
    ! if (trace_use_dull) call da_trace_entry("da_tp_to_qs")
 
    !---------------------------------------------------------------------------
    ! [1.0] initialise:
    !---------------------------------------------------------------------------
    ! write(*,*) "t_kelvin",t_kelvin
    ! write(*,*) "t ",t
    t_c = t - t_kelvin
    
    !---------------------------------------------------------------------------
    ! [2.0] Calculate saturation vapour pressure:
    !---------------------------------------------------------------------------
 
    es = es_alpha * exp( es_beta * t_c / ( t_c + es_gamma ) )
     
    !---------------------------------------------------------------------------
    ! [3.0] Calculate saturation specific humidity:
    !---------------------------------------------------------------------------
 
    qs = rd_over_rv * es / ( p - rd_over_rv1 * es )
 
    ! if (trace_use_dull) call da_trace_exit("da_tp_to_qs")
 
 end subroutine da_tp_to_qs
 





subroutine da_sfc_wtq (psfc, tg, ps, ts, qs, us, vs, &
    hs, roughness, xland, dx, u10, v10, t2, q2, regime, psim, psih,  &
    has_lsm, regime_wrf, qsfc_wrf, znt_wrf, ust_wrf, mol_wrf, hfx, qfx, pblh, ims, ime, jms, jme)
 
    !---------------------------------------------------------------------------
    ! Purpose: Calculate the  10m wind, 2m temperature and moisture based on the
    ! similarity theory/
    !
    !  The unit for pressure   : psfc, ps          is Pa.
    !  The unit for temperature: tg, ts, t2        is K.
    !  The unit for moisture   : qs, q2            is kg/kg.
    !  The unit for wind       : us, vs, u10, v10  is m/s.
    !  The unit for height     : hs, roughness     is m.
    !  xland and regime are dimensionless. (BK: unitless?)
    !
    ! History: Nov 2010 - improve calculation consistency with WRF model (Eric Chiang)
    !          Jul 2015 - further improvement on consistency
    !
    ! Reference:
    ! ---------
    !
    !  input Variables:
    ! 
    !   psfc, tg               : surface pressure and ground temperature
    !   ps, ts, qs, us, vs, hs : model variable at lowlest half sigma level
    !   dx  (m)                : horizontal resolution
    !
    !
    !  Constants:
    !
    !   hs                     : height at the lowest half sigma level
    !   roughness              : roughness
    !   xland                  : land-water-mask (=2 water, =1 land)
    !
    !  output Variables:
    !
    !   regime                 : PBL regime
    !   u10, v10               : 10-m high observed wind components
    !   t2 , q2                : 2-m high observed temperature and mixing ratio
    !
    !---------------------------------------------------------------------------
    !  (BK 2022 made psim and psih outputs.)
    !                      psim  : mechanical psi at lowlest sigma level
    !                      psim2 : mechanical psi at 2m 
    !                      psimz : mechanical psi at 10m 
    !
    !---------------------------------------------------------------------------
 
    implicit none
 
    real, dimension( ims:ime, jms:jme ),    intent (in)  :: ps , ts , qs , us, vs
    real, dimension( ims:ime, jms:jme ),    intent (in)  :: psfc, tg
    real, dimension( ims:ime, jms:jme ),    intent (in)  :: hs, roughness , xland
    ! integer, dimension(:,:), intent(in) :: xland
    real, dimension( ims:ime, jms:jme ),    intent (out) :: regime
    real, dimension( ims:ime, jms:jme ),    intent (out) :: psim, psih
    real, dimension( ims:ime, jms:jme ),    intent (out) :: u10, v10, t2, q2
    logical,                        intent(in), optional :: has_lsm
    real,                           intent(in), optional :: regime_wrf, qsfc_wrf, znt_wrf,  mol_wrf
    real,                           intent(in), optional ::  pblh
    real,    dimension( ims:ime, jms:jme ), intent(in), optional :: hfx, qfx, ust_wrf
    integer,                                 intent(in) :: ims, ime, jms, jme
 
    ! logical :: use_table = .true.
    logical :: use_ust_wrf = .false.
    logical :: vconv_wrf
    integer :: nn, nz, n2, i, j
    real    :: rr, rz, r2
    real    :: cqs2, chs2, rho, rhox, fluxc, visc, restar, z0t, z0q
 
    ! h10 is the height of 10m where the wind observed
    ! h2  is the height of 2m where the temperature and 
    !        moisture observed.
 
    real, parameter :: h10 = 10.0, h2 = 2.0
    
    ! Default roughness over the land
 
    real, parameter :: zint0 = 0.01 
    
    ! Von Karman constant
 
    real, parameter :: k_kar = 0.4
    
    ! Working variables
 
    real :: Vc2, Va2, V2, vc, wspd
    real :: rib, rcp, xx, yy, cc
    real :: psiw, psiz, mol, ust, hol, holz, hol2
    real :: psimz, psim2, psihz, psih2  !psim, psih, ! BK 2022 now output vars for YSU
    real :: psit, psit2, psiq, psiq2
    real :: gzsoz0, gz10oz0, gz2oz0
    real :: eg, qg, tvg, tvs, tvs2
    real :: ths, thg, thvs, thvg, thvs2, vsgd, vsgd2, dx
    real :: zq0, z0, gas_constant

    real, parameter :: ka = 2.4E-5

    ! if (trace_use_dull) call da_trace_entry("da_sfc_wtq")
    gas_constant=Rd
    rcp = gas_constant/cp


    do j=jms,jme ! BK added these loops to make 2d
        do i=ims,ime

            ! 1 Compute the roughness length based upon season and land use 
        
            ! 1.1 Define the roughness length
        
            z0 = roughness(i,j)
        
            if (z0 < 0.0001) z0 = 0.0001
        
            if ( present(znt_wrf) ) then
            if ( znt_wrf > 0.0 ) then
                z0 = znt_wrf
            end if
            end if
        
            ! 1.2 Define the rouhgness length for moisture
        
            if (xland(i,j) .ge. 1.5) then
            zq0 = z0
            else
            zq0 =  zint0
            end if
        
            ! 1.3 Define the some constant variable for psi
        
            gzsoz0 = log(hs(i,j)/z0)
        
            gz10oz0 = log(h10/z0)
        
            gz2oz0 = log(h2/z0)
        
        
            ! 2. Calculate the virtual temperature
        
            ! 2.1 Compute Virtual temperature on the lowest half sigma level
        
            tvs  = ts(i,j) * (1.0 + 0.608 * qs(i,j))
        
            ! 2.2 Convert ground virtual temperature assuming it's saturated
            ! write(*,*)"tg(i,j)", tg(i,j)

            call da_tp_to_qs(tg(i,j), psfc(i,j), eg, qg) !output qg is specific humidity
            qg = qg*(1.0-qg) !hcl convert to mixing ratio
            if ( present(qsfc_wrf) ) then
            if ( qsfc_wrf > 0.0 ) then
                qg = qsfc_wrf
            end if
            endif
        
            tvg  = tg(i,j) * (1.0 + 0.608 * qg)
        
            ! 3.  Compute the potential temperature
        
            ! 3.1 Potential temperature on the lowest half sigma level
        
            ths  = ts(i,j) * (1000.0 / (ps(i,j)/100.0)) ** rcp
        
            ! 3.2 Potential temperature at the ground
        
            thg  = tg(i,j) * (1000.0 / (psfc(i,j)/100.0)) ** rcp
        
            ! 4. Virtual potential temperature
        
            ! 4.1 Virtual potential temperature on the lowest half sigma level
        
            thvs = tvs * (1000.0 / (ps(i,j)/100.0)) ** rcp
        
            ! 4.2 Virtual potential temperature at ground
        
            thvg = tvg * (1000.0 / (psfc(i,j)/100.0)) ** rcp
        
        
            ! 5.  BULK RICHARDSON NUMBER AND MONI-OBUKOV LENGTH
        
            ! 5.1 Velocity
            
            !     Wind speed:
        
            Va2 =   us(i,j)*us(i,j) + vs(i,j)*vs(i,j)
            !  
            !     Convective velocity:
        
            vconv_wrf = .false.
            if ( present(hfx) .and. present(qfx) .and. present(pblh) ) then
            ! calculating vconv over land following wrf method
            if ( pblh > 0.0 ) then
                vconv_wrf = .true.
            end if
            end if
        
            if (thvg >= thvs) then
            ! prior to V3.7, Vc2 = 4.0 * (thvg - thvs)
            Vc2 = thvg - thvs
            else
            Vc2 = 0.0
            end if
            if ( xland(i,j) < 1.5 ) then !land
            if ( vconv_wrf ) then
                ! following the calculation as in module_sf_sfclay.F
                rhox = psfc(i,j)/(gas_constant*tvg)
                fluxc = max(hfx(i,j)/rhox/cp+0.608*tvg*qfx(i,j)/rhox, 0.0)
                vc = (gravity/tg(i,j)*pblh*fluxc)**0.33
                vc2 = vc*vc
            end if
            end if
        
            ! Calculate Mahrt and Sun low-res correction         ! Add by Eric Chiang ( July 2010 )
            vsgd = 0.32 * (max(dx/5000.-1.,0.))**0.33            ! Add by Eric Chiang ( July 2010 )
            vsgd2 = vsgd * vsgd                                  ! Add by Eric Chiang ( July 2010 )
            
            V2  = Va2 + Vc2 + vsgd2                              ! Add by Eric Chiang ( July 2010 )  
            wspd = sqrt(v2)
            wspd = max(wspd,0.1)
            v2 = wspd*wspd
        
            ! 5.2 Bulk richardson number
        
            rib = (gravity * hs(i,j) / ths) * (thvs - thvg) / V2
        
            ! if previously unstable, do not let into regime 1 and 2
            if ( present(mol_wrf) ) then
            if ( mol_wrf < 0.0 ) rib = min(rib, 0.0)
            end if
        
            !  Calculate   ust, m/L (mol), h/L (hol)
        
            psim(i,j) = 0.0
            psih(i,j) = 0.0
        
            ! Friction speed
        
            if ( present(ust_wrf) ) then
            if ( ust_wrf(i,j) > 0.0 ) then
                use_ust_wrf = .true.
                ust = ust_wrf(i,j)
            end if
            end if
            if ( .not. use_ust_wrf ) then
            !ust = 0.0001  !init value as in phys/module_physics_init.F
            ust = k_kar * sqrt(v2) /(gzsoz0 - psim(i,j))
            end if
        
            ! Heat flux factor
        
            if ( present(mol_wrf) ) then
            mol = mol_wrf
            else
            mol = k_kar * (ths - thg)/(gzsoz0 - psih(i,j))
            !mol = 0.0
            end if
        
            ! set regimes based on rib
            if (rib .GE. 0.2) then
            ! Stable conditions (REGIME 1)
            regime(i,j) = 1.1
            else if ((rib .LT. 0.2) .AND. (rib .GT. 0.0)) then
            ! Mechanically driven turbulence (REGIME 2)
            regime(i,j) = 2.1
            else if (rib .EQ. 0.0) then
            ! Unstable Forced convection (REGIME 3)
            regime(i,j) = 3.1
            else 
            ! Free convection (REGIME 4)
            regime(i,j) = 4.1
            end if
        
            if ( present(regime_wrf) ) then
            if ( regime_wrf > 0.0 ) then
                regime(i,j) = regime_wrf
            end if
            end if
        
            ! 6.  CALCULATE PSI BASED UPON REGIME
        
            !if (rib .GE. 0.2) then
            if ( nint(regime(i,j)) == 1 ) then
            ! 6.1 Stable conditions (regime(i,j) 1)
            !     ---------------------------
            regime(i,j) = 1.1
            psim(i,j) = -10.0*gzsoz0
            psim(i,j) = max(psim(i,j),-10.0)
            psimz = h10/hs(i,j) * psim(i,j)
            psimz = max(psimz,-10.0)
            psim2 = h2 /hs(i,j) * psim(i,j)
            psim2 = max(psim2,-10.0)
            psih(i,j) = psim(i,j)
            psihz = psimz
            psih2 = psim2
        
            !else if ((rib .LT. 0.2) .AND. (rib .GT. 0.0)) then
            else if ( nint(regime(i,j)) == 2 ) then
        
            ! 6.2 Mechanically driven turbulence (regime(i,j) 2)
        
            regime(i,j) = 2.1
            psim(i,j) = (-5.0 * rib) * gzsoz0 / (1.1 - 5.0*rib)
            psim(i,j) = max(psim(i,j),-10.0)
            psimz = h10/hs(i,j) * psim(i,j)
            psimz = max(psimz,-10.0)
            psim2 = h2 /hs(i,j) * psim(i,j)
            psim2 = max(psim2,-10.0)
            psih(i,j) = psim(i,j)
            psihz = psimz
            psih2 = psim2
        
            !else if (rib .EQ. 0.0) then
            else if ( nint(regime(i,j)) == 3 ) then
            ! 6.3 Unstable Forced convection (regime(i,j) 3)
        
            regime(i,j) = 3.1
            psim(i,j) = 0.0
            psimz = 0.0
            psim2 = 0.0
            psih(i,j) = psim(i,j)
            psihz = psimz
            psih2 = psim2
        
            else 
            ! 6.4 Free convection (regime(i,j) 4)
            regime(i,j) = 4.1
        
            cc = 2.0 * atan(1.0)
        
            ! Ratio of PBL height to Monin-Obukhov length
        
            if (ust .LT. 0.01) then
                hol = rib * gzsoz0
            else
                hol = k_kar * gravity * hs(i,j) * mol / (ths * ust * ust)
            end if
        
            ! 6.4.2  Calculate n, nz, R, Rz
        
            holz = (h10 / hs(i,j)) * hol
            hol2 = (h2 / hs(i,j)) * hol
        
            hol = min(hol,0.0)
            hol = max(hol,-9.9999)
        
            holz = min(holz,0.0)
            holz = max(holz,-9.9999)
        
            hol2 = min(hol2,0.0)
            hol2 = max(hol2,-9.9999)
        
            ! 6.4.3 Calculate Psim & psih(i,j)
        
            !    if ( use_table ) then
            !       ! Using the look-up table:
            !       nn = int(-100.0 * hol)
            !       rr = (-100.0 * hol) - nn
            !       psim = psimtb(nn) + rr * (psimtb(nn+1) - psimtb(nn))
            !       psih(i,j) = psihtb(nn) + rr * (psihtb(nn+1) - psihtb(nn))
            !    else
                ! Using the continuous function:
                xx = (1.0 - 16.0 * hol) ** 0.25
                yy = log((1.0+xx*xx)/2.0)
                psim(i,j) = 2.0 * log((1.0+xx)/2.0) + yy - 2.0 * atan(xx) + cc
                psih(i,j) = 2.0 * yy
            !    end if
        
            !    if ( use_table ) then
            !       ! Using the look-up table:
            !       nz = int(-100.0 * holz)
            !       rz = (-100.0 * holz) - nz
            !       psimz = psimtb(nz) + rz * (psimtb(nz+1) - psimtb(nz))
            !       psihz = psihtb(nz) + rz * (psihtb(nz+1) - psihtb(nz))
            !    else
                ! Using the continuous function:
                xx = (1.0 - 16.0 * holz) ** 0.25
                yy = log((1.0+xx*xx)/2.0)
                psimz = 2.0 * log((1.0+xx)/2.0) + yy - 2.0 * atan(xx) + cc
                psihz = 2.0 * yy
            !    end if
        
            !    if ( use_table ) then
            !       ! Using the look-up table:
            !       n2 = int(-100.0 * hol2)
            !       r2 = (-100.0 * hol2) - n2
            !       psim2 = psimtb(n2) + r2 * (psimtb(n2+1) - psimtb(n2))
            !       psih2 = psihtb(n2) + r2 * (psihtb(n2+1) - psihtb(n2))
            !    else
                ! Using the continuous function:
                xx = (1.0 - 16.0 * hol2) ** 0.25
                yy = log((1.0+xx*xx)/2.0)
                psim2 = 2.0 * log((1.0+xx)/2.0) + yy - 2.0 * atan(xx) + cc
                psih2 = 2.0 * yy
            !    end if
        
            ! 6.4.4 Define the limit value for psim & psih(i,j)
        
            psim(i,j) = min(psim(i,j),0.9*gzsoz0)
            psimz = min(psimz,0.9*gz10oz0)
            psim2 = min(psim2,0.9*gz2oz0)
            psih(i,j) = min(psih(i,j),0.9*gzsoz0)
            psihz = min(psihz,0.9*gz10oz0)
            psih2 = min(psih2,0.9*gz2oz0)
            end if  ! regime(i,j)
        
            ! 7.  Calculate psi for wind, temperature and moisture
        
            psiw = gzsoz0 - psim(i,j)
            psiz = gz10oz0 - psimz
            psit = max(gzsoz0-psih(i,j), 2.0)
            psit2 = gz2oz0 - psih2
        
            if ( .not. use_ust_wrf ) then
            ! re-calculate ust since psim(i,j) is now available
            ust = k_kar * sqrt(v2) /(gzsoz0 - psim(i,j))
            end if
        
            psiq  = log(k_kar*ust*hs(i,j)/ka + hs(i,j) / zq0) - psih(i,j)
            psiq2 = log(k_kar*ust*h2/ka + h2 / zq0) - psih2
        
            !V3.7, as in module_sf_sfclay.F
            if ( xland(i,j) >= 1.5 ) then !water
            visc = (1.32+0.009*(ts(i,j)-273.15))*1.e-5
            restar = ust*z0/visc
            z0t = (5.5e-5)*(restar**(-0.60))
            z0t = min(z0t,1.0e-4)
            z0t = max(z0t,2.0e-9)
            z0q = z0t
            psiq  = max(log((hs(i,j)+z0q)/z0q)-psih(i,j),  2.)
            psit  = max(log((hs(i,j)+z0t)/z0t)-psih(i,j),  2.)
            psiq2 = max(log((2.+z0q)/z0q)-psih2, 2.)
            psit2 = max(log((2.+z0t)/z0t)-psih2, 2.)
            end if
        
            ! 8.  Calculate 10m wind, 2m temperature and moisture
        
            u10(i,j) = us(i,j) * psiz / psiw
            v10(i,j) = vs(i,j) * psiz / psiw
            t2(i,j) = (thg + (ths - thg)*psit2/psit)*((psfc(i,j)/100.0)/1000.0)**rcp
            q2(i,j) = qg + (qs(i,j) - qg)*psiq2/psiq 
        
            if ( present(has_lsm) ) then
            if ( has_lsm ) then
                !cqs2: 2m surface exchange coefficient for moisture
                !chs2: 2m surface exchange coefficient for heat
                cqs2 = ust*k_kar/psiq2
                if (xland(i,j) .ge. 1.5) then
                    !water
                    chs2 = ust*k_kar/psit2
                else
                    !land
                    chs2 = cqs2 !as in subroutine lsm in phys/module_sf_noahdrv.F
                end if

                !re-calculate T2/Q2 as in module_sf_sfcdiags.F
                rho  = psfc(i,j)/(gas_constant*tg(i,j))
                if ( cqs2 < 1.e-5 ) then
                    q2(i,j) = qg
                else
                    if ( present(qfx) ) then
                        q2(i,j) = qg - qfx(i,j)/(rho*cqs2)
                    end if
                end if
                if ( chs2 < 1.e-5 ) then
                    t2(i,j) = tg(i,j)
                else
                    if ( present(hfx) ) then
                        t2(i,j) = tg(i,j) - hfx(i,j)/(rho*cp*chs2)
                    end if
                end if
            end if
            end if

    ! if (trace_use_dull) call da_trace_exit("da_sfc_wtq")
        enddo
    enddo
    end subroutine da_sfc_wtq


    ! !--------------------------------------------------------
    ! !
    ! ! from old ICAR model (v1) - currently notused
    ! !
    ! !--------------------------------------------------------

    ! subroutine calc_surface_stuff
    !             ! ----- start surface layer calculations usually done by surface layer scheme ----- !
    !     write(*,*) "start surface layer calculations"
    !     if (options%physics%boundarylayer==kPBL_SIMPLE) then
    !         write(*,*) "calculate surface layer based on log wind profile"
    !         ! ! temporary constant
    !         ! ! use log-law of the wall to convert from first model level to surface
    !         ! currw = karman / log((domain%z(2:nx-1,1,2:ny-1)-domain%terrain(2:nx-1,2:ny-1)) &
    !         !         / domain%znt(2:nx-1,2:ny-1))
    !         ! ! use log-law of the wall to convert from surface to 10m height
    !         ! lastw = log(10.0 / domain%znt(2:nx-1,2:ny-1)) / karman
    !         ! domain%ustar(2:nx-1,2:ny-1) = domain%Um(2:nx-1,1,2:ny-1) * currw
    !         ! domain%u10(2:nx-1,2:ny-1) = domain%ustar(2:nx-1,2:ny-1) * lastw
    !         ! domain%ustar(2:nx-1,2:ny-1) = domain%Vm(2:nx-1,1,2:ny-1) * currw
    !         ! domain%v10(2:nx-1,2:ny-1) = domain%ustar(2:nx-1,2:ny-1) * lastw
        
    !         ! ! now calculate master ustar based on U and V combined in quadrature
    !         ! domain%ustar(2:nx-1,2:ny-1) = sqrt(domain%Um(2:nx-1,1,2:ny-1)**2 & 
    !         !                             + domain%Vm(2:nx-1,1,2:ny-1)**2) * currw
    !         ! ! counter is just a variable helping me to detect how much rounds
    !         ! ! this subroutine went through
    !         ! write(*,*) "Counter: ", counter
    !         ! counter = counter + 1
    !         ! write(*,*) "Counter: ", counter
    !     elseif (options%physics%boundarylayer==kPBL_YSU) then
    !         ! start surface layer calculations introduced by Patrik Bohlinger
    !         write(*,*) "calculate surface layer based on monin-obukhov similarity theory"

    !         ! ----- start temporary solution ----- !
    !         ! use log-law of the wall to convert from first model level to surface
    !         currw = karman / log((domain%z(2:nx-1,1,2:ny-1)-domain%terrain(2:nx-1,2:ny-1)) &
    !                / domain%znt(2:nx-1,2:ny-1))
    !         ! use log-law of the wall to convert from surface to 10m height
    !         lastw = log(10.0 / domain%znt(2:nx-1,2:ny-1)) / karman
    !         ! calculate ustar = horizontal wind speed scale
    !         if (counter==1) then
    !             domain%ustar_new(2:nx-1,2:ny-1) = domain%ustar(2:nx-1,2:ny-1)
    !         endif
    !         ! preventing ustar from being smaller than 0.1 as it could be under
    !         ! very stable conditions, Jiminez et al. 2012
    !         where(domain%ustar_new(2:nx-1,2:ny-1) < 0.1)
    !             domain%ustar_new(2:nx-1,2:ny-1) = 0.1
    !         endwhere
    !         !domain%ustar(2:nx-1,2:ny-1) = domain%Um(2:nx-1,1,2:ny-1) * currw
    !         domain%u10(2:nx-1,2:ny-1) = domain%ustar_new(2:nx-1,2:ny-1) * lastw
    !         domain%ustar(2:nx-1,2:ny-1) = domain%Vm(2:nx-1,1,2:ny-1) * currw
    !         domain%v10(2:nx-1,2:ny-1) = domain%ustar_new(2:nx-1,2:ny-1) * lastw
    !         !now calculate master ustar based on U and V combined in quadrature
    !         domain%wspd3d(2:nx-1,1:nz,2:ny-1) = sqrt(domain%Um(2:nx-1,1:nz,2:ny-1)**2 & 
    !                                             + domain%Vm(2:nx-1,1:nz,2:ny-1)**2) 
    !         ! added by Patrik Bohlinger in case we need this later for YSU (some variables seem to 
    !         ! be 3D in articles like the YSU paper Hong et al. 2006)
    !         domain%wspd(2:nx-1,2:ny-1) = sqrt(domain%Um(2:nx-1,1,2:ny-1)**2 & 
    !                                    + domain%Vm(2:nx-1,1,2:ny-1)**2) 
    !         ! added by Patrik Bohlinger since we need this as input for YSU
    !          domain%ustar(2:nx-1,2:ny-1) = domain%wspd(2:nx-1,2:ny-1) * currw
    !         ! ----- end temporary solution ----- !

    !         ! compute z above ground used for estimating indices for 
    !         domain%z_agl(2:nx-1,2:ny-1) = (domain%z(2:nx-1,1,2:ny-1)-domain%terrain(2:nx-1,2:ny-1)) 
    !         !added by Patrik in case we need this later for YSU
    !         ! calculate the Bulk-Richardson number Rib
    !         domain%thv(2:nx-1,2:ny-1) = domain%th(2:nx-1,1,2:ny-1) & 
    !                                     *(1+0.608*domain%qv(2:nx-1,1,2:ny-1)*1000)    
    !         ! should domain%qv be multiplied by 1000? Did it since domain%qv is in kg/kg and not in g/kg 
    !         ! normally should be specific humidity and not mixing ratio domain%qv but for first order approach does not matter
    !         domain%thv3d(2:nx-1,1:nz,2:ny-1) = domain%th(2:nx-1,1:nz,2:ny-1) & 
    !                                            * (1+0.608*domain%qv(2:nx-1,1:nz,2:ny-1)*1000)   !thv 3D
    !         domain%thvg(2:nx-1,2:ny-1) = (domain%t2m(2:nx-1,2:ny-1)/domain%pii(2:nx-1,1,2:ny-1)) &
    !                                      *(1+0.608*domain%qv(2:nx-1,1,2:ny-1)*1000)
    !         ! t2m should rather be used than skin_t
    !         domain%thg(2:nx-1,2:ny-1) = domain%t2m(2:nx-1,2:ny-1)/domain%pii(2:nx-1,1,2:ny-1) 
    !         ! t2m should rather be used than skin_t

    !         ! variables described for YSU but probably not needed to be calculated outside of the scheme:
    !         !domain%wstar(2:nx-1,2:ny-1) = domain%ustar(2:nx-1,2:ny-1) / domain%psim(2:nx-1,2:ny-1) ! wstar = vertical wind speed scale
    !         !domain%thT(2:nx-1,2:ny-1) = propfact * (virtual heat flux)/ domain%wstar ! virtual temperature excess
    !         !domain%thvg(2:nx-1,2:ny-1) = domain%thv(2:nx-1,2:ny-1) ! for init thvg=thv since thT = 0, 
    !         !t2m should rather be used than skin_t, b=proportionality factor=7.8, Hong et al, 2006

    !         ! find value of pbl heights for wspd3d
    !         !domain%PBLh(2:nx-1,2:ny-1) = Rib_cr * domain%thv(2:nx-1,2:ny-1) * domain%wspd(2:nx-1,2:ny-1)**2 & 
    !         !                             / gravity * (domain%thv(2:nx-1,2:ny-1) - domain%thvg(2:nx-1,2:ny-1)) !U^2 and thv are from height PBLh in equation

    !         !write(*,*) "max min domain%pbl_height: ", maxval(domain%pbl_height), minval(domain%pbl_height)
    !         !write(*,*) "max min domain%PBLh: ", maxval(domain%PBLh), minval(domain%PBLh) ! introduced the 
    !         !PBLh variabel to not overwrite pbl_height and compare new with old calculations as the pbl 
    !         !height is one of the most crucial factors of the non-local surface layer calculations needed by the YSU-scheme

    !         ! Constraint to prevent Rib from becoming too high a lower limit of 0.1 is
    !         ! applied in for the original surface layer formulation Jiminez et al 2012
    !         where(domain%wspd(2:nx-1,2:ny-1) < 0.1)
    !             domain%wspd(2:nx-1,2:ny-1) = 0.1
    !         endwhere

    !         domain%Rib(2:nx-1,2:ny-1) = gravity/domain%th(2:nx-1,1,2:ny-1) * domain%z_agl(2:nx-1,2:ny-1) &
    !                                     * (domain%thv(2:nx-1,2:ny-1) - domain%thvg(2:nx-1,2:ny-1)) &
    !                                     / domain%wspd(2:nx-1,2:ny-1)**2
    !         ! From Jiminez et al. 2012, from what height should the theta variables really be, Rib is a function of height? To my understanding the appropriate height is the lower most level.

    !         ! calculate the integrated similarity functions
    !         where(domain%Rib(2:nx-1,2:ny-1) >= 0.2)
    !             !regime = 1, very stable night time conditions
    !             domain%psim(2:nx-1,2:ny-1) = -10*log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1))
    !             domain%psim10(2:nx-1,2:ny-1) = -10*log(10/domain%znt(2:nx-1,2:ny-1))
    !             domain%psim2m(2:nx-1,2:ny-1) = -10*log(2/domain%znt(2:nx-1,2:ny-1))
    !             domain%psih(2:nx-1,2:ny-1) = domain%psim(2:nx-1,2:ny-1)
    !             domain%psih2m(2:nx-1,2:ny-1) = domain%psim2m(2:nx-1,2:ny-1)
    !             !impose constraints
    !             where (domain%Rib(2:nx-1,2:ny-1) >= 0.2 .and. domain%psim(2:nx-1,2:ny-1) < -10.)
    !                 domain%psim(2:nx-1,2:ny-1) = -10.
    !             endwhere
    !             where (domain%Rib(2:nx-1,2:ny-1) >= 0.2 .and. domain%psih(2:nx-1,2:ny-1) < -10.)
    !                 domain%psih(2:nx-1,2:ny-1) = -10.
    !             endwhere
    !             where (domain%Rib(2:nx-1,2:ny-1) >= 0.2 .and. domain%psim10(2:nx-1,2:ny-1) < -10.)
    !                 domain%psim10(2:nx-1,2:ny-1) = -10.
    !             endwhere
    !             where (domain%Rib(2:nx-1,2:ny-1) >= 0.2 .and. domain%psih2m(2:nx-1,2:ny-1) < -10.)
    !                 domain%psih2m(2:nx-1,2:ny-1) = -10.
    !             endwhere
    !             where (domain%Rib(2:nx-1,2:ny-1) >= 0.2 .and. domain%psim2m(2:nx-1,2:ny-1) < -10.)
    !                 domain%psim2m(2:nx-1,2:ny-1) = -10.
    !             endwhere
    !         elsewhere (domain%Rib(2:nx-1,2:ny-1) < 0.2 .and. domain%Rib(2:nx-1,2:ny-1) > 0.0)
    !             !regime = 2, damped mechanical turbulence
    !             domain%psim(2:nx-1,2:ny-1) = -5*domain%Rib(2:nx-1,2:ny-1)*log(domain%z_agl(2:nx-1,2:ny-1) & 
    !                                          /domain%znt(2:nx-1,2:ny-1))/(1.1-5*domain%Rib(2:nx-1,2:ny-1))
    !             domain%psim10(2:nx-1,2:ny-1) = -5*domain%Rib(2:nx-1,2:ny-1)*log(10/domain%znt(2:nx-1,2:ny-1)) &
    !                                          /(1.1-5*domain%Rib(2:nx-1,2:ny-1)) ! Should maybe compute Rib at 10m as well?
    !             domain%psim2m(2:nx-1,2:ny-1) = -5*domain%Rib(2:nx-1,2:ny-1)*log(2/domain%znt(2:nx-1,2:ny-1)) &
    !                                          /(1.1-5*domain%Rib(2:nx-1,2:ny-1)) ! Should maybe compute Rib at 2m as well?
    !             domain%psih(2:nx-1,2:ny-1) = domain%psim(2:nx-1,2:ny-1)
    !             domain%psih2m(2:nx-1,2:ny-1) = domain%psim2m(2:nx-1,2:ny-1)
    !             !impose constraints
    !             where (domain%Rib(2:nx-1,2:ny-1) < 0.2 .and. &
    !                     domain%Rib(2:nx-1,2:ny-1) > 0.0 .and. &
    !                     domain%psim(2:nx-1,2:ny-1) < -10.)
    !                 domain%psim(2:nx-1,2:ny-1) = -10.
    !             endwhere
    !             where (domain%Rib(2:nx-1,2:ny-1) < 0.2 .and. &
    !                     domain%Rib(2:nx-1,2:ny-1) > 0.0 .and. &
    !                     domain%psih(2:nx-1,2:ny-1) < -10.)
    !                 domain%psih(2:nx-1,2:ny-1) = -10.
    !             endwhere
    !             where (domain%Rib(2:nx-1,2:ny-1) < 0.2 .and. &
    !                     domain%Rib(2:nx-1,2:ny-1) > 0.0 .and. &
    !                     domain%psim10(2:nx-1,2:ny-1) < -10.)
    !                 domain%psim10(2:nx-1,2:ny-1) = -10.
    !             endwhere
    !             where (domain%Rib(2:nx-1,2:ny-1) < 0.2 .and. &
    !                     domain%Rib(2:nx-1,2:ny-1) > 0.0 .and. &
    !                     domain%psih2m(2:nx-1,2:ny-1) < -10.)
    !                 domain%psih2m(2:nx-1,2:ny-1) = -10.
    !             endwhere
    !             where (domain%Rib(2:nx-1,2:ny-1) < 0.2 .and. &
    !                     domain%Rib(2:nx-1,2:ny-1) > 0.0 .and. &
    !                     domain%psim2m(2:nx-1,2:ny-1) < -10.)
    !                 domain%psim2m(2:nx-1,2:ny-1) = -10.
    !             endwhere
    !         elsewhere (domain%Rib(2:nx-1,2:ny-1).eq.0.0)
    !             !regime = 3, forced convection
    !             domain%psim(2:nx-1,2:ny-1) = 0.0
    !             domain%psim10(2:nx-1,2:ny-1) = 0.0
    !             domain%psih(2:nx-1,2:ny-1) = 0.0
    !             domain%psih2m(2:nx-1,2:ny-1) = 0.0
    !         elsewhere (domain%Rib(2:nx-1,2:ny-1) < 0)
    !             !regime = 4, free convection
    !             ! constraints
    !             where (domain%Rib(2:nx-1,2:ny-1) < 0. .and. domain%zol(2:nx-1,2:ny-1) < -9.9999)
    !                 domain%zol(2:nx-1,2:ny-1) = -9.9999
    !             endwhere
    !             where (domain%Rib(2:nx-1,2:ny-1) < 0. .and. domain%zol(2:nx-1,2:ny-1) > 0.)
    !                 domain%zol(2:nx-1,2:ny-1) = 0.
    !             endwhere
    !             where (domain%Rib(2:nx-1,2:ny-1) < 0. .and. domain%zol10(2:nx-1,2:ny-1) < -9.9999)
    !                 domain%zol10(2:nx-1,2:ny-1) = -9.9999
    !             endwhere
    !             where (domain%Rib(2:nx-1,2:ny-1) < 0. .and. domain%zol10(2:nx-1,2:ny-1) > 0.)
    !                 domain%zol10(2:nx-1,2:ny-1) = 0.
    !             endwhere
    !             where (domain%Rib(2:nx-1,2:ny-1) < 0. .and. domain%zol2m(2:nx-1,2:ny-1) < -9.9999)
    !                 domain%zol2m(2:nx-1,2:ny-1) = -9.9999              
    !             endwhere
    !             where (domain%Rib(2:nx-1,2:ny-1) < 0. .and. domain%zol2m(2:nx-1,2:ny-1) > 0.)
    !                 domain%zol2m(2:nx-1,2:ny-1) = 0.
    !             endwhere
    !             domain%psix(2:nx-1,2:ny-1) = (1.-16.*(domain%zol(2:nx-1,2:ny-1)))**0.25
    !             domain%psix10(2:nx-1,2:ny-1) = (1.-16.*(domain%zol10(2:nx-1,2:ny-1)))**0.25
    !             domain%psix2m(2:nx-1,2:ny-1) = (1.-16.*(domain%zol2m(2:nx-1,2:ny-1)))**0.25
    !             domain%psim(2:nx-1,2:ny-1) = 2.*log((1.+domain%psix(2:nx-1,2:ny-1))/2.) & 
    !                                          + log((1.+domain%psix(2:nx-1,2:ny-1)**2.)/2.) &
    !                                          - 2.*atan(domain%psix(2:nx-1,2:ny-1))+pi/2.
    !             domain%psim10(2:nx-1,2:ny-1) = 2.*log((1.+domain%psix10(2:nx-1,2:ny-1))/2.) &
    !                                          + log((1.+domain%psix10(2:nx-1,2:ny-1)**2.)/2.) &
    !                                          - 2.*atan(domain%psix10(2:nx-1,2:ny-1))+pi/2.
    !             domain%psih(2:nx-1,2:ny-1) = 2.*log((1.+domain%psix(2:nx-1,2:ny-1)**2.)/2.)
    !             domain%psih2m(2:nx-1,2:ny-1) = 2.*log((1.+domain%psix2m(2:nx-1,2:ny-1)**2.)/2.)
    !         endwhere
    !         ! constrain psim and psih also for 2m and 10m
    !         where (domain%psim(2:nx-1,2:ny-1) > 0.9 * log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1)))
    !             domain%psim(2:nx-1,2:ny-1) = 0.9 * log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1))
    !         endwhere
    !         where (domain%psim2m(2:nx-1,2:ny-1) > 0.9 * log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1)))
    !             domain%psim2m(2:nx-1,2:ny-1) = 0.9 * log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1))
    !         endwhere
    !         where (domain%psim10(2:nx-1,2:ny-1) > 0.9 * log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1)))
    !             domain%psim10(2:nx-1,2:ny-1) = 0.9 * log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1))
    !         endwhere
    !         where (domain%psih(2:nx-1,2:ny-1) > 0.9 * log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1)))
    !             domain%psih(2:nx-1,2:ny-1) = 0.9 * log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1))
    !         endwhere
    !         where (domain%psih2m(2:nx-1,2:ny-1) > 0.9 * log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1)))
    !             domain%psih2m(2:nx-1,2:ny-1) = 0.9 * log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1))
    !         endwhere

    !         ! calculate thstar = temperature scale
    !         domain%thstar(2:nx-1,2:ny-1) = karman*(domain%th(2:nx-1,1,2:ny-1)-domain%thg) &
    !                                        / log(domain%z_agl(2:nx-1,2:ny-1) & 
    !                                        / domain%znt(2:nx-1,2:ny-1))-domain%psih(2:nx-1,2:ny-1)
    !         ! Averaging ustar with ustar from previous time step to suppress
    !         ! large oscillations
    !         domain%ustar_tmp(2:nx-1,2:ny-1) = karman*domain%wspd(2:nx-1,2:ny-1)/(log(domain%z_agl(2:nx-1,2:ny-1) & 
    !                                           / domain%znt(2:nx-1,2:ny-1))-domain%psim(2:nx-1,2:ny-1))
    !         domain%ustar_new(2:nx-1,2:ny-1) = (domain%ustar_tmp(2:nx-1,2:ny-1) + domain%ustar_new(2:nx-1,2:ny-1))/2.0
    !         ! preventing ustar from being smaller than 0.1 as it could be under
    !         ! very stable conditions, Jiminez et al. 2012
    !         where(domain%ustar_new(2:nx-1,2:ny-1) < 0.1)
    !             domain%ustar_new(2:nx-1,2:ny-1) = 0.1
    !         endwhere
    !         ! calculate the Monin-Obukhov  stability parameter zol (z over l)
    !         ! using ustar from the similarity theory
    !         domain%zol(2:nx-1,2:ny-1) = (karman*gravity*domain%z_agl(2:nx-1,2:ny-1))/domain%th(2:nx-1,1,2:ny-1) &
    !                                     * domain%thstar(2:nx-1,2:ny-1)/(domain%ustar_new(2:nx-1,2:ny-1) &
    !                                     * domain%ustar_new(2:nx-1,2:ny-1))
    !         domain%zol10(2:nx-1,2:ny-1) = (karman*gravity*10)/domain%th(2:nx-1,1,2:ny-1) * domain%thstar(2:nx-1,2:ny-1) &
    !                                     / (domain%ustar_new(2:nx-1,2:ny-1)*domain%ustar_new(2:nx-1,2:ny-1))
    !         domain%zol2m(2:nx-1,2:ny-1) = (karman*gravity*2)/domain%th(2:nx-1,1,2:ny-1) * domain%thstar(2:nx-1,2:ny-1) &
    !                                     / (domain%ustar_new(2:nx-1,2:ny-1)*domain%ustar_new(2:nx-1,2:ny-1))
    !         ! calculate pblh over l using ustar and thstar from the similarity theory
    !         domain%hol(2:nx-1,2:ny-1) = (karman*gravity*domain%PBLh(2:nx-1,2:ny-1))/domain%th(2:nx-1,1,2:ny-1) &
    !                                     * domain%thstar(2:nx-1,2:ny-1)/(domain%ustar_new(2:nx-1,2:ny-1) &
    !                                     * domain%ustar_new(2:nx-1,2:ny-1))
    !         ! arbitrary variables
    !         domain%gz1oz0(2:nx-1,2:ny-1)=log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1))
    !         ! calculating dtmin
    !         domain%dtmin = domain%dt / 60.0
    !         ! p_top as a scalar, choosing just minimum from ptop as a start
    !         p_top = minval(domain%ptop)
    !         ! compute the dimensionless bulk coefficent for momentum, heat and
    !         ! moisture
    !         !domain%exch_m(2:nx-1,2:ny-1) = (karman**2)/(log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1)) &
    !         !- domain%psim(2:nx-1,2:ny-1))**2
    !         !domain%exch_h(2:nx-1,2:ny-1) = (karman**2)/((log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1))&
    !         !- domain%psim(2:nx-1,2:ny-1))*(log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1))-domain%psih(2:nx-1,2:ny-1)))
    !         !domain%exch_q(2:nx-1,2:ny-1) = (karman**2)/((log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1)) &
    !         !- domain%psim(2:nx-1,2:ny-1)) * (log(domain%rho(2:nx-1,1,2:ny-1)*cp*karman*domain%ustar_new(2:nx-1,2:ny-1) &
    !         !*domain%z_agl(2:nx-1,2:ny-1)/cs)-psih(2:nx-1,2:ny-1)))

    !         ! counter is just a variable helping me to detect how much rounds this subroutine went through
    !         write(*,*) "Counter: ", counter
    !         counter = counter + 1
    !         write(*,*) "Counter: ", counter
    !         ! end surface layer calculations introduced by Patrik Bohlinger
    !     endif        
    !     write(*,*) "end surface layer calculations"
    !     ! ----- end sfc layer calculations ----- !
    ! end subroutine


end module mod_pbl_utilities