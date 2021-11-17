!>------------------------------------------------
!! Contains type definitions for a variety of model data strucutres
!! Also defines model constants (e.g. gravity, and MAXFILELENGTH)
!!
!! <pre>
!! General Field Definitions
!!
!! ---- 3D fields ---- NX x NZ x NY
!! u        = wind in east direction                        [m/s]
!! v        = wind in north direction                       [m/s]
!! w        = wind in vertical direction                    [m/s] (possibly scaled by dx/dz)
!!
!! p        = pressure                                      [pa]
!! th       = potential temperature                         [K]
!!
!! qv       = water vapor (mixing ratio)                    [kg/kg]
!! cloud    = cloud water                                   [kg/kg]
!! ice      = cloud ice                                     [kg/kg]
!! qrain    = rain mixing ratio                             [kg/kg]
!! qsnow    = snow mixing ratio                             [kg/kg]
!! qgrau    = graupel mixing ratio                          [kg/kg]
!! nice     = ice number concentration                      [1/cm^3]
!! nrain    = rain number concentration                     [1/cm^3]
!!
!! rho      = dry air density                               [kg/m^3]
!! pii      = exner function                                []
!!
!! ---- 2D fields ---- NX x NY
!!      ---- moisture fluxes ----
!! rain     = rain+crain+snow+graupel                       [mm]
!! crain    = convective rain at surface                    [mm]
!! snow     = snow at surface                               [mm]
!! graupel  = graupel at surface                            [mm]
!!
!!      ---- energy fluxes ----
!! sensible_heat = Sensible heat flux from surface          [W/m^2]
!! latent_heat   = Latent heat flux from surface            [W/m^2]
!! pbl_height    = Height of the planetary boundary layer   [m]
!!
!!      ---- Radiation variables ----
!! cloudfrac    = Cloud fraction                            [0-1]
!! swdown       = Shortwave down at land surface            [W/m^2]
!! lwdown       = Longwave down at land surface             [W/m^2]
!! lwup         = Lonwave up from the land surface          [W/m^2]
!!
!! ---- Land Surface variables ----
!!   3D fields ---- NX x NZ x NY (NZ = number of soil layers)
!! soil_t       = 3D Soil temperature                       [K]
!! soil_vwc     = 3D Soil volumetric water content          [m^3/m^3]
!!
!!   3D fields ---- NX x NY x N_Times (typically 1 or 12)
!! vegfrac      = vegetation cover fraction                 [%]
!!
!!   2D fields ---- NX x NY
!! skin_t       = Land surface skin temperature             [K]
!! soil_tdeep   = Temperature at the soil column bottom     [K]
!! snow_swe     = Snow water equivalent on the land surface [mm]
!! soil_totalmoisture = Soil column total water content     [mm]
!! soil_type    = Soil type (index for SOILPARM.TBL)        [1-n]
!! veg_type     = Vegetation type (index for VEGPARM.TBL)   [1-n]
!! landmask     = Map of Land vs Water grid cells           [0,1,2]
!!
!! ---- NOTE ----
!! dX_dt variables are the increment in boundary conditions between internal model time steps
!! some of these are 2d, some are 3d
!!
!! ---- model structure ----
!! terrain  = surface elevation                             [m]
!! z        = model layer height (at mid point)             [m]
!! dz       = Model layer thickness (between mass levels)   [m]
!! dz_inter = Layer thickness (between interface levels)    [m]
!!
!! sintheta = sine of the angle between grid and geographic coords   [-]
!! costheta = cosine of the angle between grid and geographic coords [-]
!! fzs      = buffered FFT(terrain) for linear wind calculations
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!>------------------------------------------------
module data_structures
    use, intrinsic :: iso_c_binding ! needed for fftw compatible complex types
    use time_object,        only : Time_type
    use time_delta_object,  only : time_delta_t
    use icar_constants           ! Many constants including things like fixed string lengths
    implicit none


! ------------------------------------------------
!   various data structures for use in geographic interpolation routines
! ------------------------------------------------
    ! contains the location of a specific grid point
    type position
        integer::x,y
    end type position
    ! contains location of surrounding 4 grid cells
    type fourpos
        integer::x(4),y(4)
    end type fourpos

    ! a geographic look up table for spatial interpolation, from x,y with weight w
    type geo_look_up_table
        ! x,y index positions, [n by m by 4] where there are 4 surrounding low-res points
        ! for every high resolution point grid point to interpolate to
        integer,allocatable, dimension(:,:,:)   :: x, y
        ! weights to use for each of the 4 surrounding gridpoints.  Sum(over axis 3) must be 1.0
        real,   allocatable, dimension(:,:,:)   :: w
    end type geo_look_up_table

    ! ------------------------------------------------
    ! A look up table for vertical interpolation. from z with weight w
    ! ------------------------------------------------
    type vert_look_up_table
        ! z index positions for all x,y,z points (x 2 for above and below z levels)
        integer,allocatable, dimension(:,:,:,:) :: z

        ! weights to use for each of the two surrounding points.  Sum (over axis 1) must be 1.0
        real,   allocatable, dimension(:,:,:,:) :: w
    end type vert_look_up_table

    ! ------------------------------------------------
    ! generic interpolable type so geo interpolation routines will work on winds, domain, or boundary conditions.
    ! ------------------------------------------------
    type interpolable_type
        ! all interpolables must have position (lat, lon, z)
        real, allocatable, dimension(:,:) :: lat,lon
        real, allocatable, dimension(:,:,:) :: z

        ! these are the look up tables that describe how to interpolate vertically (vert_lut) and horizontally (geolut)
        type(vert_look_up_table)::vert_lut
        type(geo_look_up_table)::geolut

        ! used to keep track of whether or not a particular error has been printed yet for this structure
        logical :: dx_errors_printed=.False.
        logical :: dy_errors_printed=.False.
    end type interpolable_type


    ! ------------------------------------------------
    ! type to contain external wind fields, only real addition is nfiles... maybe this could be folded in elsewhere?
    ! ------------------------------------------------
    type, extends(interpolable_type) :: wind_type
        !
        real, allocatable, dimension(:,:,:) :: u,v
        type(interpolable_type)             :: u_geo,v_geo
        real, allocatable, dimension(:,:)   :: terrain,dzdx,dzdy
        real :: dx
        integer :: nfiles
    end type wind_type

    ! ------------------------------------------------
    ! generic linearizable type so we can add linear wind field to domain or remove it from low-res (BC) U/V
    ! Because of the need to compute the Brunt-Vaisala frequency, this now includes many atmospheric fields
    ! ------------------------------------------------
    type, extends(interpolable_type) :: linearizable_type
        ! linear theory computes u,v at z.  Trying rho to mitigate boussinesq approx...
        real, allocatable, dimension(:,:,:) :: u,v,dz,rho
        ! these are needed to compute Brunt Vaisalla Frequency...
        real, allocatable, dimension(:,:,:) :: th                   ! potential temperature         [K]
        real, allocatable, dimension(:,:,:) :: p                    ! pressure                      [Pa]
        real, allocatable, dimension(:,:,:) :: pii                  ! exner function                [-]
        real, allocatable, dimension(:,:,:) :: qv                   ! water vapor mixing ratio      [kg/kg]
        ! used to determine if moist or dry brunt vaisala frequency is used
        real, allocatable, dimension(:,:,:) :: cloud,ice,qsnow,qrain ! hydrometeor mixing ratios    [kg/kg]
        real, allocatable, dimension(:,:,:) :: nsquared             ! Brunt-Vaisala frequency       [1/s^2]
        type(interpolable_type)             :: u_geo,v_geo          ! types to store the staggered grid lat/lon & geoluts
        real, allocatable, dimension(:,:)   :: terrain              ! ground surface                [m]
        real                                :: dx                   ! the horizontal grid spacing   [m]
        real, allocatable, dimension(:,:)   :: dzdx,dzdy            ! change in terrain / dx        [m/m]
        real, allocatable, dimension(:,:)   :: linear_mask          ! weights to multiply linear wind by (default=1)
        real, allocatable, dimension(:,:)   :: nsq_calibration      ! calibration parameter to multiply brunt-vaisala frequency by [0-]
        complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: fzs ! FFT(terrain)

        real, allocatable, dimension(:,:)   :: froude               ! store a distributed map of froude numbers

    end type linearizable_type


    ! ------------------------------------------------
    ! Data type to hold all of the array temporaries required by the lineary theory calculations
    ! e.g. k and l wave number arrays
    ! ------------------------------------------------
    type linear_theory_type
        real,                       allocatable, dimension(:,:) :: sig, k, l, kl
        complex(C_DOUBLE_COMPLEX),  allocatable, dimension(:,:) :: denom, msq, mimag, m, ineta

        complex(C_DOUBLE_COMPLEX),  pointer,     dimension(:,:) :: uhat, vhat
        complex(C_DOUBLE_COMPLEX),  pointer,     dimension(:,:) :: u_perturb, v_perturb
        complex(C_DOUBLE_COMPLEX),  pointer,     dimension(:,:) :: u_accumulator, v_accumulator
        type(C_PTR) :: uh_aligned_data, up_aligned_data, ua_aligned_data
        type(C_PTR) :: vh_aligned_data, vp_aligned_data, va_aligned_data

        type(C_PTR) :: uplan, vplan
    end type linear_theory_type

    ! ------------------------------------------------
    ! Tendency terms output by various physics subroutines
    ! ------------------------------------------------
    type tendencies_type
        ! 3D atmospheric field tendencies
        ! These are used by various physics parameterizations
        real,   allocatable, dimension(:,:,:) :: th,qv,qc,qi,u,v,qr,qs

        ! advection and pbl tendencies that need to be saved for the cumulus scheme
        real, allocatable, dimension(:,:,:) :: qv_adv,qv_pbl
        real, allocatable, dimension(:,:,:) :: th_lwrad, th_swrad
    end type tendencies_type

    ! ------------------------------------------------
    ! All fields needed in the domain defined in detail above
    ! ------------------------------------------------
    type, extends(linearizable_type) :: domain_type
        ! 3D atmospheric fields
        real, allocatable, dimension(:,:,:) :: w,ur,vr,wr   ! w, and u,v,w * density
        real, allocatable, dimension(:,:,:) :: u_cu, v_cu, w_cu            ! wind from convective processes
        real, allocatable, dimension(:,:,:) :: u_block, v_block, w_block   ! wind from blocking process
        real, allocatable, dimension(:,:,:) :: w_real       ! real space w on the mass grid (including U,V*dz/dx component) for output
        real, allocatable, dimension(:,:,:) :: nice,nrain   ! number concentration for ice and rain
        real, allocatable, dimension(:,:,:) :: nsnow,ngraupel! number concentration for snow and graupel (used in mp_morrison)
        real, allocatable, dimension(:,:,:) :: qgrau        ! graupel mass mixing ratio
        real, allocatable, dimension(:,:,:) :: p_inter      ! pressure on the vertical interfaces (p[:,1,:]=psfc)
        real, allocatable, dimension(:,:,:) :: z_inter      ! z height on interface levels
        real, allocatable, dimension(:)     :: z_layers, z_interface_layers
        real, allocatable, dimension(:,:,:) :: dz_inter     ! dz between interface levels
        real, allocatable, dimension(:,:,:) :: mut          ! mass in a given cell ? (pbot-ptop) used in some physics schemes
        ! 3D soil field
        real, allocatable, dimension(:,:,:) :: soil_t       ! soil temperature (nx, nsoil, ny)  [K]
        real, allocatable, dimension(:,:,:) :: soil_vwc     ! soil volumetric water content     [-]

        ! 2D fields, primarily fluxes to/from the land surface
        ! surface pressure and model top pressure
        real, allocatable, dimension(:,:)   :: psfc, ptop
        ! precipitation fluxes
        real,    allocatable, dimension(:,:):: rain, crain, snow, graupel   ! accumulated : total precip, convective precip, snow, and graupel
        real,    allocatable, dimension(:,:):: current_rain, current_snow   ! current time step rain and snow
        integer, allocatable, dimension(:,:):: rain_bucket, crain_bucket, snow_bucket, graupel_bucket  ! buckets to preserve accuracy in all variables

        ! radiative fluxes (and cloud fraction)
        real, allocatable, dimension(:,:)   :: swdown                       ! shortwave down at the surface         [W/m^2]
        real, allocatable, dimension(:,:)   :: lwdown                       ! longwave down at the surface          [W/m^2]
        real, allocatable, dimension(:,:)   :: cloudfrac                    ! total column cloud fraction           [-]
        real, allocatable, dimension(:,:)   :: lwup                         ! longwave up at/from the surface       [W/m^2]

        ! turbulent and ground heat fluxes
        real, allocatable, dimension(:,:)   :: sensible_heat, latent_heat, ground_heat       !  [W/m^2]

        ! domain parameters
        real, allocatable, dimension(:,:)   :: landmask             ! store the land-sea mask                       [-]
        real, allocatable, dimension(:,:)   :: sintheta, costheta   ! rotations about the E-W, N-S grid             [-]
        real, allocatable, dimension(:)     :: ZNU, ZNW             ! = (p-p_top)/(psfc-ptop),  (p_inter-p_top)/(psfc-ptop)

        ! land surface state and parameters primarily used by Noah LSM
        real, allocatable, dimension(:,:)   :: soil_tdeep           ! specified deep soil temperature               [K]
        real, allocatable, dimension(:,:)   :: skin_t               ! land surface skin temperature                 [K]
        real, allocatable, dimension(:,:)   :: soil_totalmoisture   ! total column soil moisture                    [kg/m^2]
        real, allocatable, dimension(:,:)   :: snow_swe             ! snow water equivalent                         [kg/m^2]
        real, allocatable, dimension(:,:)   :: canopy_water         ! canopy water content                          [kg/m^2]
        real, allocatable, dimension(:,:,:) :: vegfrac              ! specified monthly vegetation fraction         [%?]
        integer, allocatable, dimension(:,:):: soil_type            ! soil type index into SOILPARM.TBL
        integer, allocatable, dimension(:,:):: veg_type             ! vegetation type index into VEGPARM.TBL

        ! ocean surface state
        real, allocatable, dimension(:,:)   :: sst                  ! Sea surface temperature (from forcing data)

        ! surface and PBL parameter
        real, allocatable, dimension(:,:)   :: znt                  ! surface (background?) roughness length        [m]
        real, allocatable, dimension(:,:)   :: ustar                ! surface shear velocity u*                     [m/s]
        real, allocatable, dimension(:,:)   :: pbl_height           ! height of the PBL (only used with PBL=1 & LSM=1)
        real, allocatable, dimension(:,:)   :: u10, v10             ! 10m height u and v winds                      [m/s]
        real, allocatable, dimension(:,:)   :: t2m, q2m             ! 2m height air temperature                     [K]
                                                                    ! and water vapor mixing ratio                  [kg/kg]

        real, allocatable, dimension(:,:)   :: terrain_blocking     ! smoothed terrain delta for froude num calc.   [m]
        logical :: blocking_initialized                             ! flag to mark that the terrain_blocking field has been initialized

        ! current model time step length (should this be somewhere else?)
        real :: dt
        ! current model time (seconds from options%time_zero)
        type(Time_type) :: model_time

        ! online bias correction data
        real, allocatable, dimension(:,:,:) :: rain_fraction        ! seasonally varying fraction to multiple rain  [-]

        ! model specific fields
        real, allocatable, dimension(:,:,:) :: Um, Vm ! U and V on mass coordinates
        real, allocatable, dimension(:,:,:) :: T      ! real T (not potential)

        ! store all of the physics tendency terms
        type(tendencies_type) :: tend
    end type domain_type

    ! ------------------------------------------------
    ! boundary conditions type, must be linearizable so we can remove low res linear wind field
    ! ------------------------------------------------
    type, extends(linearizable_type) :: bc_type
        ! store the full high-res 3D grid for the next time step to compute dXdt fields
        ! includes high res versions of low res terrain and z
        type(domain_type)::next_domain

        ! store the timestep to the next input
        type(time_delta_t) :: dt

        ! dX_dt variables are the change in variable X between two forcing time steps
        ! wind and pressure dX_dt fields applied to full 3d grid, others applied only to boundaries
        real, allocatable, dimension(:,:,:) :: du_dt,dv_dt,dp_dt,dth_dt,dqv_dt,dqc_dt,dqi_dt,dqr_dt,dqs_dt,dqg_dt
        ! sh, lh, and pblh fields are only 2d.
        ! These are only used with LSM option 1 and are derived from forcing file
        real, allocatable, dimension(:,:)   :: dsh_dt,dlh_dt,dpblh_dt
        real, allocatable, dimension(:,:)   :: drain_dt
        ! change in shortwave and longwave at surface if read from forcing
        real, allocatable, dimension(:,:)   :: dsw_dt, dlw_dt
        ! change in sst if read from forcing file
        real, allocatable, dimension(:,:)   :: dsst_dt

        ! store the low resolution versiof of terrain and atmospheric elevations
        real,allocatable,dimension(:,:)     :: lowres_terrain
        real,allocatable,dimension(:,:,:)   :: lowres_z

        ! if we are using an external wind field, store them here temporarily...
        ! does this need to be separate from next_domain other than the nfiles attribute?
        type(wind_type)::ext_winds
    end type bc_type
end module data_structures
