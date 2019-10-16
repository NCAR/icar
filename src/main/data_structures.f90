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
        real, allocatable, dimension(:,:,:) :: du_dt,dv_dt,dp_dt,dth_dt,dqv_dt,dqc_dt
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

    ! ------------------------------------------------
    ! type to store integer options for each physics package
    ! ------------------------------------------------
    type physics_type
        integer::microphysics
        integer::advection
        integer::boundarylayer
        integer::landsurface
        integer::watersurface
        integer::radiation
        integer::convection
        integer::windtype
    end type physics_type


    ! ------------------------------------------------
    ! store Microphysics sensitivity options
    ! ------------------------------------------------
!! ++ trude
    type mp_options_type
        real    :: Nt_c
        real    :: TNO
        real    :: am_s
        real    :: rho_g
        real    :: av_s, bv_s, fv_s, av_i
        real    :: av_g, bv_g
        real    :: Ef_si, Ef_rs, Ef_rg, Ef_ri
        real    :: C_cubes, C_sqrd
        real    :: mu_r
        real    :: t_adjust
        logical :: Ef_rw_l, EF_sw_l

        integer :: update_interval  ! maximum number of seconds between updates
        integer :: top_mp_level     ! top model level to process in the microphysics
        real    :: local_precip_fraction    ! fraction of grid cell precip to keep local vs distributing to surrounding
    end type mp_options_type
!! -- trude

    ! ------------------------------------------------
    ! store Blocked flow options
    ! ------------------------------------------------
    type block_options_type
        real    :: blocking_contribution  ! fractional contribution of flow blocking perturbation that is added [0-1]
        real    :: smooth_froude_distance ! distance (m) over which Froude number is smoothed
        integer :: n_smoothing_passes     ! number of times the smoothing window is applied
        real    :: block_fr_max           ! max froude no at which flow is only partially blocked above, no blocking
        real    :: block_fr_min           ! min froude no at which flow is only partially blocked below, full blocking
        logical :: block_flow             ! switch to use or not use the flow blocking parameterization
    end type block_options_type

    ! ------------------------------------------------
    ! store Linear Theory options
    ! ------------------------------------------------
    type lt_options_type
        integer :: buffer                   ! number of grid cells to buffer around the domain MUST be >=1
        integer :: stability_window_size    ! window to average nsq over
        real    :: max_stability            ! limits on the calculated Brunt Vaisala Frequency
        real    :: min_stability            ! these may need to be a little narrower.
        logical :: variable_N               ! Compute the Brunt Vaisala Frequency (N^2) every time step
        logical :: smooth_nsq               ! Smooth the Calculated N^2 over vert_smooth vertical levels
        integer :: vert_smooth              ! number of model levels to smooth winds over in the vertical

        real    :: N_squared                ! static Brunt Vaisala Frequency (N^2) to use
        real    :: linear_contribution      ! fractional contribution of linear perturbation to wind field (e.g. u_hat multiplied by this)
        logical :: remove_lowres_linear     ! attempt to remove the linear mountain wave from the forcing low res model
        real    :: rm_N_squared             ! static Brunt Vaisala Frequency (N^2) to use in removing linear wind field
        real    :: rm_linear_contribution   ! fractional contribution of linear perturbation to wind field to remove from the low-res field

        real    :: linear_update_fraction   ! fraction of linear perturbation to add each time step
        logical :: spatial_linear_fields    ! use a spatially varying linear wind perturbation
        logical :: linear_mask              ! use a spatial mask for the linear wind field
        logical :: nsq_calibration          ! use a spatial mask to calibrate the nsquared (brunt vaisala frequency) field

        ! Look up table generation parameters
        real    :: dirmax, dirmin           ! minimum and maximum directions to use in the LUT (typically 0 and 2*pi)
        real    :: spdmax, spdmin           ! minimum and maximum wind speeds to use in the LU (typically 0 and ~30)
        real    :: nsqmax, nsqmin           ! minimum and maximum brunt_vaisalla frequencies (typically ~1e-8 and 1e-3)
        integer :: n_dir_values, n_nsq_values, n_spd_values ! number of LUT bins for each parameter
        real    :: minimum_layer_size       ! Minimum vertical step to permit when computing LUT.
                                            ! If model layers are thicker, substepping will be used.

        logical :: read_LUT, write_LUT      ! options to read the LUT from disk (or write it)
        character(len=MAXFILELENGTH) :: u_LUT_Filename  ! u LUT filename to write
        character(len=MAXFILELENGTH) :: v_LUT_Filename  ! v LUT filename to write
        logical :: overwrite_lt_lut         ! if true any existing LUT file will be over written

    end type lt_options_type

    ! ------------------------------------------------
    ! store Advection options
    ! ------------------------------------------------
    type adv_options_type
        logical :: boundary_buffer          ! buffer to smooth a bit to cut down on oscillations at the border if FCT is not used
        logical :: flux_corrected_transport ! use Flux Corrected Transport (FCT) to maintain stability and prevent any wild oscllations
        integer :: mpdata_order             ! accuracy order for MP_DATA advection scheme.
    end type adv_options_type


    ! ------------------------------------------------
    ! store Online Bias Correction options
    ! ------------------------------------------------
    type bias_options_type
        character(len=MAXFILELENGTH):: filename             ! file containing bias correction data
        character(len=MAXVARLENGTH) :: rain_fraction_var    ! name of variable containing the fraction to multiply rain by
    end type bias_options_type

    ! ------------------------------------------------
    ! store Land Surface Model options
    ! ------------------------------------------------
    type lsm_options_type
        character (len=MAXVARLENGTH) :: LU_Categories   ! land use categories to read from VEGPARM.tbl (e.g. "USGS")
        integer :: update_interval                      ! minimum time to let pass before recomputing LSM ~300s (it may be longer)  [s]
        ! the following categories will be set by default if an known LU_Category is used
        integer :: urban_category                       ! LU index value that equals "urban"
        integer :: ice_category
        integer :: water_category
        ! use monthly vegetation fraction data, not just a single value
        logical :: monthly_vegfrac
    end type lsm_options_type


    ! ------------------------------------------------
    ! store all model options
    ! ------------------------------------------------
    type options_type
        character (len=MAXVARLENGTH) :: version,comment

        ! file names
        character (len=MAXFILELENGTH) :: init_conditions_file, linear_mask_file, nsq_calibration_file
        character (len=MAXFILELENGTH), dimension(:), allocatable::boundary_files,ext_wind_files
        character (len=MAXFILELENGTH) :: output_file,restart_file,output_file_frequency

        ! variable names from init/BC/wind/... files
        character (len=MAXVARLENGTH) :: landvar,latvar,lonvar,uvar,ulat,ulon,vvar,vlat,vlon, &
                                        hgt_hi,lat_hi,lon_hi,ulat_hi,ulon_hi,vlat_hi,vlon_hi, &
                                        pvar,pbvar,tvar,qvvar,qcvar,qivar,qrvar,qsvar,qgvar,hgtvar, &
                                        shvar,lhvar,pblhvar,zvar,zbvar,&
                                        soiltype_var, soil_t_var,soil_vwc_var,soil_deept_var, &
                                        vegtype_var,vegfrac_var, linear_mask_var, nsq_calibration_var, &
                                        swdown_var, lwdown_var, &
                                        sst_var, rain_var, time_var, &
                                        qnivar, qnrvar    ! jh - added as optional fields

        ! Filenames for files to read various physics options from
        character(len=MAXFILELENGTH) :: mp_options_filename, lt_options_filename, adv_options_filename, &
                                        lsm_options_filename, bias_options_filename, block_options_filename
        character(len=MAXFILELENGTH) :: calendar


        ! various boolean options
        logical :: debug                ! outputs a little more information at runtime (not much at present)
        logical :: interactive          ! set to true if running at the commandline to see %complete printed
        logical :: ideal                ! this is an ideal simulation, forcing will be held constant
        logical :: readz                ! read atmospheric grid elevations from file
        logical :: readdz               ! read atm model layer thicknesses from namelist
        logical :: external_winds       ! read a high res 3d wind field from an external file (e.g. a high res WRF run)
        logical :: mean_winds           ! use only a mean wind field across the entire model domain
        logical :: mean_fields          ! use only a mean forcing field across the model boundaries
        logical :: restart              ! this is a restart run, read model conditions from a restart file
        logical :: qv_is_relative_humidity! if true the input water vapor is assumed to be relative humidity instead of mixing ratio
        logical :: qv_is_spec_humidity  ! if true the input water vapor is assumed to be specific humidity instead of mixing ratio
        logical :: t_is_potential       ! if true the input temperature is interpreted as potential temperature
        logical :: z_is_geopotential    ! if true the z variable is interpreted as geopotential height
        logical :: z_is_on_interface    ! if true the z variable is interpreted as residing at model level interfaces
        logical :: advect_density       ! properly incorporate density into the advection calculations.
                                        ! Doesn't play nice with linear winds
        logical :: high_res_soil_state  ! read the soil state from the high res input file not the low res file
        logical :: surface_io_only      ! just output surface variables to speed up run and thin output

        integer :: buffer               ! buffer to remove from all sides of the high res grid supplied
        ! various integer parameters/options
        integer :: ntimesteps           ! total number of time steps to be simulated (from the first forcing data)
        integer :: nz                   ! number of model vertical levels
        integer :: ext_winds_nfiles     ! number of extrenal wind filenames to read from namelist
        type(Time_type) :: restart_time ! Date of the restart time step
        ! integer :: restart_step         ! step in forcing data to begin running
        integer :: restart_date(6)      ! date to initialize from (y,m,d, h,m,s)
        integer :: restart_step_in_file ! step in restart file to initialize from

        ! various real parameters/options
        real :: dx                      ! grid cell width [m]
        real :: dxlow                   ! forcing model grid cell width [m]
        real :: in_dt                   ! time step between forcing inputs [s]
        type(time_delta_t) :: input_dt  ! store in_dt as a time delta object
        real :: out_dt                  ! time step between output [s]
        type(time_delta_t) :: output_dt ! store out_dt as a time delta object
        real :: outputinterval          ! time steps per output
        real :: inputinterval           ! time steps per input
        real :: smooth_wind_distance    ! distance over which to smooth the forcing wind field (m)
        logical :: time_varying_z       ! read in a new z coordinate every time step and interpolate accordingly
        real :: cfl_reduction_factor    ! amount to multiple CFL by to improve stability (typically 1)
        integer :: cfl_strictness       ! CFL method 1=3D from 1D*sqrt(3), 2=ave.3D wind*sqrt(3), 3=sum.3D wind, 4=opt3 * sqrt(3), 5 = sum(max.3d)

        ! date/time parameters
        type(Time_type) :: initial_time ! Date of the first forcing time step
        type(Time_type) :: start_time   ! Date to start running the model
        type(Time_type) :: end_time     ! End point for the model simulation

        real :: t_offset                ! offset to temperature because WRF outputs potential temperature-300
        real, allocatable, dimension(:)::dz_levels ! model layer thicknesses to be read from namelist
        logical :: use_agl_height       ! interpolate from forcing to model layers using Z above ground level, not sea level

        ! defines which physics package to be used.
        type(physics_type)::physics

        ! physics parameterization options
        logical :: use_mp_options
        type(mp_options_type)::mp_options

        logical :: use_lt_options
        type(lt_options_type) :: lt_options

        logical :: use_block_options
        type(block_options_type) :: block_options

        logical :: use_adv_options
        type(adv_options_type) :: adv_options

        logical :: use_lsm_options
        type(lsm_options_type) :: lsm_options

        logical :: use_bias_correction
        type(bias_options_type) :: bias_options

        integer :: warning_level        ! level of warnings to issue when checking options settings 0-10.
                                        ! 0  = Don't print anything
                                        ! 1  = print serious warnings
        ! (DEFAULT if debug=True)       ! 2  = print all warnings
                                        ! 3-4 ... nothing specified equivalent to 2
        ! (DEFAULT if debug=False)      ! 5  = Stop for options that are likely to break the model (print all warnings)
                                        ! 6-8... nothing specified equivalent to 5
                                        ! 9  = stop on serious warnings only
                                        ! 10 = stop on all warnings
    end type options_type
end module data_structures
