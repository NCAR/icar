!>------------------------------------------------
!! Defines model constants (e.g. gravity, and MAXFILELENGTH)
!!
!!------------------------------------------------
module icar_constants

    implicit none

    character(len=5) :: kVERSION_STRING = "2.0"

    ! string lengths
    integer, parameter :: kMAX_FILE_LENGTH = 1024
    integer, parameter :: kMAX_DIM_LENGTH  = 1024
    integer, parameter :: kMAX_NAME_LENGTH = 1024
    integer, parameter :: kMAX_ATTR_LENGTH = 1024

    !>--------------------------------------------
    ! list of integer constants to be used when accessing various arrays that track variable allocation, usage, etc. requests
    !
    ! NOTE: IF YOU ADD TO THIS LIST BE SURE TO ADD AN INTEGER TO THE kVARS STRUCTURE CONSTRUCTOR BELOW IT!
    ! This could be transitioned to an enum... but then one can't just "use, only:kVARS"...
    ! enum, bind(C)
    !   enumerator ::  u, v, w,...
    ! end enum
    ! --------------------------------------------
    type var_constants_type
        SEQUENCE    ! technically SEQUENCE just requires the compiler leave them in order,
                    ! but it can also keep compilers (e.g. ifort) from padding for alignment,
                    ! as long as there is no padding we can test last_var = sizeof(kVARS)

        integer :: u
        integer :: v
        integer :: w
        integer :: w_real
        integer :: pressure
        integer :: pressure_interface
        integer :: potential_temperature
        integer :: temperature
        integer :: water_vapor
        integer :: cloud_water
        integer :: cloud_number_concentration
        integer :: cloud_ice
        integer :: ice_number_concentration
        integer :: rain_in_air
        integer :: rain_number_concentration
        integer :: snow_in_air
        integer :: snow_number_concentration
        integer :: graupel_in_air
        integer :: graupel_number_concentration
        integer :: precipitation
        integer :: convective_precipitation
        integer :: external_precipitation
        integer :: snowfall
        integer :: graupel
        integer :: snowfall_ground
        integer :: rainfall_ground
        integer :: exner
        integer :: nsquared
        integer :: density
        integer :: z
        integer :: z_interface
        integer :: dz
        integer :: dz_interface
        integer :: cloud_fraction
        integer :: shortwave
        integer :: shortwave_direct
        integer :: shortwave_diffuse
        integer :: longwave
        integer :: vegetation_fraction
        integer :: vegetation_fraction_max
        integer :: vegetation_fraction_out
        integer :: veg_type
        integer :: mass_leaf
        integer :: mass_root
        integer :: mass_stem
        integer :: mass_wood
        integer :: soil_type
        integer :: soil_texture_1
        integer :: soil_texture_2
        integer :: soil_texture_3
        integer :: soil_texture_4
        integer :: soil_sand_and_clay
        integer :: soil_carbon_stable
        integer :: soil_carbon_fast
        integer :: lai
        integer :: sai
        integer :: crop_category
        integer :: crop_type
        integer :: date_planting
        integer :: date_harvest
        integer :: growing_season_gdd
        integer :: irr_frac_total
        integer :: irr_frac_sprinkler
        integer :: irr_frac_micro
        integer :: irr_frac_flood
        integer :: irr_eventno_sprinkler
        integer :: irr_eventno_micro
        integer :: irr_eventno_flood
        integer :: irr_alloc_sprinkler
        integer :: irr_alloc_micro
        integer :: irr_alloc_flood
        integer :: irr_evap_loss_sprinkler
        integer :: irr_amt_sprinkler
        integer :: irr_amt_micro
        integer :: irr_amt_flood
        integer :: evap_heat_sprinkler
        integer :: mass_ag_grain
        integer :: growing_degree_days
        integer :: plant_growth_stage
        integer :: net_ecosystem_exchange
        integer :: gross_primary_prod
        integer :: net_primary_prod
        integer :: apar
        integer :: photosynthesis_total
        integer :: stomatal_resist_total
        integer :: stomatal_resist_sun
        integer :: stomatal_resist_shade
        integer :: gecros_state
        integer :: canopy_water
        integer :: canopy_water_ice
        integer :: canopy_water_liquid
        integer :: canopy_vapor_pressure
        integer :: canopy_temperature
        integer :: canopy_fwet
        integer :: veg_leaf_temperature
        integer :: ground_surf_temperature
        integer :: frac_between_gap
        integer :: frac_within_gap
        integer :: ground_temperature_bare
        integer :: ground_temperature_canopy
        integer :: sensible_heat
        integer :: latent_heat
        integer :: u_10m
        integer :: v_10m
        integer :: ustar
        integer :: coeff_momentum_drag
        integer :: coeff_heat_exchange
        integer :: surface_rad_temperature
        integer :: temperature_2m
        integer :: humidity_2m
        integer :: temperature_2m_veg
        integer :: temperature_2m_bare
        integer :: mixing_ratio_2m_veg
        integer :: mixing_ratio_2m_bare
        integer :: surface_pressure
        integer :: rad_absorbed_total
        integer :: rad_absorbed_veg
        integer :: rad_absorbed_bare
        integer :: rad_net_longwave
        integer :: longwave_up
        integer :: ground_heat_flux
        integer :: evap_canopy
        integer :: evap_soil_surface
        integer :: transpiration_rate
        integer :: ch_veg
        integer :: ch_veg_2m
        integer :: ch_bare
        integer :: ch_bare_2m
        integer :: ch_under_canopy
        integer :: ch_leaf
        integer :: sensible_heat_veg
        integer :: sensible_heat_bare
        integer :: sensible_heat_canopy
        integer :: evap_heat_veg
        integer :: evap_heat_bare
        integer :: evap_heat_canopy
        integer :: transpiration_heat
        integer :: ground_heat_veg
        integer :: ground_heat_bare
        integer :: net_longwave_veg
        integer :: net_longwave_bare
        integer :: net_longwave_canopy
        integer :: runoff_surface
        integer :: runoff_subsurface
        integer :: soil_totalmoisture
        integer :: soil_deep_temperature
        integer :: water_table_depth
        integer :: water_aquifer
        integer :: storage_gw
        integer :: storage_lake
        integer :: roughness_z0
        integer :: snow_water_equivalent
        integer :: snow_water_eq_prev
        integer :: snow_albedo_prev
        integer :: snow_temperature
        integer :: snow_layer_depth
        integer :: snow_layer_ice
        integer :: snow_layer_liquid_water
        integer :: snow_age_factor
        integer :: snow_height
        integer :: snow_nlayers
        integer :: soil_water_content
        integer :: eq_soil_moisture
        integer :: smc_watertable_deep
        integer :: recharge
        integer :: recharge_deep
        integer :: soil_temperature
        integer :: skin_temperature
        integer :: sst
        integer :: land_mask
        integer :: terrain
        integer :: latitude
        integer :: longitude
        integer :: u_latitude
        integer :: u_longitude
        integer :: v_latitude
        integer :: v_longitude
        integer :: tend_qv_adv
        integer :: tend_qv_pbl
        integer :: tend_qv
        integer :: tend_th
        integer :: tend_qc
        integer :: tend_qi
        integer :: tend_qs
        integer :: tend_qr
        integer :: tend_u
        integer :: tend_v
        integer :: znu
        integer :: znw
        integer :: re_cloud
        integer :: re_ice
        integer :: re_snow
        integer :: out_longwave_rad
        integer :: longwave_cloud_forcing
        integer :: shortwave_cloud_forcing
        integer :: land_emissivity
        integer :: temperature_interface
        integer :: cosine_zenith_angle
        integer :: tend_swrad
        integer :: last_var
    end type var_constants_type


    type(var_constants_type) :: kVARS = var_constants_type(   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  &
                                                             11,  12,  13,  14,  15,  16,  17,  18,  19,  20,  &
                                                             21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  &
                                                             31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  &
                                                             41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  &
                                                             51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  &
                                                             61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  &
                                                             71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  &
                                                             81,  82,  83,  84,  85,  86,  87,  88,  89,  90,  &
                                                             91,  92,  93,  94,  95,  96,  97,  98,  99, 100,  &
                                                            101, 102, 103, 104, 105, 106, 107, 108, 109, 110,  &
                                                            111, 112, 113, 114, 115, 116, 117, 118, 119, 120,  &
                                                            121, 122, 123, 124, 125, 126, 127, 128, 129, 130,  &
                                                            131, 132, 133, 134, 135, 136, 137, 138, 139, 140,  &
                                                            141, 142, 143, 144, 145, 146, 147, 148, 149, 150,  &
                                                            151, 152, 153, 154, 155, 156, 157, 158, 159, 160,  &
                                                            161, 162, 163, 164, 165, 166, 167, 168, 169, 170,  &
                                                            171, 172, 173, 174, 175, 176, 177, 178, 179, 180,  &
                                                            181, 182, 183, 184, 185, 186, 187, 188, 189, 190,  &
                                                            191, 192, 193, 194, 195, 196, 197, 198, 199, 200)

    integer, parameter :: kINTEGER_BITS     = storage_size(kINTEGER_BITS)
    integer, parameter :: kMAX_STORAGE_VARS = storage_size(kVARS) / kINTEGER_BITS

    ! Initial number of output variables for which pointers are created
    integer, parameter :: kINITIAL_VAR_SIZE= 128

    ! Maximum number of dimensions
    ! Note this is defined in NetCDF, though not enforced (a file can have more than 1024 dimensions)
    integer, parameter :: kMAX_DIMENSIONS  = 1024

!>------------------------------------------------
!! Model constants (mostly string lengths)
!! ------------------------------------------------
    integer, parameter :: MAXFILELENGTH      =   1024  ! maximum file name length
    integer, parameter :: MAXVARLENGTH       =   1024  ! maximum variable name length
    integer, parameter :: MAXLEVELS          =    500  ! maximum number of vertical layers (should typically be ~10-20)
    integer, parameter :: MAX_NUMBER_FILES   =  50000  ! maximum number of permitted input files (probably a bit extreme)
    integer, parameter :: MAXSTRINGLENGTH    =   1024  ! maximum length of other strings (e.g. netcdf attributes)
    integer, parameter :: kMAX_STRING_LENGTH =   1024  ! maximum length of other strings (e.g. netcdf attributes)


!>------------------------------------------------
!!  Default width of coarray halos, ideally might be physics dependant (e.g. based on advection spatial order)
!! ------------------------------------------------
    integer,parameter :: kDEFAULT_HALO_SIZE = 1

!>------------------------------------------------
!! Value to accept for difference between real numbers should be as a fraction but then have to test for non-zero...
!! For some variables (cloud ice) 1e-6 is not that small, for others (pressure) it might be below precision...
!! ------------------------------------------------
    real,   parameter :: kSMALL_VALUE = 1e-6

    integer, parameter :: kMAINTAIN_LON      = 0
    integer, parameter :: kPRIME_CENTERED    = 1
    integer, parameter :: kDATELINE_CENTERED = 2
    integer, parameter :: kGUESS_LON         = 3

! ------------------------------------------------
! Physics scheme selection definitions
!
! NB: BASIC typically means "use the data from the low res model"
!     SIMPLE typically means a relatively simple formulation written for ICAR
! These could all be switched to enums too, but this makes it easy to see what number each has for the options file...
! ------------------------------------------------
    integer, parameter :: kNO_STOCHASTIC = -9999
    integer, parameter :: kCU_TIEDTKE    = 1
    integer, parameter :: kCU_SIMPLE     = 2
    integer, parameter :: kCU_KAINFR     = 3
    integer, parameter :: kCU_NSAS       = 4
    integer, parameter :: kCU_BMJ        = 5

    integer, parameter :: kMP_THOMPSON   = 1
    integer, parameter :: kMP_SB04       = 2
    integer, parameter :: kMP_MORRISON   = 3
    integer, parameter :: kMP_WSM6       = 4
    integer, parameter :: kMP_THOMP_AER  = 5

    integer, parameter :: kPBL_BASIC     = 1
    integer, parameter :: kPBL_SIMPLE    = 2
    integer, parameter :: kPBL_YSU       = 3

    integer, parameter :: kWATER_BASIC   = 1
    integer, parameter :: kWATER_SIMPLE  = 2

    integer, parameter :: kLSM_BASIC     = 1
    integer, parameter :: kLSM_SIMPLE    = 2
    integer, parameter :: kLSM_NOAH      = 3
    integer, parameter :: kLSM_NOAHMP    = 4

    integer, parameter :: kRA_BASIC      = 1
    integer, parameter :: kRA_SIMPLE     = 2
    integer, parameter :: kRA_RRTMG      = 3

    integer, parameter :: kADV_UPWIND    = 1
    integer, parameter :: kADV_MPDATA    = 2

    integer, parameter :: kWIND_LINEAR   = 1
    integer, parameter :: kCONSERVE_MASS = 2
    integer, parameter :: kITERATIVE_WINDS = 3
    integer, parameter :: kLINEAR_ITERATIVE_WINDS = 5

    integer, parameter :: kLC_LAND       = 1
    integer, parameter :: kLC_WATER      = 2

    ! mm of accumulated precip before "tipping" into the bucket
    ! only performed on output operations
    integer, parameter :: kPRECIP_BUCKET_SIZE=100

! ------------------------------------------------
! Physical Constants
! ------------------------------------------------
    real, parameter :: LH_vaporization=2260000.0 ! J/kg
    ! could be calculated as 2.5E6 + (-2112.0)*temp_degC ?
    real, parameter :: Rd  = 287.058   ! J/(kg K) specific gas constant for dry air
    real, parameter :: Rw  = 461.5     ! J/(kg K) specific gas constant for moist air
    real, parameter :: cp  = 1012.0    ! J/kg/K   specific heat capacity of moist STP air?
    real, parameter :: gravity= 9.81   ! m/s^2    gravity
    real, parameter :: pi  = 3.1415927 ! pi
    real, parameter :: stefan_boltzmann = 5.67e-8 ! the Stefan-Boltzmann constant
    real, parameter :: karman = 0.41   ! the von Karman constant
    real, parameter :: solar_constant = 1366 ! W/m^2

    ! convenience parameters for various physics packages
    real, parameter :: rovcp = Rd/cp
    real, parameter :: rovg  = Rd/gravity

    ! from wrf module_model_constants
    ! parameters for calculating latent heat as a function of temperature for
    ! vaporization
    real, parameter ::  XLV0 = 3.15E6
    real, parameter ::  XLV1 = 2370.
    ! sublimation
    real, parameter ::  XLS0 = 2.905E6
    real, parameter ::  XLS1 = 259.532

    ! saturated vapor pressure parameters (?)
    real, parameter ::  SVP1 = 0.6112
    real, parameter ::  SVP2 = 17.67
    real, parameter ::  SVP3 = 29.65
    real, parameter ::  SVPT0= 273.15

    real, parameter ::  EP1  = Rw/Rd-1.
    real, parameter ::  EP2  = Rd/Rw

end module
