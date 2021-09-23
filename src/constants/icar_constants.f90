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
        integer :: exner
        integer :: nsquared
        integer :: density
        integer :: z
        integer :: z_interface
        integer :: dz
        integer :: dz_interface
        integer :: cloud_fraction
        integer :: shortwave
        integer :: longwave
        integer :: vegetation_fraction
        integer :: veg_type
        integer :: soil_type
        integer :: lai
        integer :: canopy_water
        integer :: sensible_heat
        integer :: latent_heat
        integer :: u_10m
        integer :: v_10m
        integer :: ustar
        integer :: temperature_2m
        integer :: humidity_2m
        integer :: surface_pressure
        integer :: longwave_up
        integer :: ground_heat_flux
        integer :: soil_totalmoisture
        integer :: soil_deep_temperature
        integer :: roughness_z0
        integer :: snow_water_equivalent
        integer :: snow_height
        integer :: soil_water_content
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
        integer :: last_var
    end type var_constants_type

    type(var_constants_type) :: kVARS = var_constants_type( 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,  &
                                                           11, 12, 13, 14, 15, 16, 17, 18, 19, 20,  &
                                                           21, 22, 23, 24, 25, 26, 27, 28, 29, 30,  &
                                                           31, 32, 33, 34, 35, 36, 37, 38, 39, 40,  &
                                                           41, 42, 43, 44, 45, 46, 47, 48, 49, 50,  &
                                                           51, 52, 53, 54, 55, 56, 57, 58, 59, 60,  &
                                                           61, 62, 63, 64, 65, 66, 67, 68, 69, 70,  &
                                                           71, 72, 73, 74, 75, 76, 77, 78, 79)

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

    integer, parameter :: kRA_BASIC      = 1
    integer, parameter :: kRA_SIMPLE     = 2

    integer, parameter :: kADV_UPWIND    = 1
    integer, parameter :: kADV_MPDATA    = 2

    integer, parameter :: kWIND_LINEAR   = 1
    integer, parameter :: kCONSERVE_MASS = 2
    integer, parameter :: kITERATIVE_WINDS = 3

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
