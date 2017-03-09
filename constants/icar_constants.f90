module icar_constants

    implicit none
    character(len=5) :: kVERSION_STRING = "0.9.5"

! ------------------------------------------------
! Model constants (string lengths)
! ------------------------------------------------
    integer,parameter :: MAXFILELENGTH    =1024   ! maximum file name length
    integer,parameter :: MAXVARLENGTH     =1024   ! maximum variable name length
    integer,parameter :: MAXLEVELS        = 500   ! maximum number of vertical layers (should typically be ~10-20)
    integer,parameter :: MAX_NUMBER_FILES = 50000 ! maximum number of permitted input files (probably a bit extreme)


! ------------------------------------------------
! Value to accept for difference between real numbers should be as a fraction but then have to test for non-zero...
! ------------------------------------------------
    real,   parameter :: kSMALL_VALUE = 1e-6

! ------------------------------------------------
! Physics scheme selection definitions
!
! NB: BASIC typically means "use the data from the low res model"
!     SIMPLE typically means a relatively simple formulation written for ICAR
! ------------------------------------------------
    integer, parameter :: kCU_TIEDTKE    = 1
    integer, parameter :: kCU_SIMPLE     = 2
    integer, parameter :: kCU_KAINFR     = 3

    integer, parameter :: kMP_THOMPSON   = 1
    integer, parameter :: kMP_SB04       = 2
    integer, parameter :: kMP_MORRISON   = 3
    integer, parameter :: kMP_WSM6       = 4

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
