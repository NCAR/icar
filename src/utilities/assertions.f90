module assertions_mod
!! Summary: Utility for runtime checking of logical assertions.
!!
!! Compile with -DNO_ASSERTIONS to turn assertions off
!!
!! Use case 1
!! ----------
!!    Pass the optional success argument & check for false return value as an indication of assertion failure:
!!
!!    use assertions_mod, only : assert,assertions
!!    if (assertions) call assert( 2 > 1, "always true inequality", success)
!!    if (error_code/=0) call my_error_handler()
!!
!! Use case 2
!! ----------
!!    Error-terminate if the assertion fails:
!!
!!    use assertions_mod, only : assert,assertions
!!    if (assertions) call assert( 2 > 1, "always true inequality")
!!
!! code contributed by Damian Rouson Sourcery Institute
!!
    implicit none
    private
    public :: assert
    public :: assertions

! Set the USE_ASSERTIONS constant below using the C preprocessor:
!
!    gfortran -cpp -DUSE_ASSERTIONS=.false. -c assertions.f90
!
! or set the corresponding ASSERTIONS variable defined in the makefile
!
!    make ASSERTIONS=on
!
! Conditioning assertion calls on this compile-time constant enables optimizing compilers
! to eliminate assertion calls during a dead-code removal phase of optimization.

    logical, parameter :: assertions=USE_ASSERTIONS

contains
    elemental impure subroutine assert(assertion,description,success)
        use iso_fortran_env, only : error_unit
        !! Report on the truth of an assertion or error-terminate on assertion failure
        implicit none
        logical, intent(in) :: assertion
        !! Most assertions will be expressions, e.g., call assert( i>0, "positive i")
        character(len=*), intent(in) :: description
        !! Brief statement of what is being asserted
        logical, intent(out), optional :: success
        !! Optional assertion result

        if (present(success)) success=assertion

        if (.not.assertion) then
            write(error_unit,*) 'Assertion "',description,'" failed on image ',this_image()
            if (.not. present(success)) error stop
        end if
    end subroutine
end module
