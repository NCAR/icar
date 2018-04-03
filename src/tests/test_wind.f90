!>------------------------------------------------------------
!!  Test running ideal wind calculation.  Doesn't do much
!!  Hasn't been updated in a while, not sure it will work.
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
program test_winds
    use init
    use data_structures
    use options_interface,  only: options_t
    use domain_interface,   only: domain_t
    use wind,               only: update_winds
    implicit none

    type(options_t) :: options
    type(domain_t)  :: domain
    integer ::i,nz,j,ny

    call options%init_options("icar_options.nml")
    call domain%init(options)

    ! write(*,*) options%parameters%init_conditions_file,options%parameters%output_file
    ! write(*,*) options%parameters%ntimesteps,options%parameters%outputinterval
    call update_winds(domain,options)

    nz = size(domain%w%data_3d,2)
    ny = size(domain%w%data_3d,1)
    do i=nz,1,-1
        write(*,*) "W"
!       write(*,"(10f15.3)") domain%w(1,i,:)
        write(*,"(10f10.3)") domain%w%data_3d(2,i,:)
        write(*,*) "U"
        write(*,"(10f10.3)") domain%u%data_3d(2,i,:)
!       write(*,*) "V"
!       do j=1,ny
!           write(*,"(10f15.3)") domain%v(j,i,:)
!       enddo
    enddo

end program test_init
