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
    use wind
    implicit none
    
    type(options_type) :: options
    type(domain_type)  :: domain
    integer ::i,nz,j,ny
    
    call init_options("icar_options.nml",options)
    call init_domain(options,domain)
    
    write(*,*) options%init_conditions_file,options%output_file
    write(*,*) options%ntimesteps,options%outputinterval
    call update_winds(domain,options)
    nz=size(domain%w,2)
    ny=size(domain%w,1)
    do i=nz,1,-1
        write(*,*) "W"
!       write(*,"(10f15.3)") domain%w(1,i,:)
        write(*,"(10f10.3)") domain%w(2,i,:)
        write(*,*) "U"
        write(*,"(10f10.3)") domain%u(2,i,:)
!       write(*,*) "V"
!       do j=1,ny
!           write(*,"(10f15.3)") domain%v(j,i,:)
!       enddo
    enddo
    
end program test_init
