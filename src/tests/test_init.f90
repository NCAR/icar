!>------------------------------------------------------------
!!  Test running initialization code.  Doesn't do much.
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
program test_init
    use init
    use data_structures

    implicit none

    type(options_type) :: options
    type(domain_type)  :: domain
    call init_options("icar_options.nml",options)
    call init_domain(options,domain)

    write(*,*) options%init_conditions_file,options%output_file
    write(*,*) options%ntimesteps,options%timestep,options%out_dt

    write(*,*) domain%u(0,0), domain%u(5,5)
    write(*,*) shape(domain%rain)
    write(*,*) shape(domain%u)
    write(*,*) shape(domain%cloud)
end program test_init
