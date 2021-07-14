!>------------------------------------------------------------
!! Test to confirm the geographic interpolation routines work
!! Hasn't been updated in a while but should still work
!! Requires two input files with XLAT,XLONG 2D variables
!!
!! Writes 5 output files
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module load_data
    use io_routines
    use data_structures
    implicit none
contains
    subroutine init_hires(filename,domain)
        implicit none
        character(len=*) ,intent(in):: filename
        type(domain_type),intent(inout)::domain
        real, dimension(:,:),allocatable::tempdata
        integer::nx,ny,nz

        call io_read2d(filename,"XLAT",domain%lat)
        call io_read2d(filename,"XLONG",domain%lon)
        ny=size(domain%lat,1)
        nx=size(domain%lat,2)
        nz=3

        allocate(domain%w(ny,nz,nx))
        domain%w=0

    end subroutine init_hires

    subroutine init_lowres(filename,bc)
        implicit none
        character(len=*), intent(in) :: filename
        type(bc_type),intent(inout)::bc
        integer::nx,ny,nz,i,j

        call io_read2d(filename,"XLAT",bc%lat)
        call io_read2d(filename,"XLONG",bc%lon)
        nx=size(bc%lat,1)
        nz=3
        ny=size(bc%lat,2)
        allocate(bc%w(nx,nz,ny))
        do i=1,nx
            do j=1,ny
                bc%w(i,:,j)=i+j
            enddo
        enddo
    end subroutine init_lowres
end module load_data

program test_geo
    use data_structures
    use io_routines
    use geo
    use load_data

    implicit none
    character(len=100)::lo_res_file="lo.nc"
    character(len=100)::hi_res_file="hi.nc"
    character(len=100)::output_file="out.nc"
    type(domain_type)::domain
    type(bc_type)::bc

    call init_hires(hi_res_file,domain)
    call init_lowres(lo_res_file,bc)

    call geo_LUT(domain,bc)
    call geo_interp(domain%w,bc%w,bc%geolut,.FALSE.)

    call io_write3d(output_file,"W",domain%w)
    call io_write3d("out_lo.nc","W",bc%w)
    call io_write3di("out_geolutx.nc","lut",bc%geolut%x)
    call io_write3di("out_geoluty.nc","lut",bc%geolut%y)
    call io_write3d("out_geolutw.nc","lut",bc%geolut%w)

end program test_geo
