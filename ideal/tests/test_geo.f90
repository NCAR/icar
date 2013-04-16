module load_data
	use io_routines
	use data_structures
	implicit none
contains
	subroutine init_hires(filename,domain)
		character(len=*) ,intent(in):: filename
		type(domain_type),intent(inout)::domain
		real(:,:),allocatable::tempdata
		integer::nx,ny,nz
		
		allocate(domain%lat(1:100,1:100))
		allocate(domain%lon(1:100,1:100))
		
		call io_read2d(filename,"XLAT",tempdata)
		nx=size(tempdata,1)
		ny=size(tempdata,2)
		domain%lat=tempdata(nx/2-50:nx/2+49,ny/2-50:ny/2+49)
		deallocate(tempdata)
		
		call io_read2d(filename,"XLONG",tempdata)
		domain%lat=tempdata(nx/2-50:nx/2+49,ny/2-50:ny/2+49)
		deallocate(tempdata)
		
		ny=100
		nx=100
		nz=3
		allocate(domain%w(ny,nz,nx))
		domain%w=0
		
	end subroutine init_hires
	
	subroutine init_lowres(filename,bc)
		character, intent(in) :: filename
		type(boundary_type),intent(inout)::bc
		
		
	end subroutine init_lowres
end module load_data

program test_geo
	use data_structures
	use io_routines
	use geo
	use load_data
	
	character::lo_res_file="lo.nc"
	character::hi_res_file="hi.nc"
	character::output_file="out.nc"
	type(domain_type)::domain
	type(boundary_type)::bc

	call init_hires(hi_res_file,domain)
	call init_lowres(lo_res_file,bc)
	
	call geo_LUT(domain,bc)
	
	
	
	
	
end program test_geo