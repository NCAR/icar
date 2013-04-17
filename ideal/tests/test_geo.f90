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
		
		allocate(domain%lat(100,100))
		allocate(domain%lon(100,100))
		
		call io_read2d(filename,"XLAT",tempdata)
		nx=size(tempdata,1)
		ny=size(tempdata,2)
		write(*,*) "domain raw lat nx,ny", nx,ny
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
		implicit none
		character, intent(in) :: filename
		type(bc_type),intent(inout)::bc

		integer::nx,ny,nz,i,j
		
		call io_read2d(filename,"XLAT",bc%lat)
		call io_read2d(filename,"XLONG",bc%lon)
		nx=size(bc%lat,1)
		nz=3
		ny=size(bc%lat,2)
		write(*,*) "BC nx,nz,ny", nx,nz,ny
		allocate(bc%w(ny,nz,nx))
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
	character::lo_res_file="lo.nc"
	character::hi_res_file="hi.nc"
	character::output_file="out.nc"
	type(domain_type)::domain
	type(bc_type)::bc
	integer::nx,ny,nz

	call init_hires(hi_res_file,domain)
	call init_lowres(lo_res_file,bc)
	
	call geo_LUT(domain,bc)
	nx=size(domain%w,1)
	nz=size(domain%w,2)
	ny=size(domain%w,3)
	call geo_interp(domain%w,bc%w,bc%geolut,.FALSE.,nx,nz,ny)
	call io_write3d(output_file,"W",domain%w)
	
	
end program test_geo