module boundary_conditions
! ----------------------------------------------------------------------------
! 	NOTE: This module was initially written to read WRF output files as input.
! 		This should serve as a basis for any additional file types, and can be 
! 		readily modified.  At some point this could be modified to check the 
! 		type if input file being specified, and call an appropriate routine. 
! 		e.g. if options.inputtype=="WRF" then call wrf_init/update()
! ----------------------------------------------------------------------------
	use io_routines
	use data_structures
	use wind
	use geo
	
	implicit none
	private
! 	these could be better stored in bc_type and initialized in init_bc?
	character (len=255),dimension(:),allocatable :: file_list
	integer::curfile,curstep
	integer::steps_in_file,nfiles
	
	public::bc_init
	public::bc_update
contains
	subroutine read_var(highres,filename,varname,geolut,curstep,boundary_only)
		implicit none
		real,dimension(:,:,:),intent(out)::highres
		character(len=*),intent(in) :: filename,varname
		type(geo_look_up_table),intent(in) :: geolut
		integer,intent(in)::curstep
		logical, intent(in) :: boundary_only
		real,dimension(:,:,:),allocatable :: inputdata,extra_data
		integer,dimension(io_maxDims)::dims
		integer::nx,ny,nz
		
! 		Read the data in
		call io_getdims(filename,varname, dims)
		call io_read3d(filename,varname,inputdata,curstep)
		
		nx=dims(1)
		ny=dims(2)
		nz=dims(3)
! 		perform destaggering and unit conversion on specific variables
		if (varname=="U") then
			nx=nx-1
			inputdata(:,:,1:nx)=(inputdata(:,:,1:nx)+inputdata(:,:,2:nx+1))/2
		else if (varname=="V") then
			ny=ny-1
			inputdata(:,1:ny,:)=(inputdata(:,1:ny,:)+inputdata(:,2:ny+1,:))/2
		else if (varname=="T") then
			inputdata=inputdata+300
		else if (varname=="P") then
			call io_read3d(filename,"PB",extra_data,curstep)
			inputdata=inputdata+extra_data
			deallocate(extra_data)
		else if (varname=="PH") then
			call io_read3d(filename,"PHB",extra_data,curstep)
			inputdata=(inputdata+extra_data)/9.8
			nz=nz-1
			inputdata(1:nz,:,:)=(inputdata(1:nz,:,:)+inputdata(2:nz+1,:,:))/2
			deallocate(extra_data)
		endif
		
! 		interpolate data onto the high resolution grid after re-arranging the dimensions. 
		call geo_interp(highres, &
						reshape(inputdata(1:nz,1:ny,1:nx),[nx,nz,ny],order=[3,1,2]), &
						geolut,boundary_only)
						
	end subroutine read_var
	
	subroutine bc_init(domain,bc,options)
		implicit none
		type(domain_type),intent(inout)::domain
		type(bc_type),intent(inout)::bc
		type(options_type),intent(in)::options
		integer,dimension(io_maxDims)::dims	!note, io_maxDims is included from io_routines.
		real,dimension(:,:,:),allocatable::inputdata
		logical :: boundary_value
		
		curfile=1
		curstep=1
		nfiles=1
		allocate(file_list(nfiles))
		file_list(curfile)=options%boundary_file
		
		call io_getdims(file_list(curfile),"P", dims)
		if (dims(1)==3) then
			steps_in_file=1
		else
			steps_in_file=dims(dims(1)+1) !dims(1) = ndims
		endif
		
		boundary_value=.False.
		call read_var(domain%u,    file_list(curfile),"U",      bc%geolut,curstep,boundary_value)
		call read_var(domain%v,    file_list(curfile),"V",      bc%geolut,curstep,boundary_value)
		call read_var(domain%p,    file_list(curfile),"P",      bc%geolut,curstep,boundary_value)
		call read_var(domain%th,   file_list(curfile),"T",      bc%geolut,curstep,boundary_value)
		call read_var(domain%qv,   file_list(curfile),"QVAPOR", bc%geolut,curstep,boundary_value)
		call read_var(domain%cloud,file_list(curfile),"QCLOUD", bc%geolut,curstep,boundary_value)
		call read_var(domain%ice,  file_list(curfile),"QICE",   bc%geolut,curstep,boundary_value)
		
		call update_winds(domain,options)
	end subroutine bc_init


	subroutine update_edges(dxdt,d1,d2)
		implicit none
		real,dimension(:,:,:), intent(inout) :: dxdt
		real,dimension(:,:,:), intent(in) :: d1,d2
		integer :: nx,ny

		nx=size(d1,1)
		ny=size(d1,3)

		dxdt(:,:ny,1)=d1(1,:,:) -d2(1,:,:)
		dxdt(:,:ny,2)=d1(nx,:,:)-d2(nx,:,:)
		dxdt(:,:nx,3)=d1(:,:,1) -d2(:,:,1)
		dxdt(:,:nx,4)=d1(:,:,ny)-d2(:,:,ny)
		dxdt(:,1,3:4)=0
		dxdt(:,nx,3:4)=0
	end subroutine update_edges
	
	
	subroutine update_dxdt(bc,domain)
		implicit none
		type(bc_type), intent(inout) :: bc
		type(domain_type), intent(in) :: domain
		
		bc%dudt=bc%next_domain%u-domain%u
		bc%dvdt=bc%next_domain%v-domain%v
		bc%dwdt=bc%next_domain%w-domain%w
		bc%dpdt=bc%next_domain%p-domain%p

		call update_edges(bc%dthdt,bc%next_domain%th,domain%th)
		call update_edges(bc%dqvdt,bc%next_domain%qv,domain%qv)
		call update_edges(bc%dqcdt,bc%next_domain%cloud,domain%cloud)
	end subroutine update_dxdt
	
	
	subroutine bc_update(domain,bc,options)
		implicit none
		type(domain_type),intent(inout)::domain
		type(bc_type),intent(inout)::bc
		type(options_type),intent(in)::options
		integer,dimension(io_maxDims)::dims	!note, io_maxDims is included from io_routines.
		real,dimension(:,:,:),allocatable::inputdata
		integer::nx,ny,nz
		logical :: use_boundary,use_interior
		
		curstep=curstep+1
		if (curstep>steps_in_file) then
			curfile=curfile+1
			curstep=1
			call io_getdims(file_list(curfile),"P", dims)
			if (dims(1)==3) then
				steps_in_file=1
			else
				steps_in_file=dims(dims(1)+1) !dims(1) = ndims
			endif
		endif
		if (curfile>nfiles) then
			stop "Ran out of files to process!"
		endif
		nx=dims(2)
		ny=dims(3)
		nz=dims(4)
		allocate(inputdata(nx,nz,ny))
		
		use_interior=.False.
		use_boundary=.True.
		call read_var(bc%next_domain%u,    file_list(curfile),"U",      bc%geolut,curstep,use_interior)
		call read_var(bc%next_domain%v,    file_list(curfile),"V",      bc%geolut,curstep,use_interior)
		call read_var(bc%next_domain%p,    file_list(curfile),"P",      bc%geolut,curstep,use_interior)
		call read_var(bc%next_domain%th,   file_list(curfile),"T",      bc%geolut,curstep,use_boundary)
		call read_var(bc%next_domain%qv,   file_list(curfile),"QVAPOR", bc%geolut,curstep,use_boundary)
		call read_var(bc%next_domain%cloud,file_list(curfile),"QCLOUD", bc%geolut,curstep,use_boundary)
		call read_var(bc%next_domain%ice,  file_list(curfile),"QICE",   bc%geolut,curstep,use_boundary)
		
		call update_winds(bc%next_domain,options)
		call update_dxdt(bc,domain)
	end subroutine bc_update
end module boundary_conditions