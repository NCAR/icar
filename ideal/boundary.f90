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
		character(len=255),intent(in) :: filename,varname
		integer,dimension(:,:,:),intent(in) :: geolut
		integer,intent(in)::curstep
		logical, intent(in) :: boundary_only
		real,dimension(:,:,:) :: inputdata
		
		call io_getdims(filename,varname, dims)
		call io_read3d(filename,varname,inputdata,curstep)
		
		nx=dims(1)
		ny=dims(2)
		nz=dims(3)
		if (varname=="U") then
			nx=nx-1
			inputdata(1:nx,:,:)=(inputdata(1:nx,:,:)+inputdata(2:nx+1,:,:))/2
		else if (varname=="V") then
			ny=ny-1
			inputdata(:,1:ny,:)=(inputdata(:,1:ny,:)+inputdata(:,2:ny+1,:))/2
		else if (varname=="T") then
			inputdata=inputdata+300
		endif
		
		call geo_interp(highres,
						reshape(inputdata(1:nx,1:ny,1:nz),[nx,nz,ny],order=[1,3,2]), &
						bc%geolut,boundary_only)
						
	end subroutine read_var
	
	subroutine bc_init(domain,bc,options)
		implicit none
		type(domain_type),intent(inout)::domain
		type(bc_type),intent(inout)::bc
		type(options_type),intent(in)::options
		integer,dimension(io_maxDims)::dims	!note, io_maxDims is included from io_routines.
		real,dimension(:,:,:)::inputdata
		
		curfile=1
		curstep=1
		nfiles=1
		allocate(file_list(nfiles))
		file_list(curfile)=options%boundary_file
		
		call io_getdims(file_list(curfile),"P", dims)
		if dims(1)==3 then
			steps_in_file=1
		else
			steps_in_file=dims(dims(1)+1) !dims(1) = ndims
		endif
		
		call read_var(domain%u,    file_list(curfile),"U",      bc%geolut,curstep,.False.)
		call read_var(domain%v,    file_list(curfile),"V",      bc%geolut,curstep,.False.)
		call read_var(domain%p,    file_list(curfile),"P",      bc%geolut,curstep,.False.)
		call read_var(domain%th,   file_list(curfile),"T",      bc%geolut,curstep,.False.)
		call read_var(domain%qv,   file_list(curfile),"QVAPOR", bc%geolut,curstep,.False.)
		call read_var(domain%cloud,file_list(curfile),"QCLOUD", bc%geolut,curstep,.False.)
		call read_var(domain%ice,  file_list(curfile),"QICE",   bc%geolut,curstep,.False.)
		
		update_winds(domain,options)
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
		call update_edges(bc%dqcdt,bc%next_domain%qc,domain%qc)
	end subroutine update_dxdt
	
	
	subroutine bc_update(domain,bc,options)
		implicit none
		type(domain_type),intent(inout)::domain,next_domain
		type(bc_type),intent(inout)::bc
		type(options_type),intent(in)::options
		integer,dimension(io_maxDims)::dims	!note, io_maxDims is included from io_routines.
		real,dimension(:,:,:)::inputdata
		
		curstep=curstep+1
		if (curstep>steps_in_file) then
			curfile=curfile+1
			curstep=1
			call io_getdims(file_list(curfile),"P", dims)
			if dims(1)==3 then
				steps_in_file=1
			else
				steps_in_file=dims(dims(1)+1) !dims(1) = ndims
			endif
		if (curfile>nfiles) then
			stop "Ran out of files to process!"
		endif
		nx=dims(2)
		ny=dims(3)
		nz=dims(4)
		allocate(inputdata(nx,nz,ny))
		
		call read_var(bc%next_domain%u,    file_list(curfile),"U",      bc%geolut,curstep,.False.)
		call read_var(bc%next_domain%v,    file_list(curfile),"V",      bc%geolut,curstep,.False.)
		call read_var(bc%next_domain%p,    file_list(curfile),"P",      bc%geolut,curstep,.False.)
		call read_var(bc%next_domain%th,   file_list(curfile),"T",      bc%geolut,curstep,.TRUE.)
		call read_var(bc%next_domain%qv,   file_list(curfile),"QVAPOR", bc%geolut,curstep,.TRUE.)
		call read_var(bc%next_domain%cloud,file_list(curfile),"QCLOUD", bc%geolut,curstep,.TRUE.)
		call read_var(bc%next_domain%ice,  file_list(curfile),"QICE",   bc%geolut,curstep,.TRUE.)
		
		update_winds(bc%next_domain,options)
		update_dxdt(bc,domain)
	end subroutine bc_init
	