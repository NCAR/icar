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
! 	manage file pointer and position in file for boundary conditions
	integer::curfile,curstep
	integer::steps_in_file,nfiles
! 	manage file pointer and position in file for external winds
	character (len=255),dimension(:),allocatable :: ext_winds_file_list
	integer::ext_winds_curfile,ext_winds_curstep
	integer::ext_winds_steps_in_file,ext_winds_nfiles

    real, parameter :: R=287.058 ! J/(kg K) specific gas constant for air
    real, parameter :: cp = 1012.0 ! specific heat capacity of moist STP air? J/kg/K
	
	public::bc_init
	public::bc_update
contains
	subroutine smooth_wind(wind,windowsize)
		real, intent(inout), dimension(:,:,:):: wind
		integer,intent(in)::windowsize
		real,allocatable,dimension(:,:,:)::inputwind
		integer::i,j,k,nx,ny,nz
		nx=size(wind,1)
		ny=size(wind,2)
		nz=size(wind,3)
		allocate(inputwind(nx,ny,nz))
		inputwind=wind
		
		!$omp parallel firstprivate(windowsize,nx,ny,nz),private(i,j,k),shared(wind,inputwind)
		!$omp do
		do k=1,nz
			do j=1+windowsize,ny-windowsize
				do i=1+windowsize,nx-windowsize
					wind(i,j,k)=sum(inputwind(i-windowsize:i+windowsize,j-windowsize:j+windowsize,k))/((windowsize*2+1)**2)
				enddo
			enddo
		enddo
		!$omp end do
		!$omp end parallel
		do i=1,windowsize
			wind(:,i,:)=wind(:,1+windowsize,:)
			wind(:,ny-i+1,:)=wind(:,ny-windowsize,:)
			wind(i,:,:)=wind(1+windowsize,:,:)
			wind(nx-i+1,:,:)=wind(nx-windowsize,:,:)
		enddo
		deallocate(inputwind)
	end subroutine smooth_wind
	
	subroutine read_var(highres,filename,varname,geolut,curstep,boundary_only,low_res_terrain)
		implicit none
		real,dimension(:,:,:),intent(inout)::highres
		character(len=*),intent(in) :: filename,varname
		type(geo_look_up_table),intent(in) :: geolut
		integer,intent(in)::curstep
		logical, intent(in) :: boundary_only
		real,optional,dimension(:,:,:) :: low_res_terrain
		real,dimension(:,:,:),allocatable :: inputdata,extra_data
		integer,dimension(io_maxDims)::dims
		integer::nx,ny,nz,subnz,i
		
! 		Read the data in
		call io_getdims(filename,varname, dims)
		call io_read3d(filename,varname,inputdata,curstep)
		
! 		note dims(1)=ndims
		nx=dims(2)
		ny=dims(3)
		nz=dims(4)
! 		if (subnz.ne.nz) then
! 			allocate(extra_data(nx,ny,nz))
! 			extra_data=inputdata
! 			deallocate(inputdata)
! 			allocate(inputdata(nx,ny,subnz))
! 			inputdata=extra_data(:,:,1:subnz)
! 			nz=subnz
! 			deallocate(extra_data)
! 		endif
! 		perform destaggering and unit conversion on specific variables
		if (varname=="U") then
			allocate(extra_data(nx,ny,nz))
			extra_data=inputdata
			deallocate(inputdata)
			nx=nx-1
			allocate(inputdata(nx,ny,nz))
			inputdata=(extra_data(1:nx,:,:)+extra_data(2:nx+1,:,:))/2
			call smooth_wind(inputdata,2)
			deallocate(extra_data)
		else if (varname=="V") then
			allocate(extra_data(nx,ny,nz))
			extra_data=inputdata
			deallocate(inputdata)
			ny=ny-1
			allocate(inputdata(nx,ny,nz))
			inputdata=(extra_data(:,1:ny,:)+extra_data(:,2:ny+1,:))/2
			call smooth_wind(inputdata,2)
			deallocate(extra_data)
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
			inputdata(:,:,1:nz)=(inputdata(:,:,1:nz)+inputdata(:,:,2:nz+1))/2
			deallocate(extra_data)
		endif
! 		interpolate data onto the high resolution grid after re-arranging the dimensions. 
		nz=min(nz,size(highres,2))
		allocate(extra_data(nx,nz,ny))
		do i=1,nz
			extra_data(:,i,:)=inputdata(:,:,i)
		enddo
! 		extra_data=reshape(inputdata,[nx,nz,ny],order=[1,3,2])
		call geo_interp(highres, &
						extra_data, &
! 						reshape(inputdata,[nx,nz,ny],order=[1,3,2]), &
						geolut,boundary_only)
		deallocate(extra_data)
		deallocate(inputdata)
						
	end subroutine read_var
	
	subroutine ext_winds_init(domain,bc,options)
		implicit none
		type(domain_type),intent(inout)::domain
		type(bc_type),intent(inout)::bc
		type(options_type),intent(in)::options
		integer,dimension(io_maxDims)::dims	!note, io_maxDims is included from io_routines.
		! MODULE variables : ext_winds_ curstep, curfile, nfiles, steps_in_file, file_list
		ext_winds_curfile=1
		if (options%restart) then
			ext_winds_curstep=options%restart_step
		else
			ext_winds_curstep=1
		endif
		ext_winds_nfiles=options%ext_winds_nfiles
		allocate(ext_winds_file_list(ext_winds_nfiles))
		ext_winds_file_list=options%ext_wind_files
		call io_getdims(ext_winds_file_list(ext_winds_curfile),"U", dims)
		if (dims(1)==3) then
			ext_winds_steps_in_file=1
		else
			ext_winds_steps_in_file=dims(dims(1)+1) !dims(1) = ndims
		endif
		
! 		ext_winds_curstep=ext_winds_curstep+1
		do while (ext_winds_curstep>ext_winds_steps_in_file)
			ext_winds_curfile=ext_winds_curfile+1
			if (ext_winds_curfile>ext_winds_nfiles) then
				stop "Ran out of files to process!"
			endif
			ext_winds_curstep=ext_winds_curstep-ext_winds_steps_in_file 
			!instead of setting=1, this way we can set an arbitrary starting point multiple files in
			call io_getdims(ext_winds_file_list(ext_winds_curfile),"U", dims)
			if (dims(1)==3) then
				ext_winds_steps_in_file=1
			else
				ext_winds_steps_in_file=dims(dims(1)+1) !dims(1) = ndims; dims(ndims+1)=ntimesteps
			endif
		enddo
		
		call read_var(domain%u,    ext_winds_file_list(ext_winds_curfile),"U",      bc%ext_winds%geolut,ext_winds_curstep,.FALSE.)
		call read_var(domain%v,    ext_winds_file_list(ext_winds_curfile),"V",      bc%ext_winds%geolut,ext_winds_curstep,.FALSE.)
		
	end subroutine ext_winds_init
	
	subroutine remove_linear_winds(domain,bc,options,filename,curstep)
        implicit none
		type(domain_type), intent(inout) :: domain
		type(bc_type), intent(inout) :: bc
		type(options_type),intent(in) :: options
		character(len=*),intent(in) :: filename
		integer,intent(in)::curstep
		integer,dimension(io_maxDims)::dims !note, io_maxDims is included from io_routines.
		real,allocatable,dimension(:,:,:)::inputdata,extra_data
		logical :: reverse_winds=.TRUE.
		integer :: nx,ny,nz,nz_output
		
		call io_getdims(filename,"U", dims)
		nx=dims(2)-1
		ny=dims(3)
		nz=dims(4)
		nz_output=size(bc%u,2)
		
		allocate(inputdata(nx,ny,nz_output))
		
! 		load low-res U data
! 		allocate(extra_data(nx+1,ny,nz))
		call io_read3d(filename,"U",extra_data,curstep)
		inputdata=(extra_data(1:nx,:,:nz_output)+extra_data(2:nx+1,:,:nz_output))/2
		call smooth_wind(inputdata,3)
		deallocate(extra_data)
		bc%u=reshape(inputdata,[nx,nz_output,ny],order=[1,3,2])
		
! 		load low-res V data
! 		allocate(extra_data(nx,ny+1,nz))
		call io_read3d(filename,"V",extra_data,curstep)
		inputdata=(extra_data(:,1:ny+1,:nz_output)+extra_data(:,2:ny+1,:nz_output))/2
		call smooth_wind(inputdata,3)
		bc%v=reshape(inputdata,[nx,nz_output,ny],order=[1,3,2])
		deallocate(extra_data,inputdata)
		
! 		remove the low-res linear wind contribution effect
		call linear_perturb(bc,reverse_winds)
! 		finally interpolate low res winds to the high resolutions grid
		call geo_interp(domain%u, bc%u,bc%geolut,.FALSE.)
		call geo_interp(domain%v, bc%v,bc%geolut,.FALSE.)
	end subroutine remove_linear_winds
	
	subroutine mean_winds(domain,filename,curstep)
		type(domain_type), intent(inout) :: domain
		character(len=*),intent(in)::filename
		integer,intent(in)::curstep
		
		real,allocatable,dimension(:,:,:)::extra_data
		integer,dimension(io_maxDims)::dims
		integer::nx,ny,nz
		
		nz=size(domain%u,2)
		
! 		load low-res U data
		call io_read3d(filename,"U",extra_data,curstep)
		domain%u=sum(extra_data(:,:,:nz))/size(extra_data)
		deallocate(extra_data)

! 		load low-res V data
		call io_read3d(filename,"V",extra_data,curstep)
		domain%v=sum(extra_data(:,:,:nz))/size(extra_data)
		deallocate(extra_data)
				
	end subroutine mean_winds
	
	subroutine load_restart_file(domain,restart_file)
		type(domain_type), intent(inout) :: domain
		character(len=*),intent(in)::restart_file
		real,allocatable,dimension(:,:,:)::inputdata
		
		call io_read3d(restart_file,"u",inputdata)
		domain%u=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"v",inputdata)
		domain%v=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"qv",inputdata)
		domain%qv=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"qc",inputdata)
		domain%cloud=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"qr",inputdata)
		domain%qrain=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"qi",inputdata)
		domain%ice=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"qs",inputdata)
		domain%qsnow=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"qg",inputdata)
		domain%qgrau=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"nr",inputdata)
		domain%nrain=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"ni",inputdata)
		domain%nice=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"p",inputdata)
		domain%p=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"th",inputdata)
		domain%th=inputdata
		deallocate(inputdata)
	end subroutine load_restart_file
	
	subroutine bc_init(domain,bc,options)
		implicit none
		type(domain_type),intent(inout)::domain
		type(bc_type),intent(inout)::bc
		type(options_type),intent(in)::options
		integer,dimension(io_maxDims)::dims	!note, io_maxDims is included from io_routines.
		real,dimension(:,:,:),allocatable::inputdata
		logical :: boundary_value
		integer::nx,ny,nz,i
		real::domainsize
		! MODULE variables : curstep, curfile, nfiles, steps_in_file, file_list
		
		curfile=1
		if (options%restart) then
			curstep=options%restart_step
		else
			curstep=1
		endif
		nfiles=options%nfiles
		allocate(file_list(nfiles))
		file_list=options%boundary_files
		call io_getdims(file_list(curfile),"P", dims)
		if (dims(1)==3) then
			steps_in_file=1
		else
			steps_in_file=dims(dims(1)+1) !dims(1) = ndims
		endif
		
		do while (curstep>steps_in_file)
			curfile=curfile+1
			if (curfile>nfiles) then
				stop "Ran out of files to process!"
			endif
			curstep=curstep-steps_in_file !instead of setting=1, this way we can set an arbitrary starting point multiple files in
			call io_getdims(file_list(curfile),"P", dims)
			if (dims(1)==3) then
				steps_in_file=1
			else
				steps_in_file=dims(dims(1)+1) !dims(1) = ndims; dims(ndims+1)=ntimesteps
			endif
		enddo
		
		if (options%restart) then
			call load_restart_file(domain,options%restart_file)
			if (options%external_winds) then
				call ext_winds_init(domain,bc,options)
			endif
		else
			boundary_value=.False.
			nx=size(domain%u,1)
			ny=size(domain%u,3)
			if (options%external_winds) then
				call ext_winds_init(domain,bc,options)
			elseif (options%remove_lowres_linear) then
				call remove_linear_winds(domain,bc,options,file_list(curfile),curstep)
			elseif (options%mean_winds) then
				call mean_winds(domain,file_list(curfile),curstep)
			else
				call read_var(domain%u,    file_list(curfile),"U",      bc%geolut,curstep,boundary_value)
				call read_var(domain%v,    file_list(curfile),"V",      bc%geolut,curstep,boundary_value)
				call smooth_wind(domain%u,8)
				call smooth_wind(domain%v,8)
			endif
			call read_var(domain%p,    file_list(curfile),"P",      bc%geolut,curstep,boundary_value)
			call read_var(domain%th,   file_list(curfile),"T",      bc%geolut,curstep,boundary_value)
			call read_var(domain%qv,   file_list(curfile),"QVAPOR", bc%geolut,curstep,boundary_value)
			call read_var(domain%cloud,file_list(curfile),"QCLOUD", bc%geolut,curstep,boundary_value)
			call read_var(domain%ice,  file_list(curfile),"QICE",   bc%geolut,curstep,boundary_value)
		
			call update_pressure(domain%p,domain%th/((100000.0/domain%p)**(R/cp)), &
								 bc%next_domain%terrain,domain%terrain)
							 
	 		nz=size(domain%th,2)
			domainsize=size(domain%th,1)*size(domain%th,3)
	 		if (options%mean_fields) then
	 			do i=1,nz
	 				domain%th(:,i,:)=sum(domain%th(:,i,:))/domainsize
	 				domain%qv(:,i,:)=sum(domain%qv(:,i,:))/domainsize
	 				domain%cloud(:,i,:)=sum(domain%cloud(:,i,:))/domainsize
	 				domain%ice(:,i,:)=sum(domain%ice(:,i,:))/domainsize
	 			enddo
	 		endif
			
			call update_winds(domain,options)
		endif

	end subroutine bc_init


	subroutine update_edges(dxdt,d1,d2)
		implicit none
		real,dimension(:,:,:), intent(inout) :: dxdt
		real,dimension(:,:,:), intent(in) :: d1,d2
		integer :: nx,nz,ny,i

		nx=size(d1,1)
		nz=size(d1,2)
		ny=size(d1,3)
		do i=1,nz
			dxdt(i,:ny,1)=d1(1,i,:) -d2(1,i,:)
			dxdt(i,:ny,2)=d1(nx,i,:)-d2(nx,i,:)
			dxdt(i,:nx,3)=d1(:,i,1) -d2(:,i,1)
			dxdt(i,:nx,4)=d1(:,i,ny)-d2(:,i,ny)
		enddo
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
	
	subroutine update_pressure(pressure,temperature,hgt_lo,hgt_hi)
		implicit none
		real,dimension(:,:,:), intent(inout) :: pressure
		real,dimension(:,:,:), intent(in) :: temperature
		real,dimension(:,:), intent(in) :: hgt_lo,hgt_hi
		real,dimension(:,:),allocatable::slp
		integer :: nx,ny,nz,i
		nx=size(pressure,1)
		nz=size(pressure,2)
		ny=size(pressure,3)
		allocate(slp(nx,ny))
		do i=1,nz
		    slp = pressure(:,i,:)/(1 - 2.25577E-5*hgt_lo)**5.25588
			pressure(:,i,:) = slp * (1 - 2.25577e-5 * hgt_hi)**5.25588
		enddo
		deallocate(slp)
	end subroutine update_pressure
	
	subroutine update_ext_winds(bc,options)
		implicit none
		type(bc_type),intent(inout)::bc
		type(options_type),intent(in)::options
		integer,dimension(io_maxDims)::dims	!note, io_maxDims is included from io_routines.
		logical :: use_boundary,use_interior
		! MODULE variables : ext_winds_ curstep, curfile, nfiles, steps_in_file, file_list
		ext_winds_curstep=ext_winds_curstep+1
		if (ext_winds_curstep>ext_winds_steps_in_file) then
			ext_winds_curfile=ext_winds_curfile+1
			ext_winds_curstep=1
			call io_getdims(ext_winds_file_list(ext_winds_curfile),"U", dims)
			if (dims(1)==3) then
				ext_winds_steps_in_file=1
			else
				ext_winds_steps_in_file=dims(dims(1)+1) !dims(1) = ndims
			endif
		endif
		
		if (ext_winds_curfile>ext_winds_nfiles) then
			stop "Ran out of files to process!"
		endif
		
		use_interior=.False.
		use_boundary=.True.
		call read_var(bc%next_domain%u,    ext_winds_file_list(ext_winds_curfile),"U",      bc%ext_winds%geolut,ext_winds_curstep,use_interior)
		call read_var(bc%next_domain%v,    ext_winds_file_list(ext_winds_curfile),"V",      bc%ext_winds%geolut,ext_winds_curstep,use_interior)
	
	end subroutine update_ext_winds
	
	
	subroutine bc_update(domain,bc,options)
		implicit none
		type(domain_type),intent(inout)::domain
		type(bc_type),intent(inout)::bc
		type(options_type),intent(in)::options
		integer,dimension(io_maxDims)::dims	!note, io_maxDims is included from io_routines.
		logical :: use_boundary,use_interior
		integer::i,nz,nx,ny
		! MODULE variables : curstep, curfile, nfiles, steps_in_file, file_list
		
		curstep=curstep+1
		do while (curstep>steps_in_file)
			curfile=curfile+1
			if (curfile>nfiles) then
				stop "Ran out of files to process!"
			endif
			curstep=curstep-steps_in_file !instead of setting=1, this way we can set an arbitrary starting point multiple files in
			call io_getdims(file_list(curfile),"P", dims)
			if (dims(1)==3) then
				steps_in_file=1
			else
				steps_in_file=dims(dims(1)+1) !dims(1) = ndims; dims(ndims+1)=ntimesteps
			endif
		enddo
		
		use_interior=.False.
		use_boundary=.True.
		if (options%external_winds) then
			call update_ext_winds(bc,options)
		elseif (options%remove_lowres_linear) then
			call remove_linear_winds(bc%next_domain,bc,options,file_list(curfile),curstep)
		elseif (options%mean_winds) then
			call mean_winds(bc%next_domain,file_list(curfile),curstep)
		else
			call read_var(bc%next_domain%u,    file_list(curfile),"U",      bc%geolut,curstep,use_interior)
			call read_var(bc%next_domain%v,    file_list(curfile),"V",      bc%geolut,curstep,use_interior)
			call smooth_wind(bc%next_domain%u,8)
			call smooth_wind(bc%next_domain%v,8)
		endif
		call read_var(bc%next_domain%p,    file_list(curfile),"P",      bc%geolut,curstep,use_interior)
		call read_var(bc%next_domain%th,   file_list(curfile),"T",      bc%geolut,curstep,use_boundary)
		call read_var(bc%next_domain%qv,   file_list(curfile),"QVAPOR", bc%geolut,curstep,use_boundary)
		call read_var(bc%next_domain%cloud,file_list(curfile),"QCLOUD", bc%geolut,curstep,use_boundary)
		call read_var(bc%next_domain%ice,  file_list(curfile),"QICE",   bc%geolut,curstep,use_boundary)
		
		call update_pressure(bc%next_domain%p,domain%th/((100000.0/domain%p)**(R/cp)), &
							 bc%next_domain%terrain,domain%terrain)
		nx=size(bc%next_domain%th,1)
		nz=size(bc%next_domain%th,2)
		ny=size(bc%next_domain%th,3)
		if (options%mean_fields) then
			do i=1,nz
				bc%next_domain%th(1,i,:)=sum(bc%next_domain%th(1,i,:))/ny
				bc%next_domain%th(nx,i,:)=sum(bc%next_domain%th(nx,i,:))/ny
				bc%next_domain%th(:,i,1)=sum(bc%next_domain%th(:,i,1))/nx
				bc%next_domain%th(:,i,ny)=sum(bc%next_domain%th(:,i,ny))/nx
				
				bc%next_domain%qv(1,i,:)=sum(bc%next_domain%qv(1,i,:))/ny
				bc%next_domain%qv(nx,i,:)=sum(bc%next_domain%qv(nx,i,:))/ny
				bc%next_domain%qv(:,i,1)=sum(bc%next_domain%qv(:,i,1))/nx
				bc%next_domain%qv(:,i,ny)=sum(bc%next_domain%qv(:,i,ny))/nx

				bc%next_domain%cloud(1,i,:)=sum(bc%next_domain%cloud(1,i,:))/ny
				bc%next_domain%cloud(nx,i,:)=sum(bc%next_domain%cloud(nx,i,:))/ny
				bc%next_domain%cloud(:,i,1)=sum(bc%next_domain%cloud(:,i,1))/nx
				bc%next_domain%cloud(:,i,ny)=sum(bc%next_domain%cloud(:,i,ny))/nx

				bc%next_domain%ice(1,i,:)=sum(bc%next_domain%ice(1,i,:))/ny
				bc%next_domain%ice(nx,i,:)=sum(bc%next_domain%ice(nx,i,:))/ny
				bc%next_domain%ice(:,i,1)=sum(bc%next_domain%ice(:,i,1))/nx
				bc%next_domain%ice(:,i,ny)=sum(bc%next_domain%ice(:,i,ny))/nx
			enddo
		endif
		
		call update_winds(bc%next_domain,options)
		call update_dxdt(bc,domain)
	end subroutine bc_update
end module boundary_conditions