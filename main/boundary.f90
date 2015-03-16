!>------------------------------------------------------------
!!
!!	Handles reading boundary conditions from the forcing file(s)
!!  Provides necessary interpolation on to the grid. 
!!
!!	Primary entry points 
!! 		bc_init      - first call only
!! 		bc_update    - all successive calls
!!
!!  Contains options to use the 
!!		mean wind field
!!		mean boundary forcing
!!		wind field smoothing
!!		removal of low resolution model linear wind field
!!
!!  bc_init loads the restart file as necessary. 
!!
!!  Both init and update compute the exner function and density fields
!!  for the forcing step and update the wind field with linear perturbations
!!  as necessary
!!
!!	Author: Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module boundary_conditions
! ----------------------------------------------------------------------------
! 	NOTE: This module was initially written to read WRF output files as input.
! 		This should serve as a basis for any additional file types, and can be 
! 		readily modified.  At some point this could be modified to check the 
! 		type if input file being specified, and call an appropriate routine. 
! 		e.g. if options%inputtype=="WRF" then call wrf_init/update()
! ----------------------------------------------------------------------------
	use data_structures
	use io_routines,            only : io_getdims, io_read3d, io_maxDims, io_read2d, io_variable_is_present, &
									   io_write3di, io_write3d
	use wind,                   only : update_winds,balance_uvw
	use linear_theory_winds,    only : linear_perturb
	use geo,                    only : geo_interp2d, geo_interp
	use vertical_interpolation, only : vinterp, vLUT_forcing
	use output,                 only : write_domain
	
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

	integer::smoothing_window=14 !(bc%dx/domain%dx) * smooth_n_coarse
	integer::smooth_n_coarse=5
	
	public::bc_init
	public::bc_update
contains
	
! Smooth an array (written for wind but will work for anything)
! only smooths over the first (x) and second or third (y) dimension
! ydim can be specified to allow working with (x,y,z) data or (x,z,y) data
! WARNING: this is a moderately complex setup to be efficient for the ydim=3 (SLOW) case
! be careful when editting
	subroutine smooth_wind(wind,windowsize,ydim)
		implicit none
		real, intent(inout), dimension(:,:,:):: wind	! 3 dimensional wind field to be smoothed
		integer,intent(in)::windowsize              	! halfwidth-1/2 of window to smooth over
														! Specified in grid cells, (+/- windowsize)
		integer,intent(in)::ydim                    	! the dimension to use for the y coordinate
														! It can be 2, or 3 (but not 1)
		real,allocatable,dimension(:,:,:)::inputwind    ! temporary array to store the input data in
		integer::i,j,k,nx,ny,nz,startx,endx,starty,endy ! various array indices/bounds
		! intermediate sums to speed up the computation
		real,allocatable,dimension(:) :: rowsums,rowmeans
		real :: cursum
		integer :: cur_n,curcol,ncols,nrows
		
		ncols=windowsize*2+1
		nx=size(wind,1)
		ny=size(wind,2) !note, this is Z for the high-res domain (ydim=3)
		nz=size(wind,3) !note, this could be the Y or Z dimension depending on ydim
		! this is ~20MB that has to be allocated deallocated every call...
		allocate(inputwind(nx,ny,nz)) ! Can't be module level because nx,ny,nz could change between calls, 
									  ! could be part of a "smoothable" object to avoid allocate-deallocating constantly
! 		nx=nx-1
! 		if (ydim==3) then
! 			nz=nz-1
! 		endif
		inputwind=wind !make a copy so we always use the unsmoothed data when computing the smoothed data
		if (((windowsize*2+1)>nx).and.(ydim==3)) then
			write(*,*) "WARNING can not operate if windowsize*2+1 is larger than nx"
			write(*,*) "NX         = ", nx
			write(*,*) "windowsize = ", windowsize
			stop
		endif
		
		!parallelize over a slower dimension (not the slowest because it is MUCH easier this way)
		! as long as the inner loops (the array operations) are over the fastest dimension we are mostly OK
		!$omp parallel firstprivate(windowsize,nx,ny,nz,ydim), &
		!$omp private(i,j,k,startx,endx,starty,endy, rowsums,rowmeans,nrows,ncols,cursum), &
		!$omp shared(wind,inputwind)
		allocate(rowsums(nx)) !this is only used when ydim=3, so nz is really ny
		allocate(rowmeans(nx)) !this is only used when ydim=3, so nz is really ny
		nrows=windowsize*2+1
		ncols=windowsize*2+1
		!$omp do schedule(static)
		do j=1,ny
			
			! so we pre-compute the sum over rows for each column in the current window
			if (ydim==3) then
				rowsums=inputwind(1:nx,j,1)*(windowsize+1)
				do i=1,windowsize
					rowsums=rowsums+inputwind(1:nx,j,i)
				enddo
			endif
			! don't parallelize over this loop because it is much more efficient to be able to assume
			! that you ran the previous row in serial for the slow (ydim=3) case
			do k=1,nz ! note this is y for ydim=3
				! ydim=3 for the main model grid which is large and takes a long time
				if (ydim==3) then
					! if we are pinned to the top edge
					if ((k-windowsize)<=1) then
						starty=1
						endy  =k+windowsize
						rowsums=rowsums-inputwind(1:nx,j,starty)+inputwind(1:nx,j,endy)
					! if we are pinned to the bottom edge
					else if ((k+windowsize)>nz) then
						starty=k-windowsize
						endy  =nz
						rowsums=rowsums-inputwind(1:nx,j,starty-1)+inputwind(1:nx,j,endy)
					! if we are in the middle (this is the most common)
					else
						starty=k-windowsize
						endy  =k+windowsize
						rowsums=rowsums-inputwind(1:nx,j,starty-1)+inputwind(:,j,endy)
					endif
					rowmeans=rowsums/nrows
					cursum=sum(rowmeans(1:windowsize))+rowmeans(1)*(windowsize+1)
				endif
				
				do i=1,nx
					if (ydim==3) then
						! if we are pinned to the left edge
						if ((i-windowsize)<=1) then
							startx=1
							endx  =i+windowsize
							cursum=cursum-rowmeans(startx)+rowmeans(endx)
						! if we are pinned to the right edge
						else if ((i+windowsize)>nx) then
							startx=i-windowsize
							endx  =nx
							cursum=cursum-rowmeans(startx-1)+rowmeans(endx)
						! if we are in the middle (this is the most common)
						else
							startx=i-windowsize
							endx  =i+windowsize
							cursum=cursum-rowmeans(startx-1)+rowmeans(endx)
						endif

						wind(i,j,k)=cursum/ncols

						! old slower way, also creates artifacts windowsize distance from the borders
! 						startx=max(1, i-windowsize)
! 						endx  =min(nx,i+windowsize)
! 						starty=max(1, k-windowsize)
! 						endy  =min(nz,k+windowsize)
! 						wind(i,j,k)=sum(inputwind(startx:endx,j,starty:endy)) &
! 									/ ((endx-startx+1)*(endy-starty+1))
					else ! ydim==2
						! ydim=2 for the input data which is a small grid, thus cheap so we still use the slow method
						!first find the current window bounds
						startx=max(1, i-windowsize)
						endx  =min(nx,i+windowsize)
						starty=max(1, j-windowsize)
						endy  =min(ny,j+windowsize)
						! then compute the mean within that window (sum/n)
						! note, artifacts near the borders (mentioned in ydim==3) don't affect this
						! because the borders *should* be well away from the domain
						! if the forcing data are not much larger than the model domain this *could* create issues
						wind(i,j,k)=sum(inputwind(startx:endx,starty:endy,k)) &
									/ ((endx-startx+1)*(endy-starty+1))
					endif
				enddo
			enddo
		enddo
		!$omp end do
		deallocate(rowmeans,rowsums)
		!$omp end parallel
		
		deallocate(inputwind)
	end subroutine smooth_wind
	
! 	generic routine to read a low res variable (varname) from a netcdf file (filename) at the current time step (curstep)
!   then interpolate it to the high res grid either in 3D or at the boundaries only (boundary_only)
! 	applies modifications specificaly for U,V,T, and P variables
	subroutine read_var(highres, filename, varname, geolut, vlut, curstep, boundary_only, options, &
		                z_lo, z_hi, time_varying_zlut)
		implicit none
		real,dimension(:,:,:),   intent(inout):: highres
		character(len=*),        intent(in)   :: filename,varname
		type(geo_look_up_table), intent(in)   :: geolut
		type(vert_look_up_table),intent(in)   :: vlut
		integer,                 intent(in)   :: curstep
		logical,                 intent(in)   :: boundary_only
		type(options_type),      intent(in)   :: options
		real,dimension(:,:,:),   intent(in), optional :: z_lo, z_hi
		type(vert_look_up_table),intent(in), optional :: time_varying_zlut
		
		! local variables
		real,dimension(:,:,:),allocatable :: inputdata,extra_data
		integer::nx,ny,nz,i
		
		! Read the data in, should be relatively fast because we are reading a low resolution forcing file
		call io_read3d(filename,varname,inputdata,curstep)
		
		nx=size(inputdata,1)
		ny=size(inputdata,2)
		nz=size(inputdata,3)
		
		! Variable specific options
		! For wind variables run a first pass of smoothing over the low res data
		if (((varname==options%vvar).or.(varname==options%uvar)).and.(.not.options%ideal)) then
			call smooth_wind(inputdata,1,2)
			
		! For Temperature, we may need to add an offset
		else if ((varname==options%tvar).and.(options%t_offset.ne.0)) then
			inputdata=inputdata+options%t_offset
		
		! For pressure, we may need to add a base pressure offset read from pbvar
		else if ((varname==options%pvar).and.(options%pbvar.ne.'')) then
			call io_read3d(filename,options%pbvar,extra_data,curstep)
			inputdata=inputdata+extra_data
			deallocate(extra_data)
		endif
		
		! if the z-axis of the input data varies over time, we need to first interpolate
		! to the "standard" z-axis so that the hi-res vlut doesn't need to change
		if ((options%time_varying_z).and.present(time_varying_zlut)) then
			allocate(extra_data(nx,ny,nz))
			extra_data=inputdata
			call vinterp(inputdata, extra_data, time_varying_zlut, axis=3)
			deallocate(extra_data)
		endif
			
		
		! just read the low res version with out interpolating for e.g. external wind data
		if ( (nx==size(highres,1)).and.(ny==size(highres,3)).and.(nz==size(highres,2)) ) then
			highres=reshape(inputdata,[nx,nz,ny],order=[1,3,2])
			deallocate(inputdata)
		else
			! interpolate data onto the high resolution grid after re-arranging the dimensions. 
			allocate(extra_data(nx,nz,ny))
			extra_data=reshape(inputdata,[nx,nz,ny],order=[1,3,2])
		
			! first interpolate to a high res grid (temporarily stored in inputdata)
			deallocate(inputdata)
			allocate(inputdata(size(highres,1),nz,size(highres,3)))
			call geo_interp(inputdata, &
							extra_data, &
							geolut,boundary_only)
		
			! Then apply vertical interpolation on that grid
			if (varname==options%pvar) then
				highres=inputdata(:,:size(highres,2),:)
! 				call update_pressure(highres,z_lo,z_hi,filename,options,geolut)
				call update_pressure(highres,z_lo,z_hi)
			else
				call vinterp(highres, inputdata, &
							 vlut,boundary_only)
			endif
			deallocate(extra_data)
			deallocate(inputdata)
		endif
						
	end subroutine read_var

! 	generic routine to read a low res variable (varname) from a netcdf file (filename) at the current time step (curstep)
!   then interpolate it to the high res grid either in 2D
! 	Primarily used for surface variables: Sensible and latent heat fluxes, PBL height, skin temperature?
	subroutine read_2dvar(highres,filename,varname,geolut,curstep,options)
		implicit none
		real,dimension(:,:),intent(inout)::highres
		character(len=*),intent(in) :: filename,varname
		type(geo_look_up_table),intent(in) :: geolut
		type(options_type),intent(in) :: options
		integer,intent(in)::curstep
	
		real,dimension(:,:),allocatable :: inputdata
	
! 		Read the data in
		call io_read2d(filename,varname,inputdata,curstep)
! 		interpolate data onto the high resolution grid
		call geo_interp2d(highres,inputdata,geolut)
		deallocate(inputdata)
					
	end subroutine read_2dvar

	
! 	rotate winds from real space back to terrain following grid (approximately)
!   assumes a simple slope transform in u and v independantly constant w/height
	subroutine rotate_ext_wind_field(domain,ext_winds)
        implicit none
        type(domain_type),intent(inout)::domain
        type(wind_type),intent(inout)::ext_winds
		integer :: nx,ny,nz,i
	
		nx=size(domain%z,1)
		nz=size(domain%z,2)
		ny=size(domain%z,3)
		do i=1,nz
			domain%u(1:nx-1,i,:)=domain%u(1:nx-1,i,:)*ext_winds%dzdx
			domain%v(:,i,1:ny-1)=domain%v(:,i,1:ny-1)*ext_winds%dzdy
		end do
	
	end subroutine rotate_ext_wind_field

	
! 	initialize the eternal winds information (filenames, nfiles, etc) and read the initial conditions
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
		call io_getdims(ext_winds_file_list(ext_winds_curfile),options%uvar, dims)
		if (dims(1)==3) then
			ext_winds_steps_in_file=1
		else
			ext_winds_steps_in_file=dims(dims(1)+1) !dims(1) = ndims
		endif
		
		do while (ext_winds_curstep>ext_winds_steps_in_file)
			ext_winds_curfile=ext_winds_curfile+1
			if (ext_winds_curfile>ext_winds_nfiles) then
				stop "Ran out of files to process!"
			endif
			ext_winds_curstep=ext_winds_curstep-ext_winds_steps_in_file 
			!instead of setting=1, this way we can set an arbitrary starting point multiple files in
			call io_getdims(ext_winds_file_list(ext_winds_curfile),options%uvar, dims)
			if (dims(1)==3) then
				ext_winds_steps_in_file=1
			else
				ext_winds_steps_in_file=dims(dims(1)+1) !dims(1) = ndims; dims(ndims+1)=ntimesteps
			endif
		enddo
		
		write(*,*) "Initial external wind file= ",ext_winds_curfile," : step= ",ext_winds_curstep
		call read_var(domain%u, ext_winds_file_list(ext_winds_curfile), options%uvar,  &
		      bc%ext_winds%u_geo%geolut, bc%u_geo%vert_lut, ext_winds_curstep, .FALSE., options)
			  
		call read_var(domain%v,    ext_winds_file_list(ext_winds_curfile),	options%vvar,  &
		      bc%ext_winds%v_geo%geolut,bc%v_geo%vert_lut,ext_winds_curstep,.FALSE.,options)
		call rotate_ext_wind_field(domain,bc%ext_winds)
	end subroutine ext_winds_init
	
! 	remove linear theory topographic winds perturbations from the low resolution wind field. 
	subroutine remove_linear_winds(domain,bc,options,filename,curstep)
        implicit none
		type(domain_type), intent(inout) :: domain
		type(bc_type), intent(inout) :: bc
		type(options_type),intent(in) :: options
		character(len=*),intent(in) :: filename
		integer,intent(in)::curstep
		real,allocatable,dimension(:,:,:)::inputdata,extra_data
		logical :: reverse_winds=.TRUE.
		integer :: nx,ny,nz
		character(len=255) :: outputfilename
		
		! first read in the low-res U and V data directly
		! load low-res U data
		call io_read3d(filename,options%uvar,extra_data,curstep)
		nx=size(extra_data,1)
		ny=size(extra_data,2)
		nz=size(extra_data,3)
		call smooth_wind(extra_data,1,2)
		bc%u=reshape(extra_data,[nx,nz,ny],order=[1,3,2])
		deallocate(extra_data)
		
		! load low-res V data
		call io_read3d(filename,options%vvar,extra_data,curstep)
		nx=size(extra_data,1)
		ny=size(extra_data,2)
		nz=size(extra_data,3)
		call smooth_wind(extra_data,1,2)
		bc%v=reshape(extra_data,[nx,nz,ny],order=[1,3,2])
		deallocate(extra_data)
		
		! remove the low-res linear wind contribution effect
		call linear_perturb(bc,options,options%vert_smooth,reverse_winds,options%advect_density)
		
		! finally interpolate low res winds to the high resolutions grid
		nx=size(domain%u,1)
		nz=size(bc%u,2)
		ny=size(domain%u,3)
		allocate(extra_data(nx,nz,ny))
		call geo_interp(extra_data, bc%u,bc%u_geo%geolut,.FALSE.)
		call vinterp(domain%u,extra_data,bc%u_geo%vert_lut)
		deallocate(extra_data)
		
		nx=size(domain%v,1)
		nz=size(bc%v,2)
		ny=size(domain%v,3)
		allocate(extra_data(nx,nz,ny))
		call geo_interp(extra_data, bc%v,bc%v_geo%geolut,.FALSE.)
		call vinterp(domain%v,extra_data,bc%v_geo%vert_lut)
		deallocate(extra_data)
	end subroutine remove_linear_winds
	
! for test cases compute the mean winds and make them constant everywhere...
	subroutine mean_winds(domain,filename,curstep,options)
		implicit none
		type(domain_type), intent(inout) :: domain
		character(len=*),intent(in)::filename
		integer,intent(in)::curstep
		type(options_type):: options
		
		real,allocatable,dimension(:,:,:)::extra_data
		integer::nx,ny,nz
		
		nz=size(domain%u,2)
		
! 		load low-res U data
		call io_read3d(filename,options%uvar,extra_data,curstep)
		domain%u=sum(extra_data(:,:,:nz))/size(extra_data(:,:,:nz))
		deallocate(extra_data)

! 		load low-res V data
		call io_read3d(filename,options%vvar,extra_data,curstep)
		domain%v=sum(extra_data(:,:,:nz))/size(extra_data(:,:,:nz))
		deallocate(extra_data)
				
	end subroutine mean_winds
	
	subroutine check_shapes_3d(data1,data2)
		implicit none
		real,dimension(:,:,:),intent(in)::data1,data2
		integer :: i
		do i=1,3
			if (size(data1,i).ne.size(data2,i)) then
				write(*,*) "Restart file 3D dimensions don't match domain"
				stop
			endif
		enddo
	end subroutine check_shapes_3d
	
! 	if we are restarting from a given point, initialize the domain from the given restart file
	subroutine load_restart_file(domain,restart_file,time_step)
		implicit none
		type(domain_type), intent(inout) :: domain
		character(len=*),intent(in)::restart_file
		integer,optional,intent(in) :: time_step
		real,allocatable,dimension(:,:,:)::inputdata
		real,allocatable,dimension(:,:)::inputdata_2d
		integer :: timeslice
		
		if (present(time_step)) then
			timeslice=time_step
		else
			timeslice=1
		endif
		
		write(*,*) "Reading atmospheric restart data"
		call io_read3d(restart_file,"u",inputdata,timeslice)
		call check_shapes_3d(inputdata,domain%u)
		domain%u=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"v",inputdata,timeslice)
		domain%v=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"qv",inputdata,timeslice)
		domain%qv=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"qc",inputdata,timeslice)
		domain%cloud=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"qr",inputdata,timeslice)
		domain%qrain=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"qi",inputdata,timeslice)
		domain%ice=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"qs",inputdata,timeslice)
		domain%qsnow=inputdata
		deallocate(inputdata)
		if (io_variable_is_present(restart_file,"qg")) then
			call io_read3d(restart_file,"qg",inputdata,timeslice)
			domain%qgrau=inputdata
			deallocate(inputdata)
		endif
		if (io_variable_is_present(restart_file,"nr")) then
			call io_read3d(restart_file,"nr",inputdata,timeslice)
			domain%nrain=inputdata
			deallocate(inputdata)
		endif
		if (io_variable_is_present(restart_file,"ni")) then
			call io_read3d(restart_file,"ni",inputdata,timeslice)
			domain%nice=inputdata
			deallocate(inputdata)
		endif
		call io_read3d(restart_file,"p",inputdata,timeslice)
		domain%p=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"th",inputdata,timeslice)
		domain%th=inputdata
		deallocate(inputdata)
		call io_read3d(restart_file,"rho",inputdata,timeslice)
		domain%rho=inputdata
		deallocate(inputdata)
		
		call io_read2d(restart_file,"rain",inputdata_2d,timeslice)
		domain%rain=inputdata_2d
		deallocate(inputdata_2d)
		
		if (io_variable_is_present(restart_file,"soil_t")) then
			write(*,*) "Reading land surface restart data"
			call io_read3d(restart_file,"soil_t",inputdata,timeslice)
			call check_shapes_3d(inputdata,domain%soil_t)
			domain%soil_t=inputdata
			deallocate(inputdata)
			call io_read3d(restart_file,"soil_w",inputdata,timeslice)
			domain%soil_vwc=inputdata
			deallocate(inputdata)
		
			call io_read2d(restart_file,"ts",inputdata_2d,timeslice)
			domain%skin_t=inputdata_2d
			deallocate(inputdata_2d)
			call io_read2d(restart_file,"hfgs",inputdata_2d,timeslice)
			domain%ground_heat=inputdata_2d
			deallocate(inputdata_2d)
			call io_read2d(restart_file,"snw",inputdata_2d,timeslice)
			domain%snow_swe=inputdata_2d
			deallocate(inputdata_2d)
			call io_read2d(restart_file,"canwat",inputdata_2d,timeslice)
			domain%canopy_water=inputdata_2d
			deallocate(inputdata_2d)
		endif		
		
	end subroutine load_restart_file
	
! 	initialize the boundary conditions (read inital conditions, etc.)
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
		
! 		in case we are using a restart file we have some trickery to do here to find the proper file to be reading from
!       and set the current time step appropriately... should probably be moved to a subroutine. 
		curfile=1
		if (options%restart) then
			curstep=options%restart_step
		else
			curstep=1
		endif
		nfiles=options%nfiles
		allocate(file_list(nfiles))
		file_list=options%boundary_files
		call io_getdims(file_list(curfile),options%pvar, dims)
		if (dims(1)==3) then
			steps_in_file=1
		else
			steps_in_file=dims(dims(1)+1) !dims(1) = ndims
		endif
		
		if (.not.options%ideal) then
			do while (curstep>steps_in_file)
				curfile=curfile+1
				if (curfile>nfiles) then
					stop "Ran out of files to process!"
				endif
				curstep=curstep-steps_in_file !instead of setting=1, this way we can set an arbitrary starting point multiple files in
				call io_getdims(file_list(curfile),options%pvar, dims)
				if (dims(1)==3) then
					steps_in_file=1
				else
					steps_in_file=dims(dims(1)+1) !dims(1) = ndims; dims(ndims+1)=ntimesteps
				endif
			enddo
		endif
! 		smoothing_window = min(max(int((bc%dx/domain%dx) * smooth_n_coarse),1),size(domain%lat,1)/10)
		if (.not.options%ideal) then
			smoothing_window = min(max(int(options%smooth_wind_distance/domain%dx),1),size(domain%lat,1)/5)
			write(*,*) "Smoothing winds over ",smoothing_window," grid cells"
		endif
! 		load the restart file
		if (options%restart) then
			call load_restart_file(domain,options%restart_file,options%restart_step_in_file)
			if (options%external_winds) then
				call ext_winds_init(domain,bc,options)
			endif
			domain%pii=(domain%p/100000.0)**(R/cp)
	        domain%rho=domain%p/(R*domain%th*domain%pii) ! kg/m^3
			call balance_uvw(domain,options)
			
			domain%model_time=(options%restart_step-1)*options%in_dt + options%time_zero
			call write_domain(domain,options,-1)
		else
! 			else load data from the first Boundary conditions file
			boundary_value=.False.
			nx=size(domain%p,1)
			ny=size(domain%p,3)
			if (options%external_winds) then
				call ext_winds_init(domain,bc,options)
! 				call smooth_wind(domain%u,1,3)
! 				call smooth_wind(domain%v,1,3)
			elseif (options%remove_lowres_linear) then
				! remove the low-res linear wind perturbation field 
				call remove_linear_winds(domain,bc,options,file_list(curfile),curstep)
				if (.not.options%ideal)then
					call smooth_wind(domain%u,smoothing_window,3)
					call smooth_wind(domain%v,smoothing_window,3)
				endif
			elseif (options%mean_winds) then
				call mean_winds(domain,file_list(curfile),curstep,options)
			else
				call read_var(domain%u, file_list(curfile),options%uvar,  &
								bc%u_geo%geolut,bc%u_geo%vert_lut,curstep,boundary_value, &
								options)
				call read_var(domain%v, file_list(curfile),options%vvar,  &
								bc%v_geo%geolut,bc%v_geo%vert_lut,curstep,boundary_value, &
								options)
				if (.not.options%ideal)then
					call smooth_wind(domain%u,smoothing_window,3)
					call smooth_wind(domain%v,smoothing_window,3)
				endif
			endif
			call read_var(domain%p,    file_list(curfile),   options%pvar,   &
							bc%geolut, bc%vert_lut, curstep, boundary_value, &
							options,   bc%lowres_z,domain%z)
			call read_var(domain%th,   file_list(curfile),   options%tvar,   &
							bc%geolut, bc%vert_lut, curstep, boundary_value, &
							options)
			call read_var(domain%qv,   file_list(curfile),   options%qvvar,  &
							bc%geolut, bc%vert_lut, curstep, boundary_value, &
							options)
			call read_var(domain%cloud,file_list(curfile),   options%qcvar,  &
							bc%geolut, bc%vert_lut, curstep, boundary_value, &
							options)
			call read_var(domain%ice,  file_list(curfile),   options%qivar,  &
							bc%geolut, bc%vert_lut, curstep, boundary_value, &
							options)

			if (options%physics%landsurface==1) then
				call read_2dvar(domain%sensible_heat,file_list(curfile),options%shvar,  bc%geolut,curstep,options)
				call read_2dvar(domain%latent_heat,  file_list(curfile),options%lhvar,  bc%geolut,curstep,options)
				call read_2dvar(domain%pbl_height,   file_list(curfile),options%pblhvar,bc%geolut,curstep,options)
				where(domain%latent_heat<0) domain%latent_heat=0
			else
				domain%sensible_heat=0
				domain%latent_heat=0
				domain%pbl_height=0
			endif
			
! 			call update_pressure(domain%p,bc%lowres_z,domain%z)
			
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

			domain%pii=(domain%p/100000.0)**(R/cp)
	        domain%rho=domain%p/(R*domain%th*domain%pii) ! kg/m^3
			call update_winds(domain,options)
		endif

	end subroutine bc_init


	subroutine update_edges(dx_dt,d1,d2)
! 		same as update_dxdt but only for the edges of the domains for 
!       fields that are calculated internally (e.g. temperature and moisture)
		implicit none
		real,dimension(:,:,:), intent(inout) :: dx_dt
		real,dimension(:,:,:), intent(in) :: d1,d2
		integer :: nx,nz,ny,i

		nx=size(d1,1)
		nz=size(d1,2)
		ny=size(d1,3)
		do i=1,nz
			dx_dt(i,:ny,1)=d1(1,i,:) -d2(1,i,:)
			dx_dt(i,:ny,2)=d1(nx,i,:)-d2(nx,i,:)
			dx_dt(i,:nx,3)=d1(:,i,1) -d2(:,i,1)
			dx_dt(i,:nx,4)=d1(:,i,ny)-d2(:,i,ny)
		enddo
		dx_dt(:,1,3:4)=0
		dx_dt(:,nx,3:4)=0
	end subroutine update_edges
	
	
	subroutine update_dxdt(bc,domain)
! 		calculate changes between the current boundary conditions and the time step boundary conditions
! 		these are used to linearly shift all fields between the two times. 
		implicit none
		type(bc_type), intent(inout) :: bc
		type(domain_type), intent(in) :: domain
		
		bc%du_dt=bc%next_domain%u-domain%u
		bc%dv_dt=bc%next_domain%v-domain%v
! 		bc%dw_dt=bc%next_domain%w-domain%w
		bc%dp_dt=bc%next_domain%p-domain%p
! 		bc%drho_dt=bc%next_domain%rho-domain%rho
		
		bc%dsh_dt  =bc%next_domain%sensible_heat-domain%sensible_heat
		bc%dlh_dt  =bc%next_domain%latent_heat-domain%latent_heat
		bc%dpblh_dt=bc%next_domain%pbl_height-domain%pbl_height

		call update_edges(bc%dth_dt,bc%next_domain%th,domain%th)
		call update_edges(bc%dqv_dt,bc%next_domain%qv,domain%qv)
		call update_edges(bc%dqc_dt,bc%next_domain%cloud,domain%cloud)
	end subroutine update_dxdt
	
	subroutine update_pressure(pressure,z_lo,z_hi)
		!adjust the pressure field for the vertical shift between the low resolution domain
		! and the high resolution domain. 
		implicit none
		real,dimension(:,:,:), intent(inout) :: pressure
		real,dimension(:,:,:), intent(in) :: z_lo,z_hi
		real,dimension(:,:,:),allocatable::slp !sea level pressure [Pa]
		integer :: nx,ny,nz,i,j
		nx=size(pressure,1)
		nz=size(pressure,2)
		ny=size(pressure,3)
		allocate(slp(nx,nz,ny))
		!$omp parallel shared(slp,pressure, z_lo,z_hi) &
		!$omp private(i,j) firstprivate(nx,ny,nz)
		!$omp do 
		do j=1,ny
			do i=1,nz
			    slp(:,i,j) = pressure(:,i,j) / (1 - 2.25577E-5 * z_lo(:,i,j))**5.25588
				pressure(:,i,j) = slp(:,i,j) * (1 - 2.25577e-5 * z_hi(:,i,j))**5.25588
			enddo
		enddo
		!$omp end do
		!$omp end parallel
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
			call io_getdims(ext_winds_file_list(ext_winds_curfile),options%uvar, dims)
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
		call read_var(bc%next_domain%u,    ext_winds_file_list(ext_winds_curfile),options%uvar, &
		              bc%ext_winds%u_geo%geolut,bc%u_geo%vert_lut,ext_winds_curstep,use_interior,options)
		call read_var(bc%next_domain%v,    ext_winds_file_list(ext_winds_curfile),options%vvar, &
		              bc%ext_winds%v_geo%geolut,bc%v_geo%vert_lut,ext_winds_curstep,use_interior,options)
		call rotate_ext_wind_field(bc%next_domain,bc%ext_winds)
	
	end subroutine update_ext_winds
	
	
	subroutine bc_update(domain,bc,options)
		implicit none
		type(domain_type),intent(inout)::domain
		type(bc_type),intent(inout)::bc
		type(options_type),intent(in)::options
		integer,dimension(io_maxDims)::dims	!note, io_maxDims is included from io_routines.
		type(bc_type) :: newbc ! just used for updating z coordinate
		real, allocatable, dimension(:,:,:) :: zbase ! may be needed to temporarily store PHB data
		logical :: use_boundary,use_interior
		integer::i,nz,nx,ny
		! MODULE variables : curstep, curfile, nfiles, steps_in_file, file_list
		
		if (.not.options%ideal) then
			curstep=curstep+1
			do while (curstep>steps_in_file)
				curfile=curfile+1
				if (curfile>nfiles) then
					stop "Ran out of files to process!"
				endif
				curstep=curstep-steps_in_file !instead of setting=1, this way we can set an arbitrary starting point multiple files in
				call io_getdims(file_list(curfile),options%pvar, dims)
				if (dims(1)==3) then
					steps_in_file=1
				else
					steps_in_file=dims(dims(1)+1) !dims(1) = ndims; dims(ndims+1)=ntimesteps
				endif
			enddo
		endif
		use_interior=.False.
		use_boundary=.True.
		
		if (options%time_varying_z) then
			! read in the updated vertical coordinate
			if (allocated(newbc%z)) then
				deallocate(newbc%z)
			endif
			call io_read3d(file_list(curfile), options%zvar, newbc%z, curstep)
			nx=size(newbc%z,1)
			ny=size(newbc%z,2)
			nz=size(newbc%z,3)
			if (options%zvar=="PH") then
				call io_read3d(file_list(curfile),"PHB", zbase, curstep)
				newbc%z=(newbc%z+zbase) / gravity
				zbase(:,:,1:nz-1)=(newbc%z(:,:,1:nz-1) + newbc%z(:,:,2:nz))/2
				newbc%z=zbase
! 				deallocate(newbc%z)
! 				allocate(newbc%z(nx,ny,nz-1))
! 				newbc%z=zbase(:,:,:nz-1)
				deallocate(zbase)
			endif
			! now simply generate a look up table to convert the current z coordinate to the original z coordinate
			call vLUT_forcing(bc,newbc)
			! set a maximum on z so we don't try to interpolate data above mass grid
			where(newbc%vert_lut%z==nz) newbc%vert_lut%z=nz-1
			nz=size(newbc%z,3)
			! generate a new high-res z dataset as well (for pressure interpolations)
			call geo_interp(bc%lowres_z, reshape(newbc%z,[nx,nz,ny],order=[1,3,2]), bc%geolut,.False.)
			
! 			call io_write3di("vLUT_tv_z1","z",newbc%vert_lut%z(1,:,:,:))
! 			call io_write3di("vLUT_tv_z2","z",newbc%vert_lut%z(2,:,:,:))
! 			call io_write3d("vLUT_tv_w1","w",newbc%vert_lut%w(1,:,:,:))
! 			call io_write3d("vLUT_tv_w2","w",newbc%vert_lut%w(2,:,:,:))
			
		endif
		
		if (options%external_winds) then
			call update_ext_winds(bc,options)
! 			call smooth_wind(bc%next_domain%u,1,3)
! 			call smooth_wind(bc%next_domain%v,1,3)
		elseif (options%remove_lowres_linear) then
			call remove_linear_winds(bc%next_domain,bc,options,file_list(curfile),curstep)
			if (.not.options%ideal)then
				call smooth_wind(bc%next_domain%u,smoothing_window,3)
				call smooth_wind(bc%next_domain%v,smoothing_window,3)
			endif
		elseif (options%mean_winds) then
			call mean_winds(bc%next_domain,file_list(curfile),curstep,options)
		else
			call read_var(bc%next_domain%u,  file_list(curfile), options%uvar, &
						  bc%u_geo%geolut,   bc%u_geo%vert_lut,  curstep,      &
						  use_interior, options, time_varying_zlut=newbc%vert_lut)
			call read_var(bc%next_domain%v,  file_list(curfile), options%vvar, &
						  bc%v_geo%geolut,   bc%v_geo%vert_lut,  curstep,      &
						  use_interior, options, time_varying_zlut=newbc%vert_lut)
			if (.not.options%ideal)then
				call smooth_wind(bc%next_domain%u,smoothing_window,3)
				call smooth_wind(bc%next_domain%v,smoothing_window,3)
			endif
		endif
		call read_var(bc%next_domain%p,       file_list(curfile), options%pvar,   &
					  bc%geolut, bc%vert_lut, curstep, use_interior, 			  &
					  options,   bc%lowres_z, domain%z, time_varying_zlut=newbc%vert_lut)
					  
		call read_var(bc%next_domain%th,      file_list(curfile), options%tvar,   &
					  bc%geolut, bc%vert_lut, curstep, use_boundary, 			  &
					  options, time_varying_zlut=newbc%vert_lut)
					  
		call read_var(bc%next_domain%qv,      file_list(curfile), options%qvvar,  &
					  bc%geolut, bc%vert_lut, curstep, use_boundary, 			  &
					  options, time_varying_zlut=newbc%vert_lut)
					  
		call read_var(bc%next_domain%cloud,   file_list(curfile), options%qcvar,  &
					  bc%geolut, bc%vert_lut, curstep, use_boundary, 			  &
					  options, time_varying_zlut=newbc%vert_lut)
					  
		call read_var(bc%next_domain%ice,     file_list(curfile), options%qivar,  &
					  bc%geolut, bc%vert_lut, curstep, use_boundary, 			  &
					  options, time_varying_zlut=newbc%vert_lut)

		if (options%physics%landsurface==1) then
			call read_2dvar(bc%next_domain%sensible_heat,file_list(curfile),options%shvar,  bc%geolut,curstep,options)
			call read_2dvar(bc%next_domain%latent_heat,  file_list(curfile),options%lhvar,  bc%geolut,curstep,options)
			call read_2dvar(bc%next_domain%pbl_height,   file_list(curfile),options%pblhvar,bc%geolut,curstep,options)
! 			NOTE, this is a kludge to prevent the model from sucking more moisture out of the lower model layer than exists
			where(domain%latent_heat<0) domain%latent_heat=0
		else
			bc%next_domain%sensible_heat=0
			bc%next_domain%latent_heat=0
			bc%next_domain%pbl_height=0
		endif
	
		nx=size(bc%next_domain%th,1)
		nz=size(bc%next_domain%th,2)
		ny=size(bc%next_domain%th,3)
		if (options%mean_fields) then
			do i=1,nz
				bc%next_domain%th(1,i,:)  = sum(bc%next_domain%th(1,i,:))  / ny
				bc%next_domain%th(nx,i,:) = sum(bc%next_domain%th(nx,i,:)) / ny
				bc%next_domain%th(:,i,1)  = sum(bc%next_domain%th(:,i,1))  / nx
				bc%next_domain%th(:,i,ny) = sum(bc%next_domain%th(:,i,ny)) / nx
				
				bc%next_domain%qv(1,i,:)  = sum(bc%next_domain%qv(1,i,:))  / ny
				bc%next_domain%qv(nx,i,:) = sum(bc%next_domain%qv(nx,i,:)) / ny
				bc%next_domain%qv(:,i,1)  = sum(bc%next_domain%qv(:,i,1))  / nx
				bc%next_domain%qv(:,i,ny) = sum(bc%next_domain%qv(:,i,ny)) / nx

				bc%next_domain%cloud(1,i,:)  = sum(bc%next_domain%cloud(1,i,:))  / ny
				bc%next_domain%cloud(nx,i,:) = sum(bc%next_domain%cloud(nx,i,:)) / ny
				bc%next_domain%cloud(:,i,1)  = sum(bc%next_domain%cloud(:,i,1))  / nx
				bc%next_domain%cloud(:,i,ny) = sum(bc%next_domain%cloud(:,i,ny)) / nx

				bc%next_domain%ice(1,i,:)  = sum(bc%next_domain%ice(1,i,:))  / ny
				bc%next_domain%ice(nx,i,:) = sum(bc%next_domain%ice(nx,i,:)) / ny
				bc%next_domain%ice(:,i,1)  = sum(bc%next_domain%ice(:,i,1))  / nx
				bc%next_domain%ice(:,i,ny) = sum(bc%next_domain%ice(:,i,ny)) / nx
			enddo
		endif
		bc%next_domain%pii=(bc%next_domain%p/100000.0)**(R/cp)
        bc%next_domain%rho=bc%next_domain%p/(R*domain%th*bc%next_domain%pii) ! kg/m^3
		
		call update_winds(bc%next_domain,options)
		
		
		call update_dxdt(bc,domain)
	end subroutine bc_update
end module boundary_conditions
