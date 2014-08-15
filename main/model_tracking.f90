module model_tracking
! 	Module to track each version of the model developement
!   Particularly the changes to the namelist required at each step
	implicit none
	character(len=1024),allocatable,dimension(:)::versionlist,deltas
contains
	
	subroutine init_model_diffs()
		implicit none
		integer::n=9
		
		allocate(versionlist(n))
		allocate(deltas(n))
		versionlist=["0.5.1","0.5.2","0.6  ","0.7  ","0.7.1","0.7.2","0.7.3","0.8  ","0.8.1"]
		deltas=[ & 
		"Earliest version in record                                                 "// & 
		"                                       ", &
		"Added dxlow and variable name definitions pvar,tvar,qvvar,qcvar,qivar,     "// &
		"      U/V:lat/lon:high/low res         ", &
		"Added variable name definitions for sensible(shvar)/latent(lhvar) heat     "// &
		"      fluxes and PBL height(pblhvar)   ", &
		"Added input interval vs output interval timestepping, ?? also removed dz   "// &
		"      and decrease_dz                  ", &
		"Added variable name definitions for zvar and landmask (landvar), added     "// &
		"      readz:boolean,x/y:min/max:integer", &
		"Removed x/y:min/max, Added dz_levels and z_info namelist                   "// &
		"                                       ", &
		"Added advect_density: boolean use                                          "// &
		"                                       ", &
		"Vertical interpolation requires zvar (can be PH / geopotential height)     "// &
		"      Also added smooth_wind_distance  ", &
		"Added proper date tracking, requires date='yyyy/mm/dd hh:mm:ss' option     "// &
		"      in namelist.                     " &
		]
		
	end subroutine init_model_diffs
	
	subroutine print_model_diffs(version)
		implicit none
		character(len=*), intent(in) :: version
		integer :: i,j
		logical :: found_a_version=.false.
		
		write(*,*) "Model changes:"
		if (.not.allocated(versionlist)) then
			call init_model_diffs()
		endif
		
		do i=1,size(versionlist)
			if (version.eq.versionlist(i)) then
				found_a_version=.true.
				if (i<6) then
					write(*,*) " (Versions <0.7.3 may not be as reliable)"
				endif
				do j=i+1,size(versionlist)
					write(*,*) "  "//trim(versionlist(j))
					write(*,*) "    "//trim(deltas(j))
					write(*,*) " "
				enddo
			endif
		enddo
		if (.not.found_a_version) then
			write(*,*) "Unable to find a matching version"
			write(*,*) "Available version history:"
			write(*,*) " (Versions <0.7.3 may not be as reliable)"
			do j=1,size(versionlist)
				write(*,*) " "
				write(*,*) " "//trim(versionlist(j))
				write(*,*) "   "//trim(deltas(j))
			enddo
		endif
	end subroutine print_model_diffs
	
	subroutine finalize_model_diffs()
		implicit none
		if (allocated(versionlist)) then
			deallocate(versionlist)
		endif
		if (allocated(deltas)) then
			deallocate(deltas)
		endif
		
	end subroutine finalize_model_diffs
end module model_tracking