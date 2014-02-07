module model_tracking
	implicit none
	character(len=255),allocatable,dimension(:)::versionlist,deltas
contains
	
	subroutine init_model_diffs()
		integer::n=2
		
		allocate(versionlist(n))
		allocate(deltas(n))
		versionlist=["0.7.2","0.7.3"]
		deltas=["","added advect_density: boolean use"]
		
	end subroutine init_model_diffs
	
	subroutine print_model_diffs(version)
		character(len=*), intent(in) :: version
		integer :: i,j
		
		write(*,*) "Model changes:"
		call init_model_diffs()
		do i=1,size(versionlist)
			if (version.eq.versionlist(i)) then
				do j=i+1,size(versionlist)
					write(*,*) "  "//trim(versionlist(j))
					write(*,*) "    "//trim(deltas(j))
				enddo
			endif
		enddo
		
	end subroutine print_model_diffs
end module model_tracking