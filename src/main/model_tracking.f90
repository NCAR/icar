!> ----------------------------------------------------------------------------
!!  A tool to track model versions and print useful changes between model versions
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module model_tracking
!   Module to track each major version of the model developement
!   Particularly the changes to the namelist required at each step
    implicit none
    character(len=1024),allocatable,dimension(:)::versionlist,deltas
contains

    !> ----------------------------------------------------------------------------
    !! Initialize the model version data structures
    !!
    !! ----------------------------------------------------------------------------
    subroutine init_model_diffs()
        implicit none
        integer :: n = 21

        allocate( versionlist(n) )
        allocate( deltas(n) )

        versionlist = [character(len=1024) :: &
                       "0.5.1","0.5.2","0.6","0.7","0.7.1","0.7.2","0.7.3","0.8","0.8.1","0.8.2", &
                       "0.9","0.9.1","0.9.2","0.9.3","0.9.4","0.9.5","1.0", "1.0.1","2.0a1","2.0a2",&
                       "2.0a3"]

        deltas = [ character(len=1024) :: &
        "Earliest version in git. ",                                                            &
        "Added dxlow and variable name definitions pvar,tvar,qvvar,qcvar,qivar,"//              &
        "      U/V:lat/lon:high/low res. ",                                                     &
        "Added variable name definitions for sensible(shvar)/latent(lhvar) heat"//              &
        "      fluxes and PBL height(pblhvar). ",                                               &
        "Added input interval vs output interval timestepping, ?? also removed dz"//            &
        "      and decrease_dz. ",                                                              &
        "Added variable name definitions for zvar and landmask (landvar), added"//              &
        "      readz:boolean,x/y:min/max:integer. ",                                            &
        "Removed x/y:min/max, Added dz_levels and z_info namelist. ",                           &
        "Added advect_density: boolean use. ",                                                  &
        "Vertical interpolation requires zvar (can be PH / geopotential height)"//              &
        "      Also added smooth_wind_distance. ",                                              &
        "Added proper date tracking, requires date='yyyy/mm/dd hh:mm:ss' option"//              &
        "      in namelist.",                                                                   &
        "Added preliminary support for running the Noah LSM. ",                                 &
        "Removed add_low_topo from options... MAJOR changes elsewhere, lots of "//              &
        "      new options (mp_options, lt_options).",                                          &
        "Added MPDATA and adv_options",                                                         &
        "Output file z-axis has been changed",                                                  &
        "Pre-1.0 release added end_date, date->forcing_start_date, forcing_file_list"   //      &
        "      lt:LUT_filename, mp:update_interval, moved vert_smooth to lt_parameters,"//      &
        "      added z_is_geopotential, and zbvar changed some defaults.",                      &
        "Added Morrison and WSM6 microphysics, and the ability to remove the low "  //          &
        "      resolution linear wind field.  Lots of smaller tweaks and bug fixes."//          &
        "      Also added online bias correction option. ",                                     &
        "Added convective wind advection and improved Linear wind LUT. ",                       &
        "Relatively stable checkpoint widely used. ",                                           &
        "Significantly improved geographic interpolation, bug fixes, better time handling. ",   &
        "Added coarray fortran support, massive overhaul internally, lots of features."//       &
        "      Temporarily removed. ",                                                          &
        "Added spatially variable dz coordinate system.",                                       &
        "Add option (and requirement) to specify output variables in namelist."                 &
        ]

    end subroutine init_model_diffs

    !> ----------------------------------------------------------------------------
    !!  Print all changes to the model since the version specified
    !!
    !! ----------------------------------------------------------------------------
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

    !> ----------------------------------------------------------------------------
    !! Deallocate module data structures
    !!
    !! ----------------------------------------------------------------------------
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
