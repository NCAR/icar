###################################################################
# Makefile rules:
#
# <default> all: 		makes icar
#			install: 	makes and installs icar in INSTALLDIR
#			clean: 		removes objects and module files
#			allclean: 	makes clean and removes executables
#			cleanall: 	alias for allclean
#			test: 		makes various unit tests (not all work)
#			icar: 		makes the primary model
#			doc: 		make doxygen documentation in docs/html (requires doxygen)
#
# Optional setting:
#	MODE =	fast:		Enables additional optimizations that are likely to break something
# 			profile:	Enables profileing
# 			debug:		Minimal debug options enabled, still optimized
# 			debugomp:	same as debug, but OpenMP is enabled
# 			debugslow:	Optimization curtailed, highly instrumented for debugging
# 			debugompslow: same as debugslow by OpenMP is enabled
#
# Note: adding -jn will parallelize the compile itself over n processors
#
# Example:
#	make clean; make install MODE=debugompslow -j4
#
###################################################################
# Variables that need to be set by the user:
#
# INSTALLDIR : default = ~/bin/
# FFTW       : location of fftw installation    default = $home/usr/local
# LIBNETCDF	 : location of netdcdf libraries	default = compiler/machine dependant /usr/local/lib
# INCNETCDF	 : location of netcdf headers		default = compiler/machine dependant /usr/local/include
#
# Dependencies: fftw (v3), netcdf (v4)
#	FFTW is available here: http://www.fftw.org/
#		FFTW is a C library with fortran headers
#	netcdf is available here: http://www.unidata.ucar.edu/software/netcdf/
#		NB Requires the same compiler be used to compile the Fortran interface as is used to compile ICAR

###################################################################
#  Specify where you want the resulting executable installed
###################################################################
ifndef INSTALLDIR
	INSTALLDIR=~/bin/
endif

# Model directories
BUILD=build/
PHYS=physics/
IO=io/
MAIN=main/
UTIL=utilities/
CONST=constants/
OBJECTS=objects/
TESTS=tests/

###################################################################
#	Various compiler specific flags, may need to edit
###################################################################
# It is also recommended that you set :
# LD_RUN_PATH=$LD_RUN_PATH:<your-netcdf-lib-path>:<your-fftw-lib-path>
# in your environment to point to the libraries you will need so the locations will be encoded in the
# compiled binary and you don't need to set LD_LIBRARY_PATH at runtime.

########################################################################################
# These are default parameters, also tries to load from environment variables
# They are overwritten with machine specific options below if known
########################################################################################
ECHO_MOVE=/bin/echo
RM=/bin/rm
CP=/bin/cp
# doxygen only required to enable "make doc"
DOXYGEN=doxygen

ifndef FC
	FC=gfortran
endif
ifndef F90
	F90=${FC}
endif

ifeq ($(PE_ENV),CRAY)
	FC=ftn
	F90=ftn
	ECHO_MOVE=/bin/mv
endif

ifndef FFTW
	FFTW=${HOME}/usr/local
endif
FFTW_PATH = ${FFTW}
LIBFFT = -L${FFTW_PATH}/lib -lm -lfftw3
INCFFT = -I${FFTW_PATH}/include

ifdef FFTW_DIR
	FFTW_PATH = ${FFTW_DIR}
	LIBFFT = -L${FFTW_DIR} -lm -lfftw3
endif
ifdef FFTW_INC
	INCFFT = -I${FFTW_INC}
endif

ifndef NETCDF
	NETCDF=/usr/local
	ifdef NETCDF_DIR
		NETCDF=${NETCDF_DIR}
	endif
endif

NCDF_PATH = ${NETCDF}
LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
INCNETCDF = -I$(NCDF_PATH)/include

ifdef NETCDFF
       LIBNETCDF += -L$(NETCDFF)/lib
       INCNETCDF += -I$(NETCDFF)/include
endif

USE_ASSERTIONS:=.false.

ifeq ($(ASSERTIONS),)
  USE_ASSERTIONS:=.true.
endif


########################################################################################
# Try to find the machine information
########################################################################################
NODENAME := $(shell uname -n)

ifndef COMPILER
        find_gfortran=$(findstring gfortran, $(filter gfortran%,$(F90)))
        ifeq (gfortran, $(find_gfortran))
                COMPILER=gnu
        endif
	ifeq ($(F90),ifort)
		COMPILER=intel
	endif
	ifeq ($(F90),ftn)
		COMPILER=cray
	endif
endif

# get GIT version info
GIT_VERSION := $(shell git describe --long --dirty --all --always | sed -e's/heads\///')

########################################################################################
#
# Once machine specific information is entered and compiler is specified,
# now we can set up compiler specific flags (may be overwritten later if MODE is set)
#
########################################################################################

# GNU fortran
ifeq ($(COMPILER), gnu)
	# to enable OpenMP use this line... OpenMP seems to hang with more than 1 thread, probably an OpenCoarrays (v1.9.1) issue with not supporting threading
	# COMP=-fopenmp -lgomp -O2 -c -ffree-line-length-none -ftree-vectorize -fimplicit-none -funroll-loops -march=native  -fno-protect-parens # -ffast-math #-flto #
	# LINK=-fopenmp -lgomp
	COMP= -fopenmp -lgomp -O2 -mfma -c -ffree-line-length-none -ftree-vectorize -fimplicit-none -funroll-loops -march=native  -fno-protect-parens # -ffast-math #-flto #
	LINK= -fopenmp -lgomp
	PREPROC=-cpp
	MODOUTPUT=-J $(BUILD)
endif

# GNU fortran using OpenSHMEM instead of MPI... experimental
# Currently works with gcc 6.3 MPT 2.15f and a special version of opencoarrays 1.8.5 compiled in the library directory below
ifeq ($(COMPILER),oshfort)
  F90=oshfort
  FC=$F90
  # link := $(compile)
  OSH_FLAGS := -fcoarray=lib -L/glade/u/home/gutmann/usr/local/opencoarrays/openshmem/ -lcaf_openshmem
  COMP= -O2 -mfma -c -ffree-line-length-none -ftree-vectorize -fimplicit-none -funroll-loops -march=native  -fno-protect-parens $(OSH_FLAGS) # -ffast-math #-flto #
  LINK= $(OSH_FLAGS)
  PREPROC=-cpp
  MODOUTPUT=-J $(BUILD)
  # opt := -fcheck=all -Wall -fbacktrace -g -cpp -DUSE_ASSERTIONS=$(USE_ASSERTIONS) -DDEBUG_OUTPUT=$(DEBUG_OUTPUT) -J $(BUILD)/ $(OSH_FLAGS)
  # link_flags := -fcheck=all -Wall -fbacktrace -g $(OSH_FLAGS)
  ifeq ($(MODE),debug)
	COMP= -g -c -ffree-line-length-none -fimplicit-none -fcheck=all -Wall -fbacktrace $(OSH_FLAGS) # -ffast-math #-flto #
	#   opt := -Ofast -cpp -DUSE_ASSERTIONS=$(USE_ASSERTIONS) -DDEBUG_OUTPUT=$(DEBUG_OUTPUT) -J $(BUILD)/ $(OSH_FLAGS)
	#   link_flags := -Ofast $(OSH_FLAGS)
  endif
endif


# Intel fortran
ifeq ($(COMPILER), intel)
	COMP=-c -u -qopenmp -liomp5 -O3 -xHost -ftz -fpe0 # -check stack,bounds -fp-stack-check
	LINK= -qopenmp -liomp5
	PREPROC=-fpp
	MODOUTPUT=-module $(BUILD)
endif

# Cray fortran
ifeq ($(COMPILER), cray)
	COMP= -h omp vector2 -O2 -c -eI -hfp0
	LINK= -fopenmp
	PREPROC= -eZ
	MODOUTPUT= -J $(BUILD) -em
endif


# Various compiling options.  Set the MODE variable with "make MODE=debugslow" etc.
ifeq ($(MODE), debugslow)
	ifeq ($(COMPILER), intel)
		COMP= -debug -debug-parameters all -traceback -ftrapuv -g -fpe0 -c -u -check all -check noarg_temp_created -CB
		LINK=
	endif
	ifeq ($(COMPILER), gnu)
		COMP= -c -g -fbounds-check -fbacktrace -finit-real=nan -ffree-line-length-none -Wconversion
		LINK=
	endif
	ifeq ($(COMPILER), cray)
		COMP=-h noomp -c -g -h develop -m 0 -R bcsp -M 399 -hfp0
		LINK=-h noomp
		PREPROC=-eZ
		MODOUTPUT=-e m -J $(BUILD)
	endif
endif
ifeq ($(MODE), debug)
	ifeq ($(COMPILER), intel)
		COMP= -debug -c -O2 -u -check all -check noarg_temp_created -traceback -fpe0 -fast-transcendentals -xhost
		LINK=
	endif
	ifeq ($(COMPILER), gnu)
		COMP= -c -O2 -g -fbounds-check -fbacktrace -finit-real=nan -ffree-line-length-none -Wconversion
		LINK=
	endif
	ifeq ($(COMPILER), cray)
		COMP=-O2 -h noomp -c -g -hfp0
		LINK=-h noomp
		PREPROC=-eZ
		MODOUTPUT=-e m -J $(BUILD)
	endif
endif
ifeq ($(MODE), debugompslow)
	ifeq ($(COMPILER), intel)
		# COMP= -openmp -liomp5 -debug -debug-parameters all -traceback -ftrapuv -g -fpe0 -c -u -check all -check noarg_temp_created -CB
		COMP= -qopenmp -liomp5 -debug -c -u	 -ftz -fpe0 -traceback -check all -check noarg_temp_created -fp-stack-check
		LINK= -qopenmp -liomp5
	endif
	ifeq ($(COMPILER), gnu)
		COMP= -fopenmp -lgomp -c -g -fbounds-check -fbacktrace -finit-real=nan -ffree-line-length-none -ffpe-trap=invalid -Wconversion
		LINK= -fopenmp -lgomp
	endif
	ifeq ($(COMPILER), cray)
		COMP=-c -g -m 0 -h develop -R bcsp -M 399 -hfp0
		LINK=
		PREPROC=-eZ
		MODOUTPUT=-e m -J $(BUILD)
	endif
endif
ifeq ($(MODE), debugomp)
	ifeq ($(COMPILER), intel)
		COMP= -qopenmp -liomp5 -debug -c -O3 -u -traceback -fpe0 -ftz -xHost # -fast-transcendentals -check all -check noarg_temp_created -fpe0
		LINK= -qopenmp -liomp5
	endif
	ifeq ($(COMPILER), gnu)
		COMP= -fopenmp -lgomp -c -O2 -g -fbounds-check -fbacktrace -finit-real=nan -ffree-line-length-none -fimplicit-none -Wconversion # -ftree-vectorize -funroll-loops -march=native  -fno-protect-parens
		LINK= -fopenmp -lgomp
		 # -Wall -Wno-unused-variable -Wno-unused-dummy-argument
	endif
	ifeq ($(COMPILER), cray)
		COMP=-O1 -c -g -hfp0
		LINK=
		PREPROC=-eZ
		MODOUTPUT=-e m -J $(BUILD)
	endif
endif

PROF=
ifeq ($(MODE), profile)
	ifeq ($(COMPILER), intel)
		PROF=-pg -debug inline-debug-info -shared-intel
		COMP=-c -u -qopenmp -liomp5 -O3 -xHost -ftz #because -fast includes -static # not available in ifort <13 -align array64byte
	endif
	ifeq ($(COMPILER), gnu)
		PROF=-g -p
	endif
endif
ifeq ($(MODE), fast) # WARNING -ipo (included in -fast) takes forever for very little gain, and this may be unstable
	ifeq ($(F90), ifort)
		COMP=-c -u -openmp -liomp5 -fast -ftz #-fast-transcendentals # not available in ifort <13: -align array64byte
	endif
endif


LINKER=$(F90)
ifeq ($(COMPILER), gnu)
	ifeq ($(CAF_MODE), single)
		CAF_FLAG=-fcoarray=single
	else
		ifeq ($(F90), caf)
			CAF_FLAG=-fcoarray=lib
		else
			ifndef CAF_DIR
				CAF_DIR=${HOME}/usr/local/lib
			endif
			CAF_FLAG=-fcoarray=lib
			CAF_LINK=-L${CAF_DIR} -lcaf_mpi
			ifneq ($(F90), ftn)
				LINKER=mpif90
			endif
		endif
	endif
endif
# Intel fortran
ifeq ($(COMPILER), intel)
	ifeq ($(CAF_MODE), single)
		CAF_FLAG=-coarray=single
	else
		CAF_FLAG=-coarray=distributed
		LINKER=mpiifort
	endif
endif

ASSERT=-DUSE_ASSERTIONS=$(USE_ASSERTIONS)

###################################################################
###################################################################
#
# Should not need to edit anything below this line
#
###################################################################
###################################################################
# copy required libraries into a directory accessible on compute nodes and set LD_RUN_PATH e.g.
# export LD_RUN_PATH=$LD_RUN_PATH:/path/to/netcdf/libraries/lib:/path/to/fftw/libraries/lib
LFLAGS=$(LINK) $(PROF) $(CAF_FLAG) $(CAF_LINK) $(LIBNETCDF) $(LIBFFT)
FFLAGS=$(COMP) $(PROF) $(CAF_FLAG) $(PREPROC) -DVERSION=\"$(GIT_VERSION)\" $(INCNETCDF) $(INCFFT) $(MODOUTPUT) $(ASSERT)

$(info $$NODENAME    = ${NODENAME})
$(info $$FC          = ${F90})
$(info $$FFTW_PATH   = ${FFTW_PATH})
$(info $$NCDF_PATH   = ${NCDF_PATH})
$(info $$GIT_VERSION = ${GIT_VERSION})
$(info $$COMP        = ${COMP})
$(info $$LINK        = ${LINK})
$(info $$MODE        = ${MODE})



OBJS=	\
		$(BUILD)init.o				\
		$(BUILD)time_step.o			\
		$(BUILD)debug_utils.o		\
		$(BUILD)icar_constants.o	\
		$(BUILD)wrf_constants.o		\
		$(BUILD)model_tracking.o	\
		$(BUILD)string.o        	\
		$(BUILD)time_delta_obj.o	\
		$(BUILD)time_obj.o      	\
		$(BUILD)time_h.o			\
		$(BUILD)timer_h.o			\
		$(BUILD)timer_obj.o			\
		$(BUILD)opt_types.o			\
		$(BUILD)time_io.o       	\
		$(BUILD)restart_h.o			\
		$(BUILD)restart.o	  		\
		$(BUILD)options_h.o			\
		$(BUILD)options_obj.o   	\
		$(BUILD)meta_data_h.o		\
		$(BUILD)meta_data_obj.o 	\
		$(BUILD)variable_h.o		\
		$(BUILD)variable_obj.o  	\
		$(BUILD)variable_dict_h.o	\
		$(BUILD)variable_dict_obj.o	\
		$(BUILD)grid_h.o			\
		$(BUILD)grid_obj.o			\
		$(BUILD)exchangeable_h.o	\
		$(BUILD)exchangeable_obj.o	\
		$(BUILD)domain_h.o			\
		$(BUILD)domain_obj.o		\
		$(BUILD)boundary_h.o		\
		$(BUILD)boundary_obj.o		\
		$(BUILD)assertions.o		\
		$(BUILD)atm_utilities.o		\
		$(BUILD)pbl_utilities.o		\
		$(BUILD)co_utilities.o		\
		$(BUILD)default_output_metadata.o\
		$(BUILD)output_h.o			\
		$(BUILD)output_obj.o		\
		$(BUILD)wind.o				\
		$(BUILD)linear_winds.o		\
		$(BUILD)fftshift.o			\
		$(BUILD)fftw.o				\
		$(BUILD)lt_lut_io.o			\
		$(BUILD)pbl_driver.o		\
		$(BUILD)pbl_simple.o		\
		$(BUILD)pbl_ysu.o			\
		$(BUILD)ra_driver.o			\
		$(BUILD)ra_simple.o			\
		$(BUILD)ra_rrtmg_lw.o       \
		$(BUILD)ra_rrtmg_sw.o       \
		$(BUILD)ra_clWRF_support.o  \
		$(BUILD)mp_driver.o			\
		$(BUILD)mp_simple.o			\
		$(BUILD)mp_thompson.o		\
		$(BUILD)mp_thompson_aer.o	\
		$(BUILD)mp_wsm6.o			\
		$(BUILD)mp_wsm3.o			\
		$(BUILD)cu_driver.o			\
		$(BUILD)cu_tiedtke.o		\
		$(BUILD)cu_nsas.o			\
		$(BUILD)cu_bmj.o			\
		$(BUILD)advection_driver.o	\
		$(BUILD)advect.o			\
		$(BUILD)adv_mpdata.o        \
		$(BUILD)lsm_driver.o		\
		$(BUILD)lsm_noahdrv.o		\
		$(BUILD)lsm_noahlsm.o		\
		$(BUILD)lsm_noahmpdrv.o  	\
		$(BUILD)lsm_noahmplsm.o  	\
		$(BUILD)lsm_noahmp_glacier.o \
		$(BUILD)lsm_noahmp_gecros.o \
		$(BUILD)water_simple.o		\
		$(BUILD)water_lake.o		\
		$(BUILD)geo_reader.o		\
		$(BUILD)io_routines.o		\
		$(BUILD)data_structures.o	\
		$(BUILD)vinterp.o			\
		$(BUILD)array_utilities.o

CAF_TEST_EXECUTABLES = 	caf_options_test 				\
						caf_init_test 					\
						caf_domain_init_test			\
						caf_bcast_test					\
						caf_write_domain_test			\
						caf_options_test				\
						caf_boundary_test				\
						caf_boundary_init_test			\
						caf_init-write_test				\
						caf_other_image_grids_test		\
						caf_one_d_decomposition_test	\
						caf_linear_winds_setup_test

TEST_EXECUTABLES= 	fftshift_test			\
	  				calendar_test			\
	  				fftw_test			\
					point_in_on_test		\
	  				array_util_test 		\
					variable_dict_test		\
					caf_threads_test		\
					time_obj_test			\
					$(CAF_TEST_EXECUTABLES)

# Pre-Coarray Tests that we may want to re-implement
# mpdata_test			\	(should be easy)
# blocking_test			\	(may be easy once linear winds is re-implemented)


###################################################################
#	User facing rules
###################################################################

icar:${OBJS} $(BUILD)driver.o
	${LINKER} $^ -o icar ${LFLAGS}
	make move_i

all:icar test

install:icar
	mkdir -p ${INSTALLDIR}
	${CP} icar ${INSTALLDIR}

move_i:
	$(ECHO_MOVE) *.i ${BUILD} 2>/dev/null || true

clean:
	${RM} $(BUILD)*.o $(BUILD)*.mod $(BUILD)*.smod *.i *.lst docs/doxygen_sqlite3.db 2>/dev/null ||:

allclean:cleanall

cleanall:clean
	${RM} icar $(TEST_EXECUTABLES) 2>/dev/null ||:

test: $(TEST_EXECUTABLES) move_i

caf_tests:	$(CAF_TEST_EXECUTABLES)

doc:
	doxygen docs/doxygenConfig

###################################################################
#	test cases
###################################################################
fftw_test: $(BUILD)test_fftw.o $(OBJS)
	${LINKER} $^ -o $@ ${LFLAGS}

fftshift_test: $(BUILD)test_fftshift.o $(OBJS)
	${LINKER} $^ -o $@ ${LFLAGS}

time_obj_test: $(BUILD)test_time_obj.o $(BUILD)time_obj.o $(BUILD)time_h.o $(BUILD)time_delta_obj.o $(BUILD)co_utilities.o
	${LINKER} $^ -o $@ ${LFLAGS}

calendar_test: $(BUILD)test_calendar.o $(BUILD)time_obj.o $(BUILD)time_h.o $(BUILD)time_delta_obj.o $(BUILD)co_utilities.o
	${LINKER} $^ -o $@ ${LFLAGS}

mpdata_test: $(BUILD)test_mpdata.o $(BUILD)adv_mpdata.o
	${LINKER} $^ -o $@ ${LFLAGS}

point_in_on_test: $(BUILD)test_point_in_on.o $(OBJS)
	${LINKER} $^ -o $@ ${LFLAGS}

blocking_test: $(BUILD)test_blocking.o $(BUILD)io_routines.o $(BUILD)winds_blocking.o \
				$(BUILD)linear_winds.o $(BUILD)fftshift.o $(BUILD)string.o 		      \
				$(BUILD)lt_lut_io.o $(BUILD)atm_utilities.o $(BUILD)array_utilities.o
	${LINKER} $^ -o $@ ${LFLAGS}

array_util_test: $(BUILD)test_array_utilities.o $(OBJS)
	${LINKER} $^ -o $@ ${LFLAGS}

variable_dict_test: $(BUILD)test_variable_dict.o $(BUILD)co_utilities.o	 	\
					$(BUILD)variable_dict_h.o	$(BUILD)variable_dict_obj.o	\
					$(BUILD)variable_h.o		$(BUILD)variable_obj.o		\
					$(BUILD)grid_h.o			$(BUILD)grid_obj.o			\
					$(BUILD)meta_data_h.o 		$(BUILD)meta_data_obj.o 	\
					$(BUILD)icar_constants.o	$(BUILD)assertions.o
	${LINKER} $^ -o $@ ${LFLAGS}


caf_options_test: $(BUILD)test_caf_options.o $(OBJS)
	${LINKER} $^ -o $@ ${LFLAGS}

caf_write_domain_test: $(BUILD)test_caf_write_domain.o	$(OBJS)
	${LINKER} $^ -o $@ ${LFLAGS}

caf_no_forcing_test: $(BUILD)test_caf_no_forcing.o							\
	 				$(BUILD)domain_h.o		 $(BUILD)domain_obj.o		\
					$(BUILD)boundary_h.o	 $(BUILD)boundary_obj.o		\
					$(BUILD)options_h.o		 $(BUILD)options_obj.o		\
					$(BUILD)opt_types.o      $(BUILD)time_io.o          \
					$(BUILD)grid_h.o		 $(BUILD)grid_obj.o			\
					$(BUILD)output_h.o		 $(BUILD)output_obj.o		\
					$(BUILD)meta_data_h.o    $(BUILD)meta_data_obj.o	\
					$(BUILD)variable_h.o	 $(BUILD)variable_obj.o		\
					$(BUILD)variable_dict_h.o $(BUILD)variable_dict_obj.o	\
					$(BUILD)time_delta_obj.o $(BUILD)time_obj.o    		\
					$(BUILD)time_h.o	\
					$(BUILD)exchangeable_h.o $(BUILD)exchangeable_obj.o \
					$(BUILD)icar_constants.o $(BUILD)io_routines.o		\
					$(BUILD)model_tracking.o $(BUILD)string.o           \
					$(BUILD)assertions.o	 $(BUILD)atm_utilities.o    \ $(BUILD)pbl_utilities.o		\
					$(BUILD)co_utilities.o   $(BUILD)default_output_metadata.o \
					$(BUILD)ra_driver.o      $(BUILD)ra_simple.o		\
					$(BUILD)pbl_driver.o     $(BUILD)pbl_simple.o		$(BUILD)pbl_ysu.o		\
					$(BUILD)mp_driver.o      $(BUILD)mp_simple.o		\
					$(BUILD)mp_thompson.o	 $(BUILD)mp_thompson_aer.o  \
					$(BUILD)mp_wsm6.o		 $(BUILD)mp_wsm3.o 			\
					$(BUILD)cu_driver.o		 $(BUILD)cu_tiedtke.o       \
					$(BUILD)cu_nsas.o		 $(BUILD)cu_bmj.o			\
					$(BUILD)advection_driver.o $(BUILD)geo_reader.o		\
					$(BUILD)advect.o		 $(BUILD)data_structures.o	\
					$(BUILD)vinterp.o		 $(BUILD)array_utilities.o	\
					$(BUILD)lsm_driver.o	\
					$(BUILD)lsm_noahdrv.o 	 $(BUILD)lsm_noahlsm.o \
					$(BUILD)lsm_noahmpdrv.o  $(BUILD)lsm_noahmplsm.o \
					$(BUILD)water_simple.o   $(BUILD)water_lake.o
	${LINKER} $^ -o $@ ${LFLAGS}

caf_domain_init_test: $(BUILD)test_caf_domain_init.o $(OBJS)
	${LINKER} $^ -o $@ ${LFLAGS}

caf_init_test: $(BUILD)test_caf_init.o $(OBJS)
	${LINKER} $^ -o $@ ${LFLAGS}

caf_init-write_test: $(BUILD)test_caf_init-write.o $(OBJS)
	${LINKER} $^ -o $@ ${LFLAGS}

caf_bcast_test: $(BUILD)test_caf_bcast.o	\
				$(BUILD)co_utilities.o $(BUILD)icar_constants.o
	${LINKER} $^ -o $@ ${LFLAGS}

caf_boundary_test: $(BUILD)test_caf_boundary.o $(OBJS)
	${LINKER} $^ -o $@ ${LFLAGS}

caf_boundary_init_test: $(BUILD)test_caf_boundary_init.o $(OBJS)
	${LINKER} $^ -o $@ ${LFLAGS}

caf_linear_winds_setup_test: $(BUILD)test_caf_linear_winds_setup.o $(OBJS)
	${LINKER} $^ -o $@ ${LFLAGS}


caf_other_image_grids_test: $(BUILD)test_caf_other_image_grids.o	\
				$(BUILD)grid_obj.o 			$(BUILD)grid_h.o		\
				$(BUILD)icar_constants.o 	$(BUILD)assertions.o
	${LINKER} $^ -o $@ ${LFLAGS}

caf_one_d_decomposition_test: $(BUILD)test_caf_one_d_decomposition.o
	${LINKER} $^ -o $@ ${LFLAGS}


caf_threads_test: $(BUILD)test_caf_threads.o
	${LINKER} $^ -o $@ ${LFLAGS}

###################################################################
#	Generic compilation rules
###################################################################

$(BUILD)%.o: $(UTIL)%.f90
	${F90} ${FFLAGS} $< -o $@

$(BUILD)%.o: $(CONST)%.f90
	${F90} ${FFLAGS} $< -o $@

$(BUILD)%.o: $(IO)%.f90
	${F90} ${FFLAGS} $< -o $@

$(BUILD)%.o: $(OBJECTS)%.f90
	${F90} ${FFLAGS} $< -o $@

$(BUILD)%.o: $(PHYS)%.f90
	${F90} ${FFLAGS} $< -o $@

$(BUILD)%.o: $(TESTS)%.f90
	${F90} ${FFLAGS} $< -o $@

$(BUILD)%.o: $(MAIN)%.f90
	${F90} ${FFLAGS} $< -o $@


###################################################################
#	driver code for original version of ICAR
###################################################################

$(BUILD)driver.o:$(MAIN)driver.f90 $(BUILD)time_step.o $(BUILD)domain_h.o $(BUILD)timer_h.o \
					$(BUILD)output_h.o $(BUILD)boundary_h.o $(BUILD)init.o $(BUILD)options_h.o \
					$(BUILD)restart_h.o

###################################################################
#	Options Object
###################################################################
$(BUILD)opt_types.o: $(OBJECTS)opt_types.f90 $(BUILD)icar_constants.o $(BUILD)time_h.o $(BUILD)time_delta_obj.o

$(BUILD)options_h.o: $(OBJECTS)options_h.f90 $(BUILD)icar_constants.o $(BUILD)opt_types.o

$(BUILD)options_obj.o: $(OBJECTS)options_obj.f90 $(BUILD)options_h.o $(BUILD)time_h.o		\
						$(BUILD)time_delta_obj.o $(BUILD)io_routines.o $(BUILD)opt_types.o	\
						$(BUILD)time_io.o $(BUILD)string.o $(BUILD)model_tracking.o			\
						$(BUILD)mp_driver.o $(BUILD)advection_driver.o $(BUILD)ra_driver.o	\
						$(BUILD)lsm_driver.o $(BUILD)pbl_driver.o $(BUILD)cu_driver.o 		\
						$(BUILD)wind.o $(BUILD)default_output_metadata.o

###################################################################
#	Primary Data Objects
###################################################################

$(BUILD)meta_data_obj.o: $(OBJECTS)meta_data_obj.f90 $(BUILD)meta_data_h.o

$(BUILD)meta_data_h.o: $(OBJECTS)meta_data_h.f90 $(BUILD)icar_constants.o

$(BUILD)variable_obj.o: $(OBJECTS)variable_obj.f90 $(BUILD)variable_h.o $(BUILD)co_utilities.o

$(BUILD)variable_h.o: $(OBJECTS)variable_h.f90 $(BUILD)meta_data_h.o $(BUILD)icar_constants.o $(BUILD)grid_h.o

$(BUILD)variable_dict_h.o: $(OBJECTS)variable_dict_h.f90 $(BUILD)variable_h.o $(BUILD)icar_constants.o

$(BUILD)variable_dict_obj.o: $(OBJECTS)variable_dict_obj.f90 $(BUILD)variable_dict_h.o

$(BUILD)grid_obj.o: $(OBJECTS)grid_obj.f90 $(BUILD)grid_h.o $(BUILD)assertions.o

$(BUILD)grid_h.o: $(OBJECTS)grid_h.f90 $(BUILD)icar_constants.o

$(BUILD)exchangeable_obj.o: $(OBJECTS)exchangeable_obj.f90 $(BUILD)exchangeable_h.o $(BUILD)assertions.o

$(BUILD)exchangeable_h.o: $(OBJECTS)exchangeable_h.f90 $(BUILD)grid_h.o $(BUILD)variable_h.o $(BUILD)assertions.o

$(BUILD)default_output_metadata.o: $(IO)default_output_metadata.f90 $(BUILD)icar_constants.o $(BUILD)variable_h.o $(BUILD)meta_data_h.o

$(BUILD)restart_h.o: $(IO)restart_h.f90 $(BUILD)domain_h.o $(BUILD)output_h.o $(BUILD)options_h.o

$(BUILD)restart.o: $(IO)restart.f90 $(BUILD)restart_h.o $(BUILD)time_h.o $(BUILD)io_routines.o $(BUILD)icar_constants.o

$(BUILD)output_h.o: $(IO)output_h.f90 $(BUILD)icar_constants.o $(BUILD)variable_h.o $(BUILD)grid_h.o $(BUILD)meta_data_h.o $(BUILD)domain_h.o

$(BUILD)output_obj.o: $(IO)output_obj.f90 $(BUILD)output_h.o $(BUILD)default_output_metadata.o $(BUILD)time_h.o


$(BUILD)boundary_h.o: $(OBJECTS)boundary_h.f90 $(BUILD)options_h.o				\
					$(BUILD)grid_h.o $(BUILD)meta_data_h.o $(BUILD)time_h.o 	\
					$(BUILD)data_structures.o $(BUILD)time_delta_obj.o			\
					$(BUILD)variable_dict_h.o

$(BUILD)boundary_obj.o: $(OBJECTS)boundary_obj.f90 $(BUILD)boundary_h.o					\
					$(BUILD)string.o $(BUILD)atm_utilities.o $(BUILD)time_io.o			\
					$(BUILD)io_routines.o $(BUILD)co_utilities.o $(BUILD)geo_reader.o	\
					$(BUILD)vinterp.o $(BUILD)array_utilities.o $(BUILD)icar_constants.o


$(BUILD)domain_h.o: $(OBJECTS)domain_h.f90 $(BUILD)options_h.o $(BUILD)exchangeable_h.o \
					$(BUILD)grid_h.o 	$(BUILD)meta_data_h.o $(BUILD)data_structures.o	\
					$(BUILD)time_h.o 	$(BUILD)boundary_h.o	 						\
					$(BUILD)time_io.o 	$(BUILD)io_routines.o

$(BUILD)domain_obj.o: $(OBJECTS)domain_obj.f90 $(BUILD)domain_h.o $(BUILD)options_h.o $(BUILD)grid_h.o 	\
					$(BUILD)string.o $(BUILD)assertions.o $(BUILD)geo_reader.o 							\
					$(BUILD)atm_utilities.o $(BUILD)io_routines.o $(BUILD)co_utilities.o 				\
					$(BUILD)array_utilities.o $(BUILD)vinterp.o

###################################################################
#	Core initial and boundary condition and time steping
###################################################################

$(BUILD)init.o:$(MAIN)init.f90 $(BUILD)data_structures.o $(BUILD)options_h.o $(BUILD)domain_h.o $(BUILD)boundary_h.o \
	 				$(BUILD)mp_driver.o $(BUILD)advection_driver.o $(BUILD)atm_utilities.o $(BUILD)wind.o $(BUILD)ra_driver.o \
					$(BUILD)pbl_driver.o $(BUILD)cu_driver.o $(BUILD)lsm_driver.o $(BUILD)io_routines.o

$(BUILD)time_step.o:$(MAIN)time_step.f90 $(BUILD)data_structures.o $(BUILD)wind.o $(BUILD)string.o 	\
					$(BUILD)advection_driver.o $(BUILD)mp_driver.o $(BUILD)domain_h.o 				\
					$(BUILD)options_h.o $(BUILD)atm_utilities.o $(BUILD)ra_driver.o $(BUILD)lsm_driver.o \
					$(BUILD)debug_utils.o $(BUILD)pbl_driver.o $(BUILD)cu_driver.o


###################################################################
#	Utility Routines
###################################################################

$(BUILD)co_utilities.o:$(UTIL)co_utilities.f90 $(BUILD)icar_constants.o

$(BUILD)time_delta_obj.o:$(UTIL)time_delta_obj.f90 $(BUILD)icar_constants.o

$(BUILD)timer_h.o:$(OBJECTS)timer_h.f90

$(BUILD)timer_obj.o:$(UTIL)time_obj.f90 $(BUILD)timer_h.o

$(BUILD)time_h.o:$(UTIL)time_h.f90 $(BUILD)time_delta_obj.o

$(BUILD)time_obj.o:$(UTIL)time_obj.f90 $(BUILD)time_h.o $(BUILD)co_utilities.o

$(BUILD)time_io.o:$(UTIL)time_io.f90 $(BUILD)time_h.o $(BUILD)data_structures.o \
					$(BUILD)icar_constants.o $(BUILD)string.o $(BUILD)io_routines.o

$(BUILD)string.o:$(UTIL)string.f90

$(BUILD)array_utilities.o:$(UTIL)array_utilities.f90

$(BUILD)atm_utilities.o:$(UTIL)atm_utilities.f90 $(BUILD)icar_constants.o 	\
					$(BUILD)data_structures.o $(BUILD)options_h.o $(BUILD)mp_thompson.o \
					$(BUILD)opt_types.o

$(BUILD)pbl_utilities.o:$(UTIL)pbl_utilities.f90 $(BUILD)icar_constants.o 	\
					$(BUILD)data_structures.o


###################################################################
#	I/O routines
###################################################################

$(BUILD)io_routines.o:$(IO)io_routines.f90 $(BUILD)data_structures.o

$(BUILD)lt_lut_io.o: $(IO)lt_lut_io.f90 $(BUILD)data_structures.o $(BUILD)io_routines.o \
				$(BUILD)string.o $(BUILD)opt_types.o $(BUILD)icar_constants.o


###################################################################
#	Interpolation Routines
###################################################################

$(BUILD)geo_reader.o:$(UTIL)geo_reader.f90 $(BUILD)data_structures.o $(BUILD)icar_constants.o

$(BUILD)vinterp.o: $(UTIL)vinterp.f90 $(BUILD)data_structures.o

###################################################################
#	Microphysics code
###################################################################
$(BUILD)mp_driver.o:$(PHYS)mp_driver.f90 $(BUILD)mp_simple.o $(BUILD)mp_thompson_aer.o $(BUILD)mp_thompson.o \
					$(BUILD)mp_wsm6.o $(BUILD)mp_wsm3.o \
					$(BUILD)data_structures.o $(BUILD)domain_h.o  $(BUILD)icar_constants.o\
					$(BUILD)wrf_constants.o

$(BUILD)mp_morrison.o:$(PHYS)mp_morrison.f90 $(BUILD)data_structures.o

$(BUILD)mp_wsm3.o:$(PHYS)mp_wsm3.f90 $(BUILD)wrf_constants.o

$(BUILD)mp_wsm6.o:$(PHYS)mp_wsm6.f90 $(BUILD)wrf_constants.o

$(BUILD)mp_thompson.o:$(PHYS)mp_thompson.f90 $(BUILD)opt_types.o

$(BUILD)mp_thompson_aer.o:$(PHYS)mp_thompson_aer.f90 $(BUILD)co_utilities.o

$(BUILD)mp_simple.o:$(PHYS)mp_simple.f90 $(BUILD)data_structures.o $(BUILD)options_h.o

###################################################################
#	Convection code
###################################################################
$(BUILD)cu_driver.o:$(PHYS)cu_driver.f90 $(BUILD)cu_tiedtke.o $(BUILD)cu_nsas.o $(BUILD)cu_bmj.o  \
					$(BUILD)data_structures.o $(BUILD)icar_constants.o $(BUILD)domain_h.o

$(BUILD)cu_tiedtke.o:$(PHYS)cu_tiedtke.f90

$(BUILD)cu_kf.o:$(PHYS)cu_kf.f90

$(BUILD)cu_nsas.o:$(PHYS)cu_nsas.f90

$(BUILD)cu_bmj.o:$(PHYS)cu_bmj.f90

###################################################################
#	Radiation code
###################################################################

$(BUILD)ra_driver.o:$(PHYS)ra_driver.f90 $(BUILD)ra_simple.o $(BUILD)ra_rrtmg_lw.o $(BUILD)ra_rrtmg_sw.o $(BUILD)data_structures.o $(BUILD)icar_constants.o $(BUILD)options_h.o $(BUILD)domain_h.o

$(BUILD)ra_rrtmg_lw.o:$(PHYS)ra_rrtmg_lw.f90 $(BUILD)ra_clWRF_support.o $(BUILD)io_routines.o

$(BUILD)ra_rrtmg_sw.o:$(PHYS)ra_rrtmg_sw.f90 $(BUILD)ra_clWRF_support.o $(BUILD)io_routines.o $(BUILD)ra_rrtmg_lw.o

$(BUILD)ra_clWRF_support.o:$(PHYS)ra_clWRF_support.f90 $(BUILD)io_routines.o

$(BUILD)ra_simple.o:$(PHYS)ra_simple.f90 $(BUILD)data_structures.o $(BUILD)atm_utilities.o $(BUILD)options_h.o $(BUILD)domain_h.o

###################################################################
#	Surface code
###################################################################
$(BUILD)lsm_driver.o: $(PHYS)lsm_driver.f90 $(BUILD)data_structures.o	\
						$(BUILD)lsm_noahdrv.o $(BUILD)lsm_noahlsm.o		\
						$(BUILD)lsm_noahmpdrv.o $(BUILD)lsm_noahmplsm.o \
						$(BUILD)water_simple.o $(BUILD)water_lake.o		\
						$(BUILD)domain_h.o $(BUILD)options_h.o			\
						$(BUILD)atm_utilities.o $(BUILD)ra_simple.o		\
						$(BUILD)opt_types.o	$(BUILD)icar_constants.o

# $(BUILD)lsm_basic.o $(BUILD)lsm_simple.o

$(BUILD)water_simple.o: $(PHYS)water_simple.f90 $(BUILD)data_structures.o $(BUILD)options_h.o

$(BUILD)water_lake.o: $(PHYS)water_lake.f90 $(BUILD)data_structures.o

$(BUILD)lsm_simple.o: $(PHYS)lsm_simple.f90 $(BUILD)data_structures.o

$(BUILD)lsm_basic.o: $(PHYS)lsm_basic.f90 $(BUILD)data_structures.o

$(BUILD)lsm_noahdrv.o: $(PHYS)lsm_noahdrv.f90 $(BUILD)lsm_noahlsm.o

$(BUILD)lsm_noahlsm.o: $(PHYS)lsm_noahlsm.f90

$(BUILD)lsm_noahmpdrv.o: $(PHYS)lsm_noahmpdrv.f90 $(BUILD)lsm_noahmplsm.o \
	          $(BUILD)lsm_noahmp_gecros.o $(BUILD)lsm_noahmp_glacier.o

$(BUILD)lsm_noahmplsm.o: $(PHYS)lsm_noahmplsm.f90 $(BUILD)lsm_noahmp_gecros.o \
	          $(BUILD)lsm_noahmp_glacier.o

$(BUILD)lsm_noahmp_gecros.o: $(PHYS)lsm_noahmp_gecros.f90

$(BUILD)lsm_noahmp_glacier.o: $(PHYS)lsm_noahmp_glacier.f90

###################################################################
#	Planetary Boundary Layer code
###################################################################
$(BUILD)pbl_driver.o: $(PHYS)pbl_driver.f90 $(BUILD)pbl_simple.o $(BUILD)pbl_ysu.o $(BUILD)domain_h.o \
					$(BUILD)options_h.o $(BUILD)lsm_driver.o $(BUILD)data_structures.o $(BUILD)atm_utilities.o $(BUILD)pbl_utilities.o

$(BUILD)pbl_simple.o: $(PHYS)pbl_simple.f90 $(BUILD)domain_h.o $(BUILD)options_h.o $(BUILD)data_structures.o

$(BUILD)pbl_ysu.o: $(PHYS)pbl_ysu.f90


###################################################################
#	Advection related code
###################################################################
$(BUILD)advection_driver.o:$(PHYS)advection_driver.f90 $(BUILD)data_structures.o $(BUILD)advect.o \
							$(BUILD)options_h.o $(BUILD)domain_h.o $(BUILD)adv_mpdata.o

$(BUILD)advect.o:$(PHYS)advect.f90 $(BUILD)data_structures.o $(BUILD)options_h.o $(BUILD)domain_h.o

$(BUILD)adv_mpdata.o:$(PHYS)adv_mpdata.f90 $(PHYS)adv_mpdata_FCT_core.f90 $(BUILD)data_structures.o $(BUILD)options_h.o $(BUILD)domain_h.o


###################################################################
#	Wind related code
###################################################################
$(BUILD)wind.o:$(PHYS)wind.f90 $(BUILD)data_structures.o $(BUILD)exchangeable_h.o \
				$(BUILD)domain_h.o $(BUILD)options_h.o $(BUILD)linear_winds.o

#########################
# $(BUILD)winds_blocking.o $(BUILD)linear_winds.o
# ${F90} ${FFLAGS} $(PHYS)wind.f90 -o $(BUILD)wind.o
#########################

$(BUILD)linear_winds.o:$(PHYS)linear_winds.f90 $(BUILD)io_routines.o $(BUILD)data_structures.o 	\
	 				   $(BUILD)fftshift.o $(BUILD)string.o $(BUILD)fftw.o						\
					   $(BUILD)atm_utilities.o $(BUILD)array_utilities.o $(BUILD)lt_lut_io.o	\
					   $(BUILD)domain_h.o $(BUILD)options_h.o $(BUILD)grid_h.o

$(BUILD)winds_blocking.o:$(PHYS)winds_blocking.f90 $(BUILD)linear_winds.o 	\
	 					$(BUILD)fftshift.o $(BUILD)fftw.o $(BUILD)array_utilities.o \
						$(BUILD)data_structures.o $(BUILD)atm_utilities.o $(BUILD)string.o

###################################################################
#	FFT code
###################################################################

$(BUILD)fftw.o:$(UTIL)fftw.f90

$(BUILD)fftshift.o:$(UTIL)fftshift.f90 $(BUILD)fftw.o


###################################################################
#	Generic data structures, used by almost everything
###################################################################
$(BUILD)data_structures.o:$(MAIN)data_structures.f90 $(BUILD)icar_constants.o $(BUILD)time_h.o

$(BUILD)wrf_constants.o:$(CONST)wrf_constants.f90

$(BUILD)icar_constants.o:$(CONST)icar_constants.f90

###################################################################
#	Keep track of model versions for user information
###################################################################
$(BUILD)model_tracking.o:$(MAIN)model_tracking.f90

$(BUILD)debug_utils.o:$(UTIL)debug_utils.f90 $(BUILD)domain_h.o $(BUILD)data_structures.o $(BUILD)string.o

###################################################################
#	Unit tests
###################################################################

$(BUILD)test_caf_other_image_grids.o: $(TESTS)test_caf_other_image_grids.f90 $(BUILD)grid_h.o

$(BUILD)test_caf_domain_init.o: $(TESTS)test_caf_domain_init.f90 $(BUILD)domain_h.o $(BUILD)options_h.o

$(BUILD)test_caf_init.o: $(TESTS)test_caf_init.f90 $(BUILD)domain_h.o $(BUILD)options_h.o $(BUILD)boundary_h.o

$(BUILD)test_caf_init-write.o: $(TESTS)test_caf_init-write.f90 $(BUILD)domain_h.o $(BUILD)options_h.o $(BUILD)boundary_h.o $(BUILD)output_h.o

$(BUILD)test_caf_boundary_init.o: $(TESTS)test_caf_boundary_init.f90 $(BUILD)boundary_h.o $(BUILD)options_h.o $(BUILD)variable_h.o

$(BUILD)test_caf_linear_winds_setup.o: $(TESTS)test_caf_linear_winds_setup.f90 $(BUILD)domain_h.o $(BUILD)options_h.o $(BUILD)linear_winds.o

$(BUILD)test_caf_boundary.o: $(TESTS)test_caf_boundary.f90 $(BUILD)domain_h.o $(BUILD)boundary_h.o

$(BUILD)test_caf_write_domain.o: $(TESTS)test_caf_write_domain.f90 $(BUILD)domain_h.o $(BUILD)options_h.o \
							$(BUILD)grid_h.o $(BUILD)output_h.o $(BUILD)atm_utilities.o

$(BUILD)test_caf_no_forcing.o: $(TESTS)test_caf_no_forcing.f90 $(BUILD)domain_h.o $(BUILD)options_h.o \
							$(BUILD)grid_h.o $(BUILD)output_h.o $(BUILD)atm_utilities.o

$(BUILD)test_caf_output.o: $(TESTS)test_caf_output.f90 $(BUILD)output_h.o $(BUILD)variable_h.o  \
							$(BUILD)default_output_metadata.o $(BUILD)icar_constants.o

$(BUILD)test_caf_bcast.o: $(TESTS)test_caf_bcast.f90 $(BUILD)co_utilities.o

$(BUILD)test_variable_dict.o: $(TESTS)test_variable_dict.f90 $(BUILD)grid_h.o \
						$(BUILD)variable_h.o $(BUILD)variable_dict_h.o

$(BUILD)test_caf_options.o: $(TESTS)test_caf_options.f90 $(BUILD)options_h.o

$(BUILD)test_caf_threads.o: $(TESTS)test_caf_threads.f90

$(BUILD)test_fftw.o: $(TESTS)test_fftw.f90 $(BUILD)fftw.o

$(BUILD)test_fftshift.o:$(TESTS)test_fftshift.f90 $(BUILD)fftshift.o

$(BUILD)test_calendar.o:$(TESTS)test_calendar.f90  $(BUILD)time_obj.o $(BUILD)time_h.o $(BUILD)time_delta_obj.o

$(BUILD)test_point_in_on.o:$(TESTS)test_point_in_on.f90 $(BUILD)geo_reader.o

$(BUILD)test_mpdata.o:$(TESTS)test_mpdata.f90 $(BUILD)adv_mpdata.o

$(BUILD)test_blocking.o:$(TESTS)test_blocking.f90 $(BUILD)winds_blocking.o $(BUILD)linear_winds.o \
						$(BUILD)data_structures.o $(BUILD)icar_constants.o $(BUILD)io_routines.o

$(BUILD)test_array_utilities.o:$(TESTS)test_array_utilities.f90 $(BUILD)array_utilities.o
