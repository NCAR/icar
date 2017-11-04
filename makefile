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
# LIBFFT	 : location of fftw libraries		default = /usr/local/lib
# INCFFT	 : location of fftw headers			default = /usr/local/include
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
RM=/bin/rm
CP=/bin/cp
# doxygen only required to enable "make doc"
DOXYGEN=doxygen

ifndef FC
	FC=gfortran
endif
F90=${FC}

ifndef FFTW
	FFTW=/usr/local
endif
FFTW_PATH = ${FFTW}
LIBFFT = -L${FFTW_PATH}/lib -lm -lfftw3
INCFFT = -I${FFTW_PATH}/include

ifndef NETCDF
	NETCDF=/usr/local
endif
NCDF_PATH = ${NETCDF}
LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
INCNETCDF = -I$(NCDF_PATH)/include

########################################################################################
# Try to find the machine information
########################################################################################
NODENAME := $(shell uname -n)
ifeq ($(NODENAME), Patthar.local)
	NODENAME=Nomad.local
endif
ifeq ($(patsubst vpn%.ucar.edu,vpn.ucar.edu,$(NODENAME)), vpn.ucar.edu)
	NODENAME=Nomad.local
endif
ifeq ($(patsubst hexagon%,hexagon,$(NODENAME)), hexagon)
	NODENAME=hexagon
endif
ifeq ($(patsubst cheyenne%,cheyenne,$(NODENAME)), cheyenne)
	NODENAME=cheyenne
endif


# on hexagon (uib computer)
ifeq ($(NODENAME), hexagon)
	F90=ftn
	FFTW_PATH = /home/gfi/pbo003/Libraries/FFTW/fftw-3.3.4
	LIBFFT = -L${FFTW_PATH}/lib -lm -lfftw3
	INCFFT = -I${FFTW_PATH}/include
	NCDF_PATH = /opt/cray/netcdf/default/cray/83
	LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
	INCNETCDF = -I$(NCDF_PATH)/include
endif

# traveling laptop / home computer
# may want to add -Wl,-no_compact_unwind to suppress some warnings
ifeq ($(NODENAME), Nomad.local)
	F90=gfortran
	LIBFFT=-L/Users/gutmann/usr/local/lib -lm -lfftw3
	INCFFT=-I/Users/gutmann/usr/local/include
	NCDF_PATH = /usr/local/
	LIBNETCDF = -L/usr/local/gfortran/lib -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
	INCNETCDF = -I$(NCDF_PATH)/include
endif
ifeq ($(NODENAME), dablam.rap.ucar.edu)
	F90=gfortran
	LIBFFT=-L/usr/local/lib -lm -lfftw3
	INCFFT=-I/usr/local/include
	NCDF_PATH = /usr/local
	LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
	INCNETCDF = -I$(NCDF_PATH)/include
endif
# hydro-c1 cluster
ifeq ($(NODENAME),hydro-c1)
	F90=ifort
	NCDF_PATH = /opt/netcdf4-intel

	# F90=gfortran
	# NCDF_PATH = /opt/netcdf4-gcc

	LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
	INCNETCDF = -I$(NCDF_PATH)/include

	LIBFFT=-L/home/gutmann/.usr/local/lib -lm -lfftw3
	INCFFT=-I/home/gutmann/.usr/local/include
endif
# on yellowstone:
ifeq ($(LMOD_FAMILY_COMPILER),gnu)
	F90=gfortran
	LIBFFT = -L/glade/u/home/gutmann/usr/local/lib -lm -lfftw3
	INCFFT = -I/glade/u/home/gutmann/usr/local/include
	NCDF_PATH=/glade/apps/opt/netcdf/4.3.0/gnu/4.8.2
	# note this works for almost all versions of gfortran EXCEPT 4.8.2... the default version
	NCDF_PATH=/glade/apps/opt/netcdf/4.3.3.1/gnu/$(GNU_MAJOR_VERSION).$(GNU_MINOR_VERSION)
	# LIBNETCDF = $(LIB_NCAR) # when netcdf includes are setup by the yellowstone module system
	# INCNETCDF = $(INC_NCAR)
	LIBNETCDF = -Wl,-rpath,$(NCDF_PATH)/lib -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf # if using a compiler for which netcdf includes are
	INCNETCDF = -I$(NCDF_PATH)/include # NOT setup correctly by the yellowstone module system
endif
ifeq ($(LMOD_FAMILY_COMPILER),intel)
	F90=ifort
	LIBFFT = -L/glade/u/home/gutmann/usr/local/lib -lm -lfftw3
	INCFFT = -I/glade/u/home/gutmann/usr/local/include
	NCDF_PATH=/glade/apps/opt/netcdf/4.3.0/intel/default
	LIBNETCDF = $(LIB_NCAR) #-L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
	INCNETCDF = $(INC_NCAR) #-I$(NCDF_PATH)/include # netcdf includes are setup by the yellowstone module system
endif
ifeq ($(LMOD_FAMILY_COMPILER),pgi)
	F90=pgf90
	LIBFFT = -L/glade/u/home/gutmann/usr/local/lib -lm -lfftw3
	INCFFT = -I/glade/u/home/gutmann/usr/local/include
	NCDF_PATH=/glade/apps/opt/netcdf/4.3.0/pgi/default
	LIBNETCDF = -rpath $(NCDF_PATH)/lib -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf # if using a compiler for which netcdf includes are
	INCNETCDF = -I$(NCDF_PATH)/include # NOT setup correctly by the yellowstone module system
endif

ifeq ($(NODENAME), cheyenne)
	F90=$(FC)
	FFTW_PATH = /glade/u/home/gutmann/usr/local
	LIBFFT = -L$(FFTW_PATH)/lib -lm -lfftw3
	INCFFT = -I$(FFTW_PATH)/include
	NCDF_PATH = $(NETCDF)
	LIBNETCDF = -L$(NETCDF)/lib -lnetcdff -lnetcdf
	INCNETCDF = -I$(NETCDF)/include
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
ifeq ($(F90), gfortran)
	COMP=-fopenmp -lgomp -O3 -c -ffree-line-length-none -ftree-vectorize -fimplicit-none -funroll-loops -march=native  -fno-protect-parens # -ffast-math #-flto #
	LINK=-fopenmp -lgomp
	PREPROC=-cpp
	MODOUTPUT=-J $(BUILD)
endif
# Intel fortran
ifeq ($(F90), ifort)
	COMP=-c -u -qopenmp -liomp5 -O3 -xHost -ftz -fpe0 # -check stack,bounds -fp-stack-check
	LINK= -qopenmp -liomp5
	PREPROC=-fpp
	MODOUTPUT=-module $(BUILD)
endif
# PGI fortran
ifeq ($(F90), pgf90)
	COMP=-O2 -mp -c -Mdclchk #-fast -O3 -mp -c -Mdclchk
	LINK=-mp
	PREPROC=-Mpreprocess
	MODOUTPUT=-module $(BUILD)
endif

# Cray fortran
ifeq ($(F90), ftn)
	COMP= -h omp vector2 -O2 -c -eI
	LINK= -fopenmp
	PREPROC= -eZ
	MODOUTPUT= -J $(BUILD) -em
endif


# Various compiling options.  Set the MODE variable with "make MODE=debugslow" etc.
ifeq ($(MODE), debugslow)
	ifeq ($(F90), ifort)
		COMP= -debug -debug-parameters all -traceback -ftrapuv -g -fpe0 -c -u -check all -check noarg_temp_created -CB
		LINK=
	endif
	ifeq ($(F90), gfortran)
		COMP= -c -g -fbounds-check -fbacktrace -finit-real=nan -ffree-line-length-none
		LINK=
	endif
	ifeq ($(F90), pgf90)
		COMP= -c -g -Mbounds -Mlist -Minfo  -Mdclchk
		LINK=
	endif
	ifeq ($(F90), ftn)
		COMP=-h noomp -c -g -m 0 -R abcsp
		LINK=-h noomp
		PREPROC=-eZ
		MODOUTPUT=-e m -J $(BUILD)
	endif
endif
ifeq ($(MODE), debug)
	ifeq ($(F90), ifort)
		COMP= -debug -c -O1 -u -check all -check noarg_temp_created -traceback -fpe0 -fast-transcendentals -xhost
		LINK=
	endif
	ifeq ($(F90), gfortran)
		COMP= -c -O1 -g -fbounds-check -fbacktrace -finit-real=nan -ffree-line-length-none
		LINK=
	endif
	ifeq ($(F90), pgf90)
		COMP= -c -gopt -O1 -Mbounds -Mlist -Minfo  -Mdclchk
		LINK=
	endif
	ifeq ($(F90), ftn)
		COMP=-O1 -h noomp -c -g
		LINK=-h noomp
		PREPROC=-eZ
		MODOUTPUT=-e m -J $(BUILD)
	endif
endif
ifeq ($(MODE), debugompslow)
	ifeq ($(F90), ifort)
		# COMP= -openmp -liomp5 -debug -debug-parameters all -traceback -ftrapuv -g -fpe0 -c -u -check all -check noarg_temp_created -CB
		COMP= -qopenmp -liomp5 -debug -c -u	-fpe0 -traceback -check all -check noarg_temp_created -fp-stack-check
		LINK= -qopenmp -liomp5
	endif
	ifeq ($(F90), gfortran)
		COMP= -fopenmp -lgomp -c -g -fbounds-check -fbacktrace -finit-real=nan -ffree-line-length-none
		LINK= -fopenmp -lgomp
	endif
	ifeq ($(F90), pgf90)
		COMP= -c -g -Mbounds -Mlist -Minfo -mp -Mdclchk
		LINK= -mp
	endif
	ifeq ($(F90), ftn)
		COMP=-c -g -m 0 -R abcsp
		LINK=
		PREPROC=-eZ
		MODOUTPUT=-e m -J $(BUILD)
	endif
endif
ifeq ($(MODE), debugomp)
	ifeq ($(F90), ifort)
		COMP= -qopenmp -liomp5 -debug -c -O3 -u -traceback -fpe0 -ftz -xHost # -fast-transcendentals -check all -check noarg_temp_created -fpe0
		LINK= -qopenmp -liomp5
	endif
	ifeq ($(F90), gfortran)
		COMP= -fopenmp -lgomp -c -O1 -g -fbounds-check -fbacktrace -finit-real=nan -ffree-line-length-none
		LINK= -fopenmp -lgomp
	endif
	ifeq ($(F90), pgf90)
		COMP= -c -g -O1 -Mbounds -Mlist -Minfo -mp -Mdclchk
		LINK= -mp
	endif
	ifeq ($(F90), ftn)
		COMP=-O1 -c -g
		LINK=
		PREPROC=-eZ
		MODOUTPUT=-e m -J $(BUILD)
	endif
endif

PROF=
ifeq ($(MODE), profile)
	ifeq ($(F90), ifort)
		PROF=-pg -debug inline-debug-info -shared-intel
		COMP=-c -u -qopenmp -liomp5 -O3 -xHost -ftz #because -fast includes -static # not available in ifort <13 -align array64byte
	endif
	ifeq ($(F90), gfortran)
		PROF=-g
	endif
endif
ifeq ($(MODE), fast) # WARNING -ipo (included in -fast) takes forever for very little gain, and this may be unstable
	ifeq ($(F90), ifort)
		COMP=-c -u -openmp -liomp5 -fast -ftz #-fast-transcendentals # not available in ifort <13: -align array64byte
	endif
endif
###################################################################
###################################################################
#
# Should not need to edit anything below this line
#
###################################################################
###################################################################
# copy required libraries into a directory accessible on compute nodes and set LD_RUN_PATH e.g.
# export LD_RUN_PATH=$LD_RUN_PATH:/path/to/netcdf/libraries/lib:/path/to/fftw/libraries/lib
LFLAGS=$(LINK) $(PROF) ${LIBNETCDF} ${LIBFFT}
FFLAGS=$(COMP) $(PROF) ${INCNETCDF} ${INCFFT} ${MODOUTPUT}


$(info $$NODENAME    = ${NODENAME})
$(info $$FC          = ${F90})
$(info $$FFTW_PATH   = ${FFTW_PATH})
$(info $$NCDF_PATH   = ${NCDF_PATH})
$(info $$GIT_VERSION = ${GIT_VERSION})
$(info $$COMP        = ${COMP})
$(info $$LINK        = ${LINK})
$(info $$MODE        = ${MODE})


# Model directories
BUILD=build/
PHYS=physics/
IO=io/
MAIN=main/
UTIL=utilities/
CONST=constants/

OBJS=	$(BUILD)driver.o 			\
		$(BUILD)init.o 				\
		$(BUILD)init_options.o 		\
		$(BUILD)model_tracking.o 	\
		$(BUILD)boundary.o 			\
		$(BUILD)time_step.o 		\
		$(BUILD)output.o 			\
		$(BUILD)io_routines.o 		\
		$(BUILD)lt_lut_io.o 		\
		$(BUILD)mp_driver.o 		\
		$(BUILD)mp_wsm6.o 			\
		$(BUILD)mp_thompson.o 		\
		$(BUILD)mp_simple.o 		\
		$(BUILD)mp_morrison.o 		\
		$(BUILD)cu_driver.o 		\
		$(BUILD)cu_tiedtke.o 		\
		$(BUILD)cu_kf.o 			\
		$(BUILD)ra_driver.o 		\
		$(BUILD)ra_simple.o 		\
		$(BUILD)lsm_driver.o 		\
		$(BUILD)lsm_simple.o 		\
		$(BUILD)lsm_basic.o 		\
		$(BUILD)lsm_noahdrv.o 		\
		$(BUILD)lsm_noahlsm.o 		\
		$(BUILD)water_simple.o 		\
		$(BUILD)pbl_driver.o 		\
		$(BUILD)pbl_simple.o 		\
		$(BUILD)pbl_ysu.o 			\
		$(BUILD)advection_driver.o 	\
		$(BUILD)adv_mpdata.o 		\
		$(BUILD)advect.o 			\
		$(BUILD)wind.o 				\
		$(BUILD)linear_winds.o 		\
		$(BUILD)winds_blocking.o	\
		$(BUILD)fftw.o 				\
		$(BUILD)fftshift.o 			\
		$(BUILD)geo_reader.o 		\
		$(BUILD)vinterp.o 			\
		$(BUILD)atm_utilities.o 	\
		$(BUILD)time.o 				\
		$(BUILD)data_structures.o 	\
		$(BUILD)icar_constants.o	\
		$(BUILD)wrf_constants.o 	\
		$(BUILD)string.o 			\
		$(BUILD)array_utilities.o 	\
		$(BUILD)debug_utils.o


TEST_EXECUTABLES= 	fftshift_test	\
	  				calendar_test	\
	  				mpdata_test		\
	  				fftw_test		\
	  				blocking_test	\
					point_in_on_test\
	  				array_util_test


###################################################################
#	User facing rules
###################################################################

icar:${OBJS}
	${F90} -o icar ${OBJS} ${LFLAGS}

all:icar test

install:icar
	${CP} icar ${INSTALLDIR}

clean:
	${RM} $(BUILD)*.o $(BUILD)*.mod *.lst docs/doxygen_sqlite3.db 2>/dev/null ||:

allclean:cleanall

cleanall: clean
	${RM} icar $(TEST_EXECUTABLES) 2>/dev/null ||:

test: $(TEST_EXECUTABLES)

doc:
	doxygen docs/doxygenConfig

###################################################################
#	test cases
###################################################################
fftw_test: $(BUILD)test_fftw.o
	${F90} $^ -o $@ ${LFLAGS}

fftshift_test: $(BUILD)test_fftshift.o $(BUILD)fftshift.o
	${F90} $^ -o $@ ${LFLAGS}

calendar_test: $(BUILD)test_calendar.o $(BUILD)time.o
	${F90} $^ -o $@ ${LFLAGS}

mpdata_test: $(BUILD)test_mpdata.o $(BUILD)adv_mpdata.o
	${F90} $^ -o $@ ${LFLAGS}

point_in_on_test: $(BUILD)test_point_in_on.o $(BUILD)geo_reader.o $(BUILD)data_structures.o
	${F90} $^ -o $@ ${LFLAGS}

blocking_test: $(BUILD)test_blocking.o $(BUILD)io_routines.o $(BUILD)winds_blocking.o \
				$(BUILD)linear_winds.o $(BUILD)fftshift.o $(BUILD)string.o 		      \
				$(BUILD)lt_lut_io.o $(BUILD)atm_utilities.o $(BUILD)array_utilities.o
	${F90} $^ -o $@ ${LFLAGS}

array_util_test: $(BUILD)test_array_utilities.o $(BUILD)array_utilities.o
	${F90} $^ -o $@ ${LFLAGS}



###################################################################
#	driver code for
###################################################################

$(BUILD)driver.o:$(MAIN)driver.f90 $(BUILD)data_structures.o $(BUILD)init.o $(BUILD)time_step.o \
					$(BUILD)output.o $(BUILD)boundary.o $(BUILD)time.o $(BUILD)string.o
	${F90} ${FFLAGS} $(MAIN)driver.f90 -o $(BUILD)driver.o


###################################################################
#	Core initial and boundary condition and time steping
###################################################################

$(BUILD)init.o:$(MAIN)init.f90 $(BUILD)data_structures.o $(BUILD)io_routines.o $(BUILD)geo_reader.o $(BUILD)vinterp.o \
					$(BUILD)mp_driver.o $(BUILD)cu_driver.o $(BUILD)pbl_driver.o $(BUILD)wind.o \
					$(BUILD)ra_driver.o $(BUILD)lsm_driver.o $(BUILD)init_options.o $(BUILD)advection_driver.o \
					$(BUILD)atm_utilities.o
	${F90} ${FFLAGS} $(MAIN)init.f90 -o $(BUILD)init.o

$(BUILD)boundary.o:$(MAIN)boundary.f90 $(BUILD)data_structures.o $(BUILD)io_routines.o $(BUILD)wind.o $(BUILD)geo_reader.o \
					$(BUILD)vinterp.o $(BUILD)output.o $(BUILD)linear_winds.o
	${F90} ${FFLAGS} $(MAIN)boundary.f90 -o $(BUILD)boundary.o

$(BUILD)time_step.o:$(MAIN)time_step.f90 $(BUILD)data_structures.o $(BUILD)wind.o $(BUILD)output.o \
					$(BUILD)advection_driver.o $(BUILD)ra_driver.o $(BUILD)lsm_driver.o $(BUILD)cu_driver.o \
					$(BUILD)pbl_driver.o $(BUILD)mp_driver.o $(BUILD)boundary.o $(BUILD)debug_utils.o
	${F90} ${FFLAGS} $(MAIN)time_step.f90 -o $(BUILD)time_step.o

$(BUILD)init_options.o:$(MAIN)init_options.f90 $(BUILD)data_structures.o  $(BUILD)io_routines.o \
					$(BUILD)model_tracking.o $(BUILD)time.o
	${F90} ${FFLAGS} $(MAIN)init_options.f90 -o $(BUILD)init_options.o


###################################################################
#	Utility Routines
###################################################################

$(BUILD)time.o:$(UTIL)time.f90
	${F90} ${FFLAGS} $(UTIL)time.f90 -o $(BUILD)time.o

$(BUILD)string.o:$(UTIL)string.f90
	${F90} ${FFLAGS} $(UTIL)string.f90 -o $(BUILD)string.o

$(BUILD)array_utilities.o:$(UTIL)array_utilities.f90
	${F90} ${FFLAGS} $(UTIL)array_utilities.f90 -o $(BUILD)array_utilities.o

$(BUILD)atm_utilities.o:$(UTIL)atm_utilities.f90 $(BUILD)icar_constants.o $(BUILD)data_structures.o
	${F90} ${FFLAGS} $(UTIL)atm_utilities.f90 -o $(BUILD)atm_utilities.o


###################################################################
#	I/O routines
###################################################################

$(BUILD)output.o:$(IO)output.f90 $(BUILD)data_structures.o $(BUILD)io_routines.o $(BUILD)time.o $(BUILD)string.o
	${F90} ${FFLAGS} ${PREPROC} -DVERSION=\"$(GIT_VERSION)\" $(IO)output.f90 -o $(BUILD)output.o

$(BUILD)io_routines.o:$(IO)io_routines.f90 $(BUILD)data_structures.o
	${F90} ${FFLAGS} $(IO)io_routines.f90 -o $(BUILD)io_routines.o

$(BUILD)lt_lut_io.o: $(IO)lt_lut_io.f90 $(BUILD)data_structures.o $(BUILD)io_routines.o $(BUILD)string.o
	${F90} ${FFLAGS} $(IO)lt_lut_io.f90 -o $(BUILD)lt_lut_io.o


###################################################################
#	Interpolation Routines
###################################################################

$(BUILD)geo_reader.o:$(UTIL)geo_reader.f90 $(BUILD)data_structures.o
	${F90} ${FFLAGS} $(UTIL)geo_reader.f90 -o $(BUILD)geo_reader.o

$(BUILD)vinterp.o: $(UTIL)vinterp.f90 $(BUILD)data_structures.o
	${F90} ${FFLAGS} $(UTIL)vinterp.f90 -o $(BUILD)vinterp.o

###################################################################
#	Microphysics code
###################################################################

$(BUILD)mp_driver.o:$(PHYS)mp_driver.f90 $(BUILD)mp_thompson.o $(BUILD)mp_simple.o \
					$(BUILD)mp_morrison.o $(BUILD)data_structures.o $(BUILD)mp_wsm6.o $(BUILD)time.o
	${F90} ${FFLAGS} $(PHYS)mp_driver.f90 -o $(BUILD)mp_driver.o

$(BUILD)mp_morrison.o:$(PHYS)mp_morrison.f90 $(BUILD)data_structures.o
	${F90} ${FFLAGS} $(PHYS)mp_morrison.f90 -o $(BUILD)mp_morrison.o

$(BUILD)mp_wsm6.o:$(PHYS)mp_wsm6.f90 $(BUILD)wrf_constants.o
	${F90} ${FFLAGS} $(PHYS)mp_wsm6.f90 -o $(BUILD)mp_wsm6.o

$(BUILD)mp_thompson.o:$(PHYS)mp_thompson.f90 $(BUILD)data_structures.o
	${F90} ${FFLAGS} $(PHYS)mp_thompson.f90 -o $(BUILD)mp_thompson.o

$(BUILD)mp_simple.o:$(PHYS)mp_simple.f90 $(BUILD)data_structures.o
	${F90} ${FFLAGS} $(PHYS)mp_simple.f90 -o $(BUILD)mp_simple.o

###################################################################
#	Convection code
###################################################################
$(BUILD)cu_driver.o:$(PHYS)cu_driver.f90 $(BUILD)cu_tiedtke.o $(BUILD)cu_kf.o $(BUILD)data_structures.o
	${F90} ${FFLAGS} $(PHYS)cu_driver.f90 -o $(BUILD)cu_driver.o

$(BUILD)cu_tiedtke.o:$(PHYS)cu_tiedtke.f90
	${F90} ${FFLAGS} $(PHYS)cu_tiedtke.f90 -o $(BUILD)cu_tiedtke.o

$(BUILD)cu_kf.o:$(PHYS)cu_kf.f90
	${F90} ${FFLAGS} $(PHYS)cu_kf.f90 -o $(BUILD)cu_kf.o

###################################################################
#	Radiation code
###################################################################

$(BUILD)ra_driver.o:$(PHYS)ra_driver.f90 $(BUILD)ra_simple.o $(BUILD)data_structures.o
	${F90} ${FFLAGS} $(PHYS)ra_driver.f90 -o $(BUILD)ra_driver.o

$(BUILD)ra_simple.o:$(PHYS)ra_simple.f90 $(BUILD)data_structures.o $(BUILD)time.o
	${F90} ${FFLAGS} $(PHYS)ra_simple.f90 -o $(BUILD)ra_simple.o

###################################################################
#	Surface code
###################################################################
$(BUILD)lsm_driver.o: $(PHYS)lsm_driver.f90 $(BUILD)data_structures.o $(BUILD)lsm_simple.o \
						$(BUILD)lsm_basic.o $(BUILD)lsm_noahdrv.o $(BUILD)lsm_noahlsm.o \
						$(BUILD)water_simple.o
	${F90} ${FFLAGS} $(PHYS)lsm_driver.f90 -o $(BUILD)lsm_driver.o

$(BUILD)water_simple.o: $(PHYS)water_simple.f90 $(BUILD)data_structures.o
	${F90} ${FFLAGS} $(PHYS)water_simple.f90 -o $(BUILD)water_simple.o

$(BUILD)lsm_simple.o: $(PHYS)lsm_simple.f90 $(BUILD)data_structures.o
	${F90} ${FFLAGS} $(PHYS)lsm_simple.f90 -o $(BUILD)lsm_simple.o

$(BUILD)lsm_basic.o: $(PHYS)lsm_basic.f90 $(BUILD)data_structures.o
	${F90} ${FFLAGS} $(PHYS)lsm_basic.f90 -o $(BUILD)lsm_basic.o

$(BUILD)lsm_noahdrv.o: $(PHYS)lsm_noahdrv.f90 $(BUILD)lsm_noahlsm.o
	${F90} ${FFLAGS} $(PHYS)lsm_noahdrv.f90 -o $(BUILD)lsm_noahdrv.o

$(BUILD)lsm_noahlsm.o: $(PHYS)lsm_noahlsm.f90
	${F90} ${FFLAGS} $(PHYS)lsm_noahlsm.f90 -o $(BUILD)lsm_noahlsm.o


###################################################################
#	Planetary Boundary Layer code
###################################################################
$(BUILD)pbl_driver.o: $(PHYS)pbl_driver.f90 $(BUILD)pbl_simple.o $(BUILD)pbl_ysu.o $(BUILD)data_structures.o
	${F90} ${FFLAGS} $(PHYS)pbl_driver.f90 -o $(BUILD)pbl_driver.o

$(BUILD)pbl_simple.o: $(PHYS)pbl_simple.f90 $(BUILD)data_structures.o
	${F90} ${FFLAGS} $(PHYS)pbl_simple.f90 -o $(BUILD)pbl_simple.o

$(BUILD)pbl_ysu.o: $(PHYS)pbl_ysu.f90
	${F90} ${FFLAGS} $(PHYS)pbl_ysu.f90 -o $(BUILD)pbl_ysu.o


###################################################################
#	Advection related code
###################################################################
$(BUILD)advection_driver.o:$(PHYS)advection_driver.f90 $(BUILD)data_structures.o $(BUILD)advect.o $(BUILD)adv_mpdata.o \
							$(BUILD)debug_utils.o
	${F90} ${FFLAGS} $(PHYS)advection_driver.f90 -o $(BUILD)advection_driver.o

$(BUILD)advect.o:$(PHYS)advect.f90 $(BUILD)data_structures.o
	${F90} ${FFLAGS} $(PHYS)advect.f90 -o $(BUILD)advect.o

$(BUILD)adv_mpdata.o:$(PHYS)adv_mpdata.f90 $(PHYS)adv_mpdata_FCT_core.f90 $(BUILD)data_structures.o
	${F90} ${FFLAGS} $(PHYS)adv_mpdata.f90 -o $(BUILD)adv_mpdata.o


###################################################################
#	Wind related code
###################################################################
$(BUILD)wind.o:$(PHYS)wind.f90 $(BUILD)data_structures.o \
				$(BUILD)winds_blocking.o $(BUILD)linear_winds.o
	${F90} ${FFLAGS} $(PHYS)wind.f90 -o $(BUILD)wind.o

$(BUILD)linear_winds.o:$(PHYS)linear_winds.f90 $(BUILD)io_routines.o $(BUILD)data_structures.o \
	 				   $(BUILD)fftshift.o $(BUILD)lt_lut_io.o $(BUILD)string.o $(BUILD)fftw.o  \
					   $(BUILD)atm_utilities.o $(BUILD)array_utilities.o
	${F90} ${FFLAGS} $(PHYS)linear_winds.f90 -o $(BUILD)linear_winds.o

$(BUILD)winds_blocking.o:$(PHYS)winds_blocking.f90 $(BUILD)linear_winds.o 	\
	 					$(BUILD)fftshift.o $(BUILD)fftw.o $(BUILD)array_utilities.o \
						$(BUILD)data_structures.o $(BUILD)atm_utilities.o $(BUILD)string.o
	${F90} ${FFLAGS} $(PHYS)winds_blocking.f90 -o $(BUILD)winds_blocking.o

###################################################################
#	FFT code
###################################################################

$(BUILD)fftw.o:$(UTIL)fftw.f90
	${F90} ${FFLAGS} $(UTIL)fftw.f90 -o $(BUILD)fftw.o

$(BUILD)fftshift.o:$(UTIL)fftshift.f90 $(BUILD)fftw.o
	${F90} ${FFLAGS} $(UTIL)fftshift.f90 -o $(BUILD)fftshift.o


###################################################################
#	Generic data structures, used by almost everything
###################################################################
$(BUILD)data_structures.o:$(MAIN)data_structures.f90 $(BUILD)icar_constants.o
	${F90} ${FFLAGS} $(MAIN)data_structures.f90 -o $(BUILD)data_structures.o

$(BUILD)wrf_constants.o:$(CONST)wrf_constants.f90
	${F90} ${FFLAGS} $(CONST)wrf_constants.f90 -o $(BUILD)wrf_constants.o

$(BUILD)icar_constants.o:$(CONST)icar_constants.f90
	${F90} ${FFLAGS} $(CONST)icar_constants.f90 -o $(BUILD)icar_constants.o

###################################################################
#	Keep track of model versions for user information
###################################################################
$(BUILD)model_tracking.o:$(MAIN)model_tracking.f90
	${F90} ${FFLAGS} $(MAIN)model_tracking.f90 -o $(BUILD)model_tracking.o

$(BUILD)debug_utils.o:$(UTIL)debug_utils.f90 $(BUILD)data_structures.o $(BUILD)string.o
	${F90} ${FFLAGS} $(UTIL)debug_utils.f90 -o $(BUILD)debug_utils.o

###################################################################
#	Unit tests
###################################################################
$(BUILD)test_fftw.o: tests/test_fftw.f90 $(BUILD)fftw.o
	${F90} ${FFLAGS} tests/test_fftw.f90 -o $(BUILD)test_fftw.o

$(BUILD)test_fftshift.o:$(BUILD)fftshift.o tests/test_fftshift.f90
	${F90} ${FFLAGS} tests/test_fftshift.f90 -o $(BUILD)test_fftshift.o

$(BUILD)test_calendar.o:$(BUILD)time.o tests/test_calendar.f90
	${F90} ${FFLAGS} tests/test_calendar.f90 -o $(BUILD)test_calendar.o

$(BUILD)test_point_in_on.o:$(BUILD)geo_reader.o tests/test_point_in_on.f90
	${F90} ${FFLAGS} tests/test_point_in_on.f90 -o $(BUILD)test_point_in_on.o

$(BUILD)test_mpdata.o:$(BUILD)adv_mpdata.o tests/test_mpdata.f90
	${F90} ${FFLAGS} tests/test_mpdata.f90 -o $(BUILD)test_mpdata.o

$(BUILD)test_blocking.o:tests/test_blocking.f90 $(BUILD)winds_blocking.o $(BUILD)linear_winds.o \
						$(BUILD)data_structures.o $(BUILD)icar_constants.o $(BUILD)io_routines.o
	${F90} ${FFLAGS} tests/test_blocking.f90 -o $(BUILD)test_blocking.o

$(BUILD)test_array_utilities.o:tests/test_array_utilities.f90 $(BUILD)array_utilities.o
	${F90} ${FFLAGS} tests/test_array_utilities.f90 -o $(BUILD)test_array_utilities.o
