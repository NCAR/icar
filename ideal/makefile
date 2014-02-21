###################################################################
# Makefile rules:
# <default>	all: makes real
# 			install: makes and installs real in INSTALLDIR
# 			clean: removes objects and module files
#		 	allclean: makes clean and removes executables
#		 	cleanall: alias for allclean
#		 	tests: makes various unit tests (does not always work)
# 			real: makes the primary model
#		 	ideal: makes an idealized version of the model (does not not always work)
# 
###################################################################
# Variables that need to be set by the user: 
# 
# INSTALLDIR : default = ~/bin/ 
# LIBFFT     : location of fftw libraries 		default = ~/usr/local/lib
# INCFFT     : location of fftw headers 		default = ~/usr/local/include
# LIBNETCDF  : location of netdcdf libraries 	default = compiler/machine dependant
# INCNETCDF  : location of netcdf headers 		default = compiler/machine dependant
# 
# Dependencies: fftw, netcdf
#	FFTW is available here: http://www.fftw.org/
# 		FFTW is a C library with fortran headers
#	netcdf is available here: http://www.unidata.ucar.edu/software/netcdf/
# 		netcdf is a fortran library and must be compiled with the same fortran
# 		compiler you are using to be compatible
###################################################################
#  Specify where you want the resulting executable installed
###################################################################
INSTALLDIR=~/bin/
# DEBUG=debugompslow

###################################################################
#   Various compiler specific flags, may need to edit
###################################################################
NODENAME := $(shell uname -n)

# on hydro-c1:
# note on hydro-c1 /usr/local is not available on compute nodes so 
# the ifort libraries must be copied to an available directory (e.g. ~/)
# it is HIGHLY recommended that you set :
# LD_RUN_PATH=$LD_RUN_PATH:/usr/local/netcdf-4.1.3/ifort-12.0.5/lib:/home/gutmann/usr/local/lib:/usr/local/netcdf3-ifort/lib
# in your environment to point to the libraries you will need so the locations will be encoded in the 
#  compiled binary and you don't need to set LD_LIBRARY_PATH at runtime. 

#default parameters, overwritten with machine/compiler specific options below if known
F90=ifort
LIBFFT=/home/gutmann/usr/local/lib
INCFFT=/home/gutmann/usr/local/include

ifeq ($(NODENAME), Nomad-2.local)
	NODENAME=Nomad.local
endif
ifeq ($(NODENAME), Patthar.local)
	NODENAME=Nomad.local
endif
ifeq ($(NODENAME), vpn44.ucar.edu)
	NODENAME=Nomad.local
endif


ifeq ($(NODENAME), Nomad.local)
	F90=gfortran
	COMP=-fopenmp -lgomp -O3 -c -fimplicit-none
	LINK=-fopenmp -lgomp
	LIBFFT=/Users/gutmann/usr/local/lib
	INCFFT=/Users/gutmann/usr/local/include
	NCDF_PATH = /Users/gutmann/usr/local/
	LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
	INCNETCDF = -I$(NCDF_PATH)/include
	
	# DEBUG:=g+$(DEBUG)
endif

ifeq ($(F90), ifort)
	COMP= -openmp -liomp5 -fast -c -u -fpe0 #-fast-transcendentals -xhost
	LINK= -openmp -liomp5

# NCDF_PATH = /usr/local/netcdf-4.1.3/ifort-12.0.5
	NCDF_PATH = /home/gutmann/.usr/local
	LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
	INCNETCDF = -I$(NCDF_PATH)/include
endif



# on yellowstone:
ifeq ($(LMOD_FAMILY_COMPILER),gnu)
	F90=gfortran
	COMP=-fopenmp -lgomp -O3 -c -fimplicit-none
	LINK=-fopenmp -lgomp
	LIBFFT=/glade/u/home/gutmann/usr/local/lib
	INCFFT=/glade/u/home/gutmann/usr/local/include
	LIBNETCDF = -lnetcdff -lnetcdf
	INCNETCDF =  # netcdf includes are setup by the yellowstone module system
endif 
ifeq ($(LMOD_FAMILY_COMPILER),intel)
	F90=ifort
	COMP= -openmp -liomp5 -fast -c -u -fpe0 -fast-transcendentals -xhost
#	COMP= -openmp -liomp5 -debug -c -fast -u -check all -check noarg_temp_created -traceback -fpe0
	LINK= -openmp -liomp5
	LIBFFT=/glade/u/home/gutmann/usr/local/lib
	INCFFT=/glade/u/home/gutmann/usr/local/include
	LIBNETCDF = -lnetcdff -lnetcdf
	INCNETCDF =  # netcdf includes are setup by the yellowstone module system
endif 
ifeq ($(LMOD_FAMILY_COMPILER),pgi)
	F90=pgf90
	COMP=-fast -mp -c
	LINK=-mp
	LIBFFT=/glade/u/home/gutmann/usr/local/lib
	INCFFT=/glade/u/home/gutmann/usr/local/include
	LIBNETCDF = -lnetcdff -lnetcdf
	INCNETCDF =  # netcdf includes are setup by the yellowstone module system
endif	


# F90=/usr/local/ifort-2013.4.183/bin/ifort
# Various debuging options.  Set the DEBUG variable with "make DEBUG=debugslow" etc.
ifeq ($(DEBUG), debugslow)
	ifeq ($(F90), ifort)
		COMP= -debug -debug-parameters all -traceback -ftrapuv -g -fpe0 -c -u -check all -check noarg_temp_created -CB
		LINK=  
	endif
	ifeq ($(F90), gfortran)
		COMP= -c -g -fbounds-check -fbacktrace -finit-real=nan
		LINK=  
	endif
endif
ifeq ($(DEBUG), debug)
	ifeq ($(F90), ifort)
		COMP= -debug -c -fast -u -check all -check noarg_temp_created -traceback -fpe0 -fast-transcendentals -xhost
		LINK=  
	endif
	ifeq ($(F90), gfortran)
		COMP= -c -g -fbounds-check -fbacktrace -finit-real=nan
		LINK=  
	endif
endif
ifeq ($(DEBUG), debugompslow)
	ifeq ($(F90), ifort)
		COMP= -openmp -liomp5 -debug -debug-parameters all -traceback -ftrapuv -g -fpe0 -c -u -check all -check noarg_temp_created -CB
		COMP= -openmp -liomp5 -debug -c -u  -fpe0 -traceback -check all -check noarg_temp_created
		LINK= -openmp -liomp5
	endif
	ifeq ($(F90), gfortran)
		COMP= -fopenmp -lgomp -c -g -fbounds-check -fbacktrace -finit-real=nan
		LINK= -fopenmp -lgomp  
	endif
endif
ifeq ($(DEBUG), debugomp)
	ifeq ($(F90), ifort)
		COMP= -openmp -liomp5 -debug -c -fast -u -traceback -fpe0 -fast-transcendentals -xhost
		LINK= -openmp -liomp5
	endif
	ifeq ($(F90), gfortran)
		COMP= -fopenmp -lgomp -c -g -fbounds-check -fbacktrace -finit-real=nan
		LINK= -fopenmp -lgomp  
	endif
endif

# ifeq ($(DEBUG), profile)
# 	COMP= -openmp -liomp5 -p -c -fast -u -fast-transcendentals -xhost
# 	LINK= -openmp -liomp5 -p
# endif

###################################################################
###################################################################
# 
# Should not need to edit anything below this line
# 
###################################################################
###################################################################
# this seems like it should be a solution to the missing library problem, but it doesn't work
# LFLAGS=$(LINK) ${LIBNETCDF} -L${LIBFFT} -openmp-link=static -static-intel -Bstatic
# instead copy required libraries into a directory accessible on compute nodes and set LD_RUN_PATH e.g.
# export LD_RUN_PATH=$LD_RUN_PATH:/path/to/libraries/lib:/home/gutmann/usr/local/lib
LFLAGS=$(LINK) ${LIBNETCDF} -L${LIBFFT}
FFLAGS=$(COMP) ${INCNETCDF} -I${INCFFT}

OBJS=advect.o io_routines.o data_structures.o init.o mp_driver.o time_step.o boundary.o model_tracking.o \
		mp_thompson.o cu_driver.o cu_tiedtke.o driver.o output.o wind.o linear_winds.o geo_reader.o lsm.o

IDEALOBJS=advect.o io_routines.o data_structures.o init.o mp_driver.o \
		mp_thompson.o ideal.o output.o wind.o linear_winds.o geo_reader.o

WINDOBJS=io_routines.o data_structures.o init.o tests/test_wind.o wind.o linear_winds.o output.o geo_reader.o

GEOOBJS=io_routines.o data_structures.o tests/test_geo.o geo_reader.o

###################################################################
#   User facing rules
###################################################################

all:swim_real

install:swim_real
	cp swim_real ${INSTALLDIR}

clean:
	rm *.o *.mod

allclean:cleanall

cleanall: clean
	rm tests/*.o tests/*.mod swim_real swim_ideal geo_test wind_test # test_init
tests: geo_test wind_test #test_init

swim_real:${OBJS}
	${F90} ${LFLAGS} ${OBJS} -o swim_real  -lm -lfftw3

real:${OBJS} 
	${F90} ${LFLAGS} ${OBJS} -o swim_real  -lm -lfftw3

ideal:${IDEALOBJS} 
	${F90} ${LFLAGS} ${IDEALOBJS} -o swim_ideal  -lm -lfftw3

###################################################################
#   test cases
###################################################################
geo_test:${GEOOBJS}
	${F90} ${LFLAGS} ${GEOOBJS} -o geo_test

wind_test:${WINDOBJS}
	${F90} ${LFLAGS} ${WINDOBJS} -o wind_test -lfftw3 -lm


###################################################################
#   driver code for 
###################################################################

driver.o:driver.f90 advect.o data_structures.o init.o mp_driver.o output.o wind.o boundary.o cu_driver.o
	${F90} ${FFLAGS} driver.f90

ideal.o:ideal.f90 advect.o data_structures.o init.o mp_driver.o output.o wind.o boundary.o
	${F90} ${FFLAGS} ideal.f90

###################################################################
#   Core initial and boundary condition and time steping
###################################################################

init.o:init.f90 data_structures.o io_routines.o geo_reader.o model_tracking.o
	${F90} ${FFLAGS} init.f90

time_step.o:time_step.f90 data_structures.o mp_driver.o wind.o advect.o output.o lsm.o cu_driver.o
	${F90} ${FFLAGS} time_step.f90

boundary.o:boundary.f90 data_structures.o io_routines.o wind.o geo_reader.o output.o linear_winds.o
	${F90} ${FFLAGS} boundary.f90


###################################################################
#   I/O routines
###################################################################

output.o:output.f90 data_structures.o io_routines.o
	${F90} ${FFLAGS} output.f90

io_routines.o:io_routines.f90 data_structures.o
	${F90} ${FFLAGS} io_routines.f90

geo_reader.o:geo_reader.f90 data_structures.o
	${F90} ${FFLAGS} geo_reader.f90

###################################################################
#   Microphysics code
###################################################################

mp_driver.o:mp_thompson.o data_structures.o mp_driver.f90
	${F90} ${FFLAGS} mp_driver.f90

mp_thompson.o:mp_thompson.f90
	${F90} ${FFLAGS} mp_thompson.f90
	
###################################################################
#   Convection code
###################################################################
cu_driver.o:cu_tiedtke.o data_structures.o cu_driver.f90
	${F90} ${FFLAGS} cu_driver.f90

cu_tiedtke.o:cu_tiedtke.f90
	${F90} ${FFLAGS} cu_tiedtke.f90

###################################################################
#   Land Surface code
###################################################################
lsm.o: lsm.f90 data_structures.o
	${F90} ${FFLAGS} lsm.f90

###################################################################
#   Various wind/advection related code
###################################################################
advect.o:advect.f90 data_structures.o
	${F90} ${FFLAGS} advect.f90

wind.o:wind.f90 linear_winds.o
	${F90} ${FFLAGS} wind.f90

linear_winds.o:linear_winds.f90 io_routines.o
	${F90} ${FFLAGS} linear_winds.f90


###################################################################
#   Generic data structures, used by almost everything
###################################################################
data_structures.o:data_structures.f90
	${F90} ${FFLAGS} data_structures.f90

model_tracking.o:model_tracking.f90
	${F90} ${FFLAGS} model_tracking.f90

###################################################################
#   Unit tests (may only work at their git tags?)
###################################################################


tests/test_wind.o:tests/test_wind.f90 wind.o linear_winds.o
	${F90} ${FFLAGS} tests/test_wind.f90 -o tests/test_wind.o

tests/test_geo.o:tests/test_geo.f90 geo_reader.o data_structures.o
	${F90} ${FFLAGS} tests/test_geo.f90 -o tests/test_geo.o

# NOTE: too many changes in data structures/init have definitely broken this test, could be fixed if necessary, 
#       but not worth the time right now. 
# test_init:tests/test_init.f90 init.o
# 	${F90} ${FFLAGS} tests/test_init.f90 io_routines.f90 data_structures.f90 init.f90 -o test_init

