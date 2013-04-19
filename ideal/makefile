INSTALLDIR=~/bin/

DEBUGFLAGS = -g
FOFLAGS= -c
# OPTFLAGS= -O3
FFLAGS=$(DEBUGFLAGS)
LFLAGS=$(FFLAGS)
F90=ifort

LIBFFT=/home/gutmann/usr/local/lib
INCFFT=/home/gutmann/usr/local/include

ifeq ($(F90), ifort)
	COMP=-fpp -openmp -liomp5 -fast
	LINK=-fpp -openmp -liomp5
	NCDF_PATH = /usr/local/netcdf-4.1.3/ifort-12.0.5
	LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
	INCNETCDF = -I$(NCDF_PATH)/include
endif
COMP=-g
link=-g
# on yellowstone:
ifeq ($(LMOD_FAMILY_COMPILER),gnu)
	F90=gfortran
	COMP=-fopenmp -lgomp -O3
	LINK=-fopenmp -lgomp
	LIBFFT=/glade/u/home/gutmann/usr/local/lib
	INCFFT=/glade/u/home/gutmann/usr/local/include
	LIBNETCDF = -lnetcdff -lnetcdf
	INCNETCDF =  
endif 
ifeq ($(LMOD_FAMILY_COMPILER),intel)
	F90=ifort
	COMP=-fpp -openmp -liomp5 -fast
	LINK=-fpp -openmp -liomp5
	LIBFFT=/glade/u/home/gutmann/usr/local/lib
	INCFFT=/glade/u/home/gutmann/usr/local/include
	LIBNETCDF = -lnetcdff -lnetcdf
	INCNETCDF =  
endif 
ifeq ($(LMOD_FAMILY_COMPILER),pgi)
	F90=pgf90
	COMP=-fast -mp
	LINK=-mp
	LIBFFT=/glade/u/home/gutmann/usr/local/lib
	INCFFT=/glade/u/home/gutmann/usr/local/include
	LIBNETCDF = -lnetcdff -lnetcdf
	INCNETCDF =  
endif	


LFLAGS=$(LINK) ${LIBNETCDF} -L${LIBFFT}
FFLAGS=$(COMP) ${INCNETCDF} -I${INCFFT}

OBJS=advect.o io_routines.o data_structures.o init.o mp_driver.o \
		mp_thompson.o driver.o output.o wind.o linear_winds.o geo_reader.o
WINDOBJS=io_routines.o data_structures.o init.o tests/test_wind.o wind.o linear_winds.o output.o geo_reader.o
GEOOBJS=io_routines.o data_structures.o tests/test_geo.o geo_reader.o

all:ideal

install:ideal
	cp swim_ideal ${INSTALLDIR}

clean:
	rm *.o *.mod

allclean:cleanall
		
cleanall: 
	rm *.o *.mod tests/*.o tests/*.mod swim_ideal geo_test wind_test # test_init
	
tests: geo_test wind_test #test_init


ideal:${OBJS} 
	${F90} ${LFLAGS} ${OBJS} -o swim_ideal  -lm -lfftw3


# test_init:tests/test_init.f90 init.o
# 	${F90} ${FFLAGS} tests/test_init.f90 io_routines.f90 data_structures.f90 init.f90 -o test_init

geo_test:${GEOOBJS}
	${F90} ${LFLAGS} ${GEOOBJS} -o geo_test

wind_test:${WINDOBJS}
	${F90} ${LFLAGS} ${WINDOBJS} -o wind_test -lfftw3 -lm



driver.o:driver.f90 advect.o data_structures.o init.o mp_driver.o output.o wind.o
	${F90} ${FOFLAGS} ${FFLAGS} driver.f90

advect.o:advect.f90 data_structures.o
	${F90} ${FOFLAGS} ${FFLAGS} advect.f90
	
io_routines.o:io_routines.f90 data_structures.o
	${F90} ${FOFLAGS} ${FFLAGS} io_routines.f90

data_structures.o:data_structures.f90
	${F90} ${FOFLAGS} ${FFLAGS} data_structures.f90

init.o:init.f90 data_structures.o io_routines.o geo_reader.o
	${F90} ${FOFLAGS} ${FFLAGS} init.f90

geo_reader.o:geo_reader.f90 data_structures.o
	${F90} ${FOFLAGS} ${FFLAGS} geo_reader.f90

mp_thompson.o:mp_thompson.f90
	${F90} ${FOFLAGS} ${FFLAGS} mp_thompson.f90

mp_driver.o:mp_thompson.o data_structures.o mp_driver.f90
	${F90} ${FOFLAGS} ${FFLAGS} mp_driver.f90

output.o:output.f90 data_structures.o io_routines.o
	${F90} ${FOFLAGS} ${FFLAGS} output.f90
	
wind.o:wind.f90 linear_winds.o
	${F90} ${FOFLAGS} ${FFLAGS} wind.f90

linear_winds.o:linear_winds.f90
	${F90} ${FOFLAGS} ${FFLAGS} linear_winds.f90



tests/test_wind.o:tests/test_wind.f90 wind.o linear_winds.o
	${F90} ${FOFLAGS} ${FFLAGS} tests/test_wind.f90 -o tests/test_wind.o

tests/test_geo.o:tests/test_geo.f90 geo_reader.o data_structures.o
	${F90} ${FOFLAGS} ${FFLAGS} tests/test_geo.f90 -o tests/test_geo.o
