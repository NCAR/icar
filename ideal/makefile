OBJS=advect.o io_routines.o data_structures.o init.o mp_driver.o \
		mp_thompson.o driver.o output.o wind.o linear_winds.o
WINDOBJS=io_routines.o data_structures.o init.o tests/test_wind.o wind.o linear_winds.o output.o

DEBUGFLAGS = -g
FOFLAGS= -c
# OPTFLAGS= -O3
FFLAGS=$(DEBUGFLAGS)
LFLAGS=$(FFLAGS)
ifeq ($(LMOD_FAMILY_COMPILER),gnu)
	F90=gfortran
	FFLAGS=-fopenmp -lgomp -O3
	LFLAGS=-fopenmp -lgomp -O3
	# FFLAGS=-fopenmp -lgomp -g
endif 
ifeq ($(LMOD_FAMILY_COMPILER),intel)
	F90=ifort
	FFLAGS=-fpp -openmp -liomp5 -fast
	LFLAGS=-fpp -openmp -liomp5
	# FFLAGS=-fpp -openmp -liomp5 -g
	# LFLAGS=-fpp -openmp -liomp5 -g
endif 
ifeq ($(LMOD_FAMILY_COMPILER),pgi)
	F90=pgf90
	# FFLAGS=-fast -mp
endif	

ideal:${OBJS} 
	${F90} ${LFLAGS} ${OBJS} -o swim_ideal

test_init:tests/test_init.f90 init.o
	${F90} ${FFLAGS} tests/test_init.f90 io_routines.f90 data_structures.f90 init.f90 -o test_init

test_winds:${WINDOBJS}
	${F90} ${LFLAGS} -L/glade/u/home/gutmann/usr/local/lib ${WINDOBJS} -o test_wind -lfftw3 -lm

driver.o:driver.f90 advect.o data_structures.o init.o mp_driver.o output.o wind.o
	${F90} ${FOFLAGS} ${FFLAGS} driver.f90

advect.o:advect.f90 data_structures.o
	${F90} ${FOFLAGS} ${FFLAGS} advect.f90
	
io_routines.o:io_routines.f90 data_structures.o
	${F90} ${FOFLAGS} ${FFLAGS} io_routines.f90

data_structures.o:data_structures.f90
	${F90} ${FOFLAGS} ${FFLAGS} data_structures.f90

init.o:init.f90 data_structures.o io_routines.o
	${F90} ${FOFLAGS} ${FFLAGS} init.f90

mp_thompson.o:mp_thompson.f90
	${F90} ${FOFLAGS} ${FFLAGS} mp_thompson.f90

mp_driver.o:mp_thompson.o data_structures.o mp_driver.f90
	${F90} ${FOFLAGS} ${FFLAGS} mp_driver.f90

output.o:output.f90 data_structures.o io_routines.o
	${F90} ${FOFLAGS} ${FFLAGS} output.f90
	
tests/test_wind.o:tests/test_wind.f90 wind.o linear_winds.o
	${F90} ${FOFLAGS} ${FFLAGS} tests/test_wind.f90 -o tests/test_wind.o

wind.o:wind.f90 linear_winds.o
	${F90} ${FOFLAGS} ${FFLAGS} wind.f90

linear_winds.o:linear_winds.f90
	${F90} ${FOFLAGS} ${FFLAGS} -I/glade/u/home/gutmann/usr/local/include linear_winds.f90


# %.o : %.f90
#         $(F90) $(FOLAGS) $(FFLAGS) $< -o $@

	
clean:
	rm *.o *.mod
