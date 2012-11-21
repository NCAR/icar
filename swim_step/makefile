#F90=gfortran
#FFLAGS=-x f95-cpp-input -fopenmp -D__OPENMP -lgomp -O3
#F2PYFLAGS=--fcompiler=gnu95 --f90flags='-x f95-cpp-input -fopenmp -D__OPENMP' -lgomp --opt="-O3" -DF2PY_REPORT_ON_ARRAY_COPY=1

F90=ifort
FFLAGS=-fpp -fast -openmp -c
INSTALLDIR=~/bin/python/lib/

F2PY=f2py
F2PYFLAGS=--fcompiler=intelem --f90flags='-fpp -openmp' -lguide --opt='-fast' -DF2PY_REPORT_ON_ARRAY_COPY=1


all: swim.so

swim.so:swim.pyf swim_all.f90
	${F2PY} -c swim.pyf ${F2PYFLAGS} swim_all.f90

swim.pyf:swim_all.f90
	${F2PY} -h swim.pyf -m swim swim_all.f90 --overwrite-signature
swim_all.f90:mp_simple.f90 advect.f90 py_mp_thompson.f90 swim_step.f90
	cat mp_simple.f90 advect.f90 py_mp_thompson.f90 swim_step.f90 >swim_all.f90
clean:
	rm *.o *.mod swim.so swim.pyf

install:swim.so
	cp swim.so ${INSTALLDIR}

#note this would be a better way to do it, but there seem to be issues with f2py
#OBJS=mp_simple.o py_mp_thompson.o advect.o
#
#swim.so:${OBJS} swim.pyf swim_step.f90
#	${F2PY} -c swim.pyf ${F2PYFLAGS} swim_step.f90
#
#swim.pyf:swim_step.f90
#	${F2PY} -h swim.pyf -m swim swim_step.f90 --overwrite-signature
#
#mp_simple.o:mp_simple.f90
#	${F90} ${FFLAGS} mp_simple.f90
#advect.o:advect.f90
#	${F90} ${FFLAGS} advect.f90
#py_mp_thompson.o:py_mp_thompson.f90
#	${F90} ${FFLAGS} py_mp_thompson.f90


