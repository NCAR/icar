#!/usr/bin/env bash

set -e
set -x

export FC=gfortran-9
# see link for size of runner
# https://docs.github.com/en/actions/using-github-hosted-runners/about-github-hosted-runners#supported-runners-and-hardware-resources
export JN=-j

if [ -z "$WORKDIR" ]; then
    export WORKDIR=$HOME/workdir
    mkdir -p $WORKDIR
fi

if [ -z "$INSTALLDIR" ]; then
    export INSTALLDIR=$HOME/installdir
    mkdir -p $INSTALLDIR
fi
export LD_LIBRARY_PATH=${INSTALLDIR}/lib:${LD_LIBRARY_PATH}


function install_szip {
    echo install_szip
    cd $WORKDIR
    wget --no-check-certificate -q http://www.hdfgroup.org/ftp/lib-external/szip/2.1.1/src/szip-2.1.1.tar.gz
    tar -xzf szip-2.1.1.tar.gz
    cd szip-2.1.1
    ./configure --prefix=$INSTALLDIR &> config.log
    make &> make.log
    make install ${JN}
}

function install_hdf5 {
    echo install_hdf5
    cd $WORKDIR
    wget --no-check-certificate -q https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.7/src/hdf5-1.10.7.tar.gz
    tar -xzf hdf5-1.10.7.tar.gz
    cd hdf5-1.10.7
    # FCFLAGS="-DH5_USE_110_API" ./configure --prefix=$INSTALLDIR &> config.log
    ./configure --prefix=$INSTALLDIR #&> config.log
    make
    # CFLAGS=-DH5_USE_110_API make
    # (CFLAGS=-DH5_USE_110_API make | awk 'NR%100 == 0')
    make install ${JN}
    export LIBDIR=${INSTALLDIR}/lib
}

function install_netcdf_c {
    echo install_netcdf_c
    cd $WORKDIR
    wget --no-check-certificate -q ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-c-4.8.0.tar.gz
    tar -xzf netcdf-c-4.8.0.tar.gz
    cd netcdf-c-4.8.0
    ./configure --prefix=$INSTALLDIR #&> config.log
    make #&> make.log
    make install
}

function install_netcdf_fortran {
    echo install_netcdf_fortran
    cd $WORKDIR
    wget --no-check-certificate -q ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.5.3.tar.gz
    tar -xzf netcdf-fortran-4.5.3.tar.gz
    cd netcdf-fortran-4.5.3
    ./configure --prefix=$INSTALLDIR #&> config.log
    make #&> make.log
    make install
}

function icar_dependencies {
    echo icar_dependencies

    sudo apt-get update
    sudo apt-get install libcurl4-gnutls-dev
    sudo apt-get install libfftw3-dev
    # Installing HDF5 currently not working for NetCDF
    # sudo apt-get install libhdf5-dev libhdf5-openmpi-dev

    export CPPFLAGS="$CPPFLAGS -I${INSTALLDIR}/include"
    export LDFLAGS="$LDFLAGS -L${INSTALLDIR}/lib"

    # Install szip (used by hdf5)
    install_szip
    # Install HDF5
    install_hdf5

    # Install NetCDF-C
    install_netcdf_c
    # Install NetCDF fortran
    install_netcdf_fortran

    # put installed bin directory in PATH
    export PATH=${INSTALLDIR}/bin:$PATH
}

function icar_install {
    echo icar_install
    cd ${GITHUB_WORKSPACE}

    export NETCDF=${INSTALLDIR}
    export FFTW=/usr
    export JN=-j

    # CAF_MODE=single tells it to compile with gfortran -fcoarray=single
    # test serial build
    make -C src clean; VERBOSE=1 make -C src     CAF_MODE=single MODE=debugslow
    # test parallel builds with different compile settings
    # make -C src clean; VERBOSE=1 make -C src ${JN} CAF_MODE=single MODE=debugslow
    # make -C src clean; VERBOSE=1 make -C src ${JN} CAF_MODE=single MODE=debug
    # make -C src clean; VERBOSE=1 make -C src ${JN} CAF_MODE=single MODE=debugompslow
    # make -C src clean; VERBOSE=1 make -C src ${JN} CAF_MODE=single MODE=debugomp
    # make -C src clean; VERBOSE=1 make -C src ${JN} CAF_MODE=single MODE=profile
    # # make -C src clean; make -C src ${JN} CAF_MODE=single MODE=fast
    # make -C src clean; VERBOSE=1 make -C src ${JN} CAF_MODE=single
    # make -C src ${JN} CAF_MODE=single test

    echo "icar install succeeded"

}

function icar_script {

    cd ${GITHUB_WORKSPACE}
    cd ./src

    cp ../run/complete_icar_options.nml ./icar_options.nml
    export OMP_NUM_THREADS=2

    ./icar
    # for i in *_test; do
    #     echo $i
    #     ./${i}
    # done

    echo "icar script succeeded"
    cd ..
}

function gen_test_run_data {
    cd ${GITHUB_WORKSPACE}/tests
    python --version
    python -m pip install --upgrade pip
    pip install setuptools numpy pandas netCDF4 xarray
    python gen_ideal_test.py
}

function execute_test_run {
    cp ${GITHUB_WORKSPACE}/src/icar ${GITHUB_WORKSPACE}/tests/
    cd ${GITHUB_WORKSPACE}/tests
    ./icar icar_options.nm
}

function icar_after_success {
  echo icar_after_success
  echo "icar build succeeded"
}

function icar_after_failure {
  echo icar_after_failure
  echo "icar build failed"
}

for func in "$@"
do
    case $func in
        icar_dependencies)
            icar_dependencies;;
        icar_install)
            icar_install;;
	icar_script)
	    icar_script;;
	gen_test_run_data)
	    gen_test_run_data;;
	execute_test_run)
	    execute_test_run;;
        *)
            echo "$func unknown"
    esac
done
