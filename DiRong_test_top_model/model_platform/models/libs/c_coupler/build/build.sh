#!/bin/bash
CPPFLAGS="-traditional-cpp -DLINUX -DNO_SHR_VMATH -DMPI_P"
NETCDFINC="  -g -I/public1/soft/netcdf/4.4.1-parallel-icc18/include "
NETCDFLIB="  -L/public1/soft/netcdf/4.4.1-parallel-icc18/lib -lnetcdff -lnetcdf "
MPIINC="  -I/public1/soft/intel/2018/compilers_and_libraries_2018.2.199/linux/mpi/intel64/include "
MPILIB="  -L/public1/soft/intel/2018/compilers_and_libraries_2018.2.199/linux/mpi/intel64/lib "

export CC=mpiicc
export CXX=mpiicpc
export FC=mpiifort
#export CPP=/usr/bin/cpp
#export CFLAGS="-O2 -DFORTRANUNDERSCORE -g"
export CFLAGS=" -DFORTRANUNDERSCORE -g"
#export CXXFLAGS="-O2 -c -DFORTRANUNDERSCORE -g"
export CXXFLAGS=" -O0 -c -DFORTRANUNDERSCORE -traceback -g"
export FFLAGS="-g -free -O0 -c -i4  -r8 -convert big_endian -assume byterecl -fp-model precise"
export INCLDIR=" ${NETCDFINC} ${MPIINC} "
export SLIBS=" ${NETCDFLIB} ${MPILIB} "
export CPPFLAGS="${CPPFLAGS} ${INCLDIR} "

make -j 8
