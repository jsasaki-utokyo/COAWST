#!/bin/bash
# build_mct.sh for building MCT using Intel compiler
export OPTIM="-O3 -ip -fp-model precise -xHost -I${HOME}/local/openmpi/4.0.4/include"
#export OPTIM="-O3 -ip -fp-model precise -xHost -I/home/app/intel/intel2019_up4/parallel_studio_xe_2019/compilers_and_libraries_2019/linux/mpi/intel64/include"

export CC=icc
export CXX=icpc
export F77=ifort
export FC=ifort
export F90=ifort
export CFLAGS=${OPTIM}
export CXXFLAGS=${OPTIM}
export FFLAGS=${OPTIM}
export FCFLAGS=${OPTIM}
export CPP="icc -E"
export CXXCPP="icpc -E"
./configure --prefix="${HOME}/local/mct-intel"
make
make install
