#!/bin/sh

module load cgal/4.6.3
module load netcdf/4.3.3.1p
module load openmpi/1.6.3
module list

export srcdir='../src/cpp'
export cpputilssrc='../submodules/cpp-utils/src'
export netcdfcxxdir='../submodules/netcdf-cxx4'

#GNU compiler on raijin.nci.org.au
#module load gcc/4.9.0 #netcdfc++ not linking with module load gcc/5.2.0 #module load gcc/6.2.0
#export cxx=g++
#export mpicxx=mpiCC
#export cxxflags='-std=c++11 -O3 -Wall -fdiagnostics-color=always
#export exedir='../bin/raijin/gnu'

#Intel compiler on raijin.nci.org.au
module load intel-cc/14.1.106
export cxx=icpc
export mpicxx=mpiCC
export cxxflags='-std=c++11 -O3 -Wall -diag-disable remark'
export exedir='../bin/raijin/intel'

mpiCC -showme
make -f intrepid2netcdf.make $1


