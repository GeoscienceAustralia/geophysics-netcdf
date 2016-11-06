#!/bin/sh

module load cgal/4.6.3
module load netcdf/4.3.3.1p
module load openmpi/1.6.3
module list

#Intel compiler on raijin.nci.org.au
module load intel-cc/14.1.106
export cxx=icpc
export mpicxx=mpiCC
export cxxflags='-std=c++11 -O3 -Wall -diag-disable remark'
export exedir='../bin/raijin/intel'

make -f $1 $2

