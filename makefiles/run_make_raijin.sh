#!/bin/sh

export srcdir='../src'
export cpputilssrc='../submodules/cpp-utils/src'
export marray_include='../submodules/marray/include/andres'

#Intel compiler on raijin.nci.org.au
module load geophysics-netcdf

export cxx=icpc
export mpicxx=mpiCC
export cxxflags='-std=c++11 -O3 -Wall -wd858 -diag-disable remark -D_GLIBCXX_USE_CXX11_ABI=0'
export exedir='../bin/raijin/intel'
module list

echo ---------------------------------------
echo cxx = $cxx
echo mpicxx = $mpicxx ... which is ...
$mpicxx -showme
echo ---------------------------------------

make -f intrepid2netcdf.make $1
#make -f aseggdf2netcdf.make $1
#make -f test_geophysics_netcdf_reader.make $1

