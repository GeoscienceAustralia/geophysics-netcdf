#!/bin/sh
echo '========================================================================'
echo '========================================================================'
echo '========================================================================'

export compiler=$1
export makemode=$2
export srcdir='../src'
export cpputilssrc='../submodules/cpp-utils/src'
export marray_include='../submodules/marray/include/andres'


if [ $compiler == 'intel' ] ; then
	echo 'Building with Intel compiler'
	module load intel-compiler
	export cxx=icpc
	export cxxflags='-std=c++17 -O3 -Wall'
	export exedir='../bin/gadi/intel'
elif [ $compiler == 'gnu' ] ; then
	echo 'Building with GCC compiler'
	module load gcc/system
	export cxx=g++
	export cxxflags='-std=c++17 -O3 -Wall -Wno-unknown-pragmas'
	export exedir='../bin/gadi/gnu'
else 
	echo 'Unknow compiler ' $compiler
	exit
fi

export mpicxx=mpiCC
export HAVE_GDAL=1
export HAVE_CGAL=0
module load openmpi/4.0.1
module load netcdf/4.7.1

if [ $HAVE_GDAL == 1 ] ; then
	echo 'Building with GDAL'
	module load gdal/3.0.2
fi

if [ $HAVE_CGAL == 1 ] ; then
	echo 'Building with CGAL'
	module load cgal
fi



module list
echo ---------------------------------------
echo cxx = $cxx
echo mpicxx = $mpicxx ... which is ...
$mpicxx -showme

echo HAVE_GDAL = $HAVE_GDAL
echo HAVE_CGAL = $HAVE_CGAL
echo ---------------------------------------

make -f intrepid2netcdf.make $makemode
make -f geophysicsnc2shape.make $makemode
#make -f aseggdf2netcdf.make $makemode
#make -f test_geophysics_netcdf_reader.make $makemode

