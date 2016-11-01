#!/bin/tcsh

module load cgal/4.6.3
module load netcdf/4.3.3.1p
module load intel-cc/14.1.106
module list

make -f $1 $2

