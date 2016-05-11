#!/bin/tcsh

module load netcdf/4.3.2
module load intel-cc/14.1.106
#module load gcc/5.2.0
module list

make -f $1 $2

