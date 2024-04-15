# geophysics-netcdf

## Description
- Header-only C++ library for handling Geoscience Australia's NetCDF based geophysical line data format which is still under development.

## Author
- Ross C Brodie, Geoscience Australia

## Dependencies
- NetCDF C library available from: https://www.unidata.ucar.edu/software/netcdf
- NetCDF C++ bindings library (netcdf-cxx4) available from: https://github.com/Unidata/netcdf-cxx4.git
- The cpp-utils utilities libary available from: https://github.com/rcb547/cpp-utils.git

## CMAKE
- The geophysics-netcdf library is not intended to be a cmake TOP_LEVEL project
- The external dependencies packages cpp-utils, netcdf-cxx4 and netcdf need to be loaded by a higher level cmake project.  See for example https://github.com/GeoscienceAustralia/ga-aem how the project is used.
