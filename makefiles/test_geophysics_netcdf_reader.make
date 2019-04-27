SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o
.DEFAULT_GOAL := allclean

includes  = -I$(srcdir)
includes += -I$(cpputilssrc)
includes += -I$(marray_include)

cxxflags  += -DUSEGLOBALSTACKTRACE
cxxflags  += -D_MPI_ENABLED

libs       =  -lnetcdf -lnetcdf_c++4 -lgdal -lCGAL_Core
executable =  $(exedir)/test_geophysics_netcdf_reader.exe

objects  = $(cpputilssrc)/general_utils.o
objects += $(cpputilssrc)/file_utils.o
objects += $(cpputilssrc)/gdal_utils.o
objects += $(cpputilssrc)/cgal_utils.o
objects += $(srcdir)/test_geophysics_netcdf_reader.o

%.o : %.cpp
	@echo ' '
	@echo 'Compiling ' $<
	$(mpicxx) -c $(includes) $(cxxflags) $< -o $@

compile: $(objects)

link: $(objects)
	mkdir -p $(exedir)
	@echo ' '
	@echo Linking
	$(mpicxx) $(objects) $(libs) -o $(executable)

clean:
	@echo ' '
	@echo Cleaning
	rm -f $(objects)
	rm -f $(executable)

all: compile link
allclean: clean compile link

