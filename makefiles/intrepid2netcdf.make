SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o
.DEFAULT_GOAL := allclean

includes   = -I$(cpputilssrc) -I$(srcdir)

#cxxflags  += -DUSEGLOBALSTACKTRACE
cxxflags   += -D_MPI_ENABLED

libs       =  -lnetcdf -lnetcdf_c++4 -lCGAL_Core
executable =  $(exedir)/intrepid2netcdf.exe

objects += $(cpputilssrc)/general_utils.o
objects += $(cpputilssrc)/file_utils.o
objects += $(cpputilssrc)/blocklanguage.o
objects += $(cpputilssrc)/cgal_utils.o
objects += $(srcdir)/geophysics_netcdf.o
objects += $(srcdir)/intrepid2netcdf.o

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

