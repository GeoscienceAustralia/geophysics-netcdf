SHELL = /bin/sh

srcdir     = ../src/cxx
objdir     = ./obj
includes   = -I$(srcdir)
cxxflags   = -std=c++11 -O3 -Wall

#cxxflags   += -DUSEGLOBALSTACKTRACE
cxxflags   += -D_MPI_ENABLED
libs       =  -lnetcdf -lnetcdf_c++4 -lCGAL_Core
executable =  $(exedir)/intrepid2netcdf.exe

objects += $(objdir)/general_utils.o
objects += $(objdir)/file_utils.o
objects += $(objdir)/blocklanguage.o
objects += $(objdir)/cgal_utils.o
objects += $(objdir)/geophysics_netcdf.o
objects += $(objdir)/intrepid2netcdf.o

$(objects): $(objdir)/%.o: $(srcdir)/%.cpp
	@echo ' '
	@echo Compiling $<
	$(mpicxx) -c $(includes) $(cxxflags) $< -o $@

compile: $(objects)

link: $(objects)
	@echo ' '
	@echo Linking
	mkdir -p $(exedir)
	$(mpicxx) $(objects) $(libs) -o $(executable)

clean:
	@echo ' '
	@echo Cleaning
	rm -f $(objects)
	rm -f $(executable)

all: compile link

allclean: clean compile link

.SUFFIXES:
.SUFFIXES: .cpp .o
.DEFAULT_GOAL := allclean

