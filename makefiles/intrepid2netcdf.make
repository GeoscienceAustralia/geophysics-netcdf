SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o
.DEFAULT_GOAL := allclean

includes  = -I$(srcdir)
includes += -I$(cpputilssrc)
includes += -I$(marray_include)


cxxflags  += -DUSEGLOBALSTACKTRACE
cxxflags  += -D_MPI_ENABLED

libs       =  -lnetcdf -lnetcdf_c++4

executable =  $(exedir)/intrepid2netcdf.exe
objects  = $(cpputilssrc)/general_utils.o
objects += $(cpputilssrc)/file_utils.o
ifeq ($(HAVE_GDAL),1)
    cxxflags += -DHAVE_GDAL
    objects  += $(cpputilssrc)/gdal_utils.o
    libs     +=  -lgdal
endif
ifeq ($(HAVE_CGAL),1)
    cxxflags += -DHAVE_CGAL
    objects  += $(cpputilssrc)/cgal_utils.o
endif
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

