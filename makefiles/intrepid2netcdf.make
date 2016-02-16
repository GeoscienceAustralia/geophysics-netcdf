SHELL = /bin/sh

include mydefaultmakesettings.make

#cxxflags   += -DUSEGLOBALSTACKTRACE
libs       =  -lnetcdf -lnetcdf_c++4
executable =  $(exedir)/intrepid2netcdf.exe

objects += $(objdir)/general_utils.o
objects += $(objdir)/file_utils.o
objects += $(objdir)/blocklanguage.o
objects += $(objdir)/intrepid2netcdf.o

$(objects): $(objdir)/%.o: $(srcdir)/%.cpp
	@echo ' '
	@echo Compiling $<
	$(CXX) -c $(includes) $(cxxflags) $< -o $@

compile: $(objects)

link: $(objects)
	@echo ' '
	@echo Linking
	mkdir -p $(exedir)
	$(CXX) $(objects) $(libs) -o $(executable)

clean:
	@echo ' '
	@echo Cleaning
	rm -f $(objects)
	rm -f $(executable)

all: compile link

allclean: clean compile link

