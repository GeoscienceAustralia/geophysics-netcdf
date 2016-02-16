srcdir     = ../src/cxx
objdir     = ./obj
includes   = -I$(srcdir)
cxxflags   = -std=c++11 -O3 -Wall

ifeq ($(CXX),g++)
	#GNU compiler on raijin.nci.org.au
	exedir     = ../bin/raijin/gnu
else
	#Intel compiler on raijin.nci.org.au
	exedir     = ../bin/raijin/intel
	cxxflags  += -diag-disable remark
endif

.SUFFIXES:
.SUFFIXES: .cpp .o
.DEFAULT_GOAL := allclean

