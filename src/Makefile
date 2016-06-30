ARCH = linux

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
ROOTGLIBS := $(shell root-config --glibs)

ifeq ($(ARCH),linux)
CXX = g++
CXXFLAGS = -g -O -Wno-deprecated -fPIC -m32 -fno-inline -Wno-write-strings
LD = g++
LDFLAGS = -O
SOFLAGS = -shared
INCLUDE = $(ROOTSYS)/include/
endif

CXXFLAGS += $(ROOTCFLAGS)
LIBS = $(ROOTLIBS) $(SYSLIBS)
GLIBS = $(ROOTGLIBS) $(SYSLIBS)

#-----------------------------------------------

CXXFILES = Environ.cxx RunInfo.cxx LevelTwo.cxx EventClass.cxx TCUbits.cxx RPscint.cxx BBCtiles.cxx KinBounds.cxx
HEADERS = Environ.h RunInfo.h LevelTwo.h EventClass.h TCUbits.h RPscint.h BBCtiles.h KinBounds.h

#-----------------------------------------------

All: RunInfo.so

RunInfo.so: LinkDef.h $(HEADERS) $(CXXFILES)
	@echo "generating dictionary $@..."
	@rootcint -f RunInfoDict.cxx -c $(HEADERS) LinkDef.h
	g++ $(CXXFILES) RunInfoDict.cxx -shared -o RunInfo.so $(CXXFLAGS) $(GLIBS)

clean:
	@rm RunInfo.so
