ROOTCONFIG   := root-config
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --glibs)

CXX           = g++
CXXFLAGS      = -O2 -Wall -fPIC

CXXFLAGS     += $(ROOTCFLAGS)
LDFLAGS      += $(ROOTLDFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)
#LIB_FLAGS +=  -lboost_python

SVCSOURCE     = *.cxx  

all: ${SVCSOURCE} 
	${CXX} ${SVCSOURCE} $(GLIBS) ${CXXFLAGS}  -o test.exe



clean:
	rm -f *.o *.out *.so *.exe 
