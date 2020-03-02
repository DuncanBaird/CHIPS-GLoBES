PROGNAME      = covey
SOURCES       = covey.C
INCLUDES      =
OBJECTS       = $(patsubst %.C, %.o, $(SOURCES))
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTCFLAGS    = $(shell root-config --cflags)
CXX           = g++

LIBS          = $(ROOTGLIBS)
LDFLAGS       =
CCFLAGS       = -g -O3 $(ROOTCFLAGS)

$(PROGNAME):    $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)

%.o : %.C $(INCLUDES)
	$(CXX) ${CCFLAGS} -c -o $@ $<

test:
	@echo $(ROOTCFLAGS)

clean:	
	-rm -f ${PROGNAME} ${OBJECTS}
