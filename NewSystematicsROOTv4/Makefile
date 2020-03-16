# GLoBES -- General LOng Baseline Experiment Simulator
# (C) 2002 - 2007,  The GLoBES Team
#
# GLoBES is mainly intended for academic purposes. Proper
# credit must be given if you use GLoBES or parts of it. Please
# read the section 'Credit' in the README file.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# The configure script should have replaced anything as required to
# obtain a working Makefile, which may serve as a template for own
# applications.
#
# This simple Makefile is for the GLoBES examples
#
# Compile example files with ``make example1'' etc.
#
# This Makefile assumes dynamic libraries, installed to either
# the default prefix /usr/local/ or to a user-defined directory 
# called ${prefix}.
#
# For linking against a specific version of GLoBES, libglobes.so can be 
# replaced by the respective library, such as libglobes.so.0.0.1



prefix = /usr/local
exec_prefix = ${prefix}
libdir = ${exec_prefix}/lib
globesconf= $(exec_prefix)/bin/globes-config

CC=gcc
CXX=g++

ROOTLIBS = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
ROOTCFLAGS = $(shell root-config --cflags)

local_CFLAGS = -g -O4 $(ROOTCFLAGS)

INCFLAGS:=$(shell $(globesconf) --include)
local_LDFLAGS:=$(shell $(globesconf) --libs) $(ROOTGLIBS)
local_LTLDFLAGS:=$(shell $(globesconf) --ltlibs)

BIN = example-tour example-static customsimulation testsimulation example2 example3 example4 example5 example6
OBJ = example-tour.o example-tour-s.o customsimulation.o testsimulation.o example2.o example3.o example4.o \
      example5.o example6.o myio.o

all: $(BIN)

example-tour: example-tour.o
	$(CC) example-tour.o -o example-tour $(LDFLAGS) $(local_LDFLAGS)

example-static: example-tour.o
	$(CC) example-tour.o -o example-static -static $(LDFLAGS) $(local_LDFLAGS) 

customsimulation: customsimulation.o myio.o
	$(CXX) customsimulation.o myio.o -o customsimulation  $(LDFLAGS) $(local_LDFLAGS)

testsimulation: testsimulation.o myio.o
	$(CXX) testsimulation.o myio.o -o testsimulation  $(LDFLAGS) $(local_LDFLAGS)

example2: example2.o myio.o
	$(CC)  example2.o myio.o -o example2  $(LDFLAGS) $(local_LDFLAGS)

example3: example3.o myio.o
	$(CC) example3.o myio.o -o example3  $(LDFLAGS) $(local_LDFLAGS)

example4: example4.o myio.o
	$(CC) example4.o myio.o -o example4  $(LDFLAGS) $(local_LDFLAGS)

example5: example5.o myio.o
	$(CC) example5.o myio.o -o example5  $(LDFLAGS) $(local_LDFLAGS)

example6: example6.o myio.o
	$(CC) example6.o myio.o -o example6  $(LDFLAGS) $(local_LDFLAGS)

%.o : %.c
	$(CC) $(CFLAGS) $(local_CFLAGS) -c $< $(INCFLAGS)
%.o : %.cc
	$(CXX) $(CFLAGS) $(local_CFLAGS) -c $< $(INCFLAGS)
.PHONY: clean
clean:
	rm -f $(BIN) $(OBJ)