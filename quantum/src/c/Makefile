# Include local, machine dependent Makefile.
# Use this for machine-dependent settings by copying over one of the files
# local_XXX.mk
include local.mk

CFLAGS_OPT=-O3 -std=c++11 -I$(EIGEN_DIR)
CFLAGS_DEBUG=-g -Wall -O0 -std=c++11 -I$(EIGEN_DIR)
LFLAGS=
CPP=g++

ifeq ($(DEBUG),True)
	CFLAGS=$(CFLAGS_DEBUG)
else
	CFLAGS=$(CFLAGS_OPT)
endif

SOURCES=$(wildcard *.hh *.cc)

MAIN=driver.x
OBJS=$(patsubst %.cc,%.o,$(wildcard *.cc))

all: $(MAIN)

action.o: \
	path.o

montecarlo.o: \
	path.o \
	sampler.o \
	quantityofinterest.o 

hmcsampler.o: \
	sampler.o

twolevelmetropolissampler.o: \
	sampler.o

driver.o: \
	harmonicoscillatoraction.o \
	quantityofinterest.o \
	hmcsampler.o \
	twolevelmetropolissampler.o \
	parameters.o \
	montecarlo.o

%.o: %.cc %.hh
	$(CPP) $(CFLAGS) -c -o $@ $<

%.o: %.cc
	$(CPP) $(CFLAGS) -c -o $@ $<

$(MAIN): $(OBJS)
	$(CPP) $(LFLAGS) -o $(MAIN) $(OBJS)

doc: $(SOURCES)
	doxygen doxygen.cfg

.phony: clean
clean:
	rm -rf *~ *.o $(MAIN)

.phony: cleandoc
cleandoc:
	rm -rf doc/
