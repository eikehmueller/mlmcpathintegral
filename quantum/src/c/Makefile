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

MAIN=driver.x
OBJS=\
	harmonicoscillatoraction.o \
	quantityofinterest.o \
	twolevelmetropolissampler.o \
	montecarlo.o \
	parameters.o \
	driver.o

all: $(MAIN)

montecarlo.o: \
	quantityofinterest.o \

driver.o: \
	harmonicoscillatoraction.o \
	quantityofinterest.o \
	twolevelmetropolissampler.o \
	parameters.o \
	montecarlo.o

%.o: %.cc %.hh
	$(CPP) $(CFLAGS) -c -o $@ $<

%.o: %.cc
	$(CPP) $(CFLAGS) -c -o $@ $<

$(MAIN): $(OBJS)
	$(CPP) $(LFLAGS) -o $(MAIN) $(OBJS)

.phony: clean
clean:
	rm -rf *~ *.o $(MAIN)
