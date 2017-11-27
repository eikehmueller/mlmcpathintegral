EIGEN_DIR=/home/em459/Library/eigen-3.3.4/include/eigen3/
CFLAGS_OPT=-O3 -std=c++11 -I$(EIGEN_DIR)
CFLAGS_DEBUG=-g -Wall -O0 -std=c++11 -I$(EIGEN_DIR)
LFLAGS=
CPP=g++

#DEBUG=True

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
	driver.o

all: $(MAIN)

montecarlo.o: \
	quantityofinterest.o \

driver.o: \
	harmonicoscillatoraction.o \
	quantityofinterest.o \
	twolevelmetropolissampler.o \
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
