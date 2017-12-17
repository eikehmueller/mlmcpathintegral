# Include local, machine dependent Makefile.
# Use this for machine-dependent settings by copying over one of the files
# local_XXX.mk
include local.mk

CFLAGS_OPT=-O3 -std=c++11 -I$(EIGEN_DIR)
CFLAGS_DEBUG=-g -Wall -O0 -std=c++11 -I$(EIGEN_DIR)
LFLAGS=
CXX=g++

SOURCES=$(wildcard *.hh *.cc)
MAIN=driver.x
OBJS=$(patsubst %.cc,%.o,$(wildcard *.cc))


ifeq ($(DEBUG),True)
	CFLAGS=$(CFLAGS_DEBUG)
else
	CFLAGS=$(CFLAGS_OPT)
endif

all: $(MAIN)

# Sort out dependencies, see
# http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/
DEPDIR := .dependencies
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
POSTCOMPILE = @mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

%.o: %.cc %.hh
	$(CXX) $(DEPFLAGS) $(CFLAGS) -c -o $@ $<
	$(POSTCOMPILE)

%.o: %.cc
	$(CXX) $(DEPFLAGS) $(CFLAGS) -c -o $@ $<
	$(POSTCOMPILE)

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(basename $(SOURCES))))

$(MAIN): $(OBJS) config.h
	$(CXX) $(LFLAGS) -o $(MAIN) $(OBJS)

doc: $(SOURCES)
	doxygen doxygen.cfg

.phony: clean
clean:
	rm -rf *~ *.o $(DEPDIR) $(MAIN)

.phony: cleandoc
cleandoc:
	rm -rf doc/
