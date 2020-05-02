# Include local, machine dependent Makefile.
# Use this for machine-dependent settings by copying over one of the files
# local_XXX.mk

COM_MESSAGE   = "Compiling"
LINK_MESSAGE  = "Linking"

ifeq ($(NOCOMCOLOR),1)
  COM_COLOR   =
  OBJ_COLOR   =
  NO_COLOR    =
else
  COM_COLOR   = \033[0;1m
  OBJ_COLOR   = \033[0;37m
  NO_COLOR    = \033[m
endif

include local.mk

CFLAGS_OPT=-O3 -std=c++11 -I$(EIGEN_INCLUDE_DIR) -I$(GSL_INCLUDE_DIR) -DOPT_BUILD
CFLAGS_DEBUG=-g -Wall -O0 -std=c++11 -I$(EIGEN_INCLUDE_DIR) -I$(GSL_INCLUDE_DIR) -DDEBUG_BUILD
LFLAGS=-L$(GSL_LIB_DIR) $(LOCAL_LFLAGS)
LLIBS=-lgsl -lgslcblas 

SOURCES=$(wildcard *.hh *.cc)
ifeq ($(MAIN),)
  MAIN=driver.x
endif
OBJS=$(patsubst %.cc,%.o,$(wildcard *.cc))


ifeq ($(DEBUG),True)
  CFLAGS=$(CFLAGS_DEBUG)
else
  CFLAGS=$(CFLAGS_OPT)
endif

ifeq ($(USE_MPI),True)
  CXX=$(MPICXX)
  CFLAGS:=$(CFLAGS) -DUSE_MPI
  $(info Compiling in MPI parallel mode. Compiler is $(CXX))
else
  CXX=g++
  $(info Compiling in sequential mode. Compiler is $(CXX))
endif

all: $(MAIN)

# Sort out dependencies, see
# http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/
DEPDIR := .dependencies
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
POSTCOMPILE = @mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

%.o: %.cc %.hh
	@printf "%b" "$(COM_COLOR)$(COM_MESSAGE) $(OBJ_COLOR)$(@)$(NO_COLOR)\n";
	@$(CXX) $(DEPFLAGS) $(CFLAGS) -c -o $@ $<
	@$(POSTCOMPILE)

%.o: %.cc
	@printf "%b" "$(COM_COLOR)$(COM_MESSAGE) $(OBJ_COLOR)$(@)$(NO_COLOR)\n";
	@$(CXX) $(DEPFLAGS) $(CFLAGS) -c -o $@ $<
	@$(POSTCOMPILE)

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(basename $(SOURCES))))

$(MAIN): $(OBJS) config.h
	@printf "%b" "$(COM_COLOR)$(LINK_MESSAGE) $(OBJ_COLOR)$(@)$(NO_COLOR)\n";
	@$(CXX) $(LFLAGS) -o $(MAIN) $(OBJS) $(LLIBS)

doc: $(SOURCES) config.h
	doxygen doxygen.cfg

.phony: clean
clean:
	rm -rf *~ *.o $(DEPDIR) $(MAIN)

.phony: cleandoc
cleandoc:
	rm -rf doc/
