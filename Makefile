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

BUILD_DIR:=build
SOURCE_DIR:=src

CFLAGS_OPT=-O3 -std=c++11 -I$(shell pwd)/$(SOURCE_DIR) -I$(EIGEN_INCLUDE_DIR) -I$(GSL_INCLUDE_DIR) -DNDEBUG -DOPT_BUILD
CFLAGS_DEBUG=-g -Wall -O0 -std=c++11 -I$(shell pwd)/$(SOURCE_DIR) -I$(EIGEN_INCLUDE_DIR) -I$(GSL_INCLUDE_DIR) -DDEBUG_BUILD
LFLAGS=-L$(GSL_LIB_DIR) $(LOCAL_LFLAGS)
LLIBS=-lgsl -lgslcblas 

HEADERS=$(shell find $(SOURCE_DIR) -name "*.hh")
SOURCES=$(shell find $(SOURCE_DIR) -name "*.cc")

ifeq ($(MAIN_QM),)
  MAIN_QM=driver_qm
endif

MAIN_QFT=driver_qft
MAIN_TESTDIST=test_distribution
MAIN_TESTSCHWINGERFILLIN=test_schwinger_fillin_distribution
MAIN_FASTBESSEL=test_fastbessel

MAINS=\
	driver_qm.o \
	driver_qft.o \
	test_distribution.o \
	test_schwinger_fillin_distribution.o \
	test_fastbessel.o
OBJS:=$(filter-out $(patsubst %,$(SOURCE_DIR)/%,$(MAINS)),$(patsubst %.cc,%.o,$(SOURCES)))
OBJS:=$(patsubst %,$(BUILD_DIR)/%,$(OBJS))

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
  $(info Compiling in sequential mode. Compiler is $(CXX))
endif

all: $(MAIN_QM) $(MAIN_QFT) $(MAIN_TESTDIST) $(MAIN_FASTBESSEL) $(MAIN_TESTSCHWINGERFILLIN)

# Sort out dependencies, see
# http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/
DEPDIR := $(BUILD_DIR)/dependencies
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
PRECOMPILE = @mkdir -p $(DEPDIR)/$(dir $<); mkdir -p $(BUILD_DIR)/$(dir $<)
POSTCOMPILE = @mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

$(BUILD_DIR)/%.o: %.cc %.hh
	@printf "%b" "$(COM_COLOR)$(COM_MESSAGE) $(OBJ_COLOR)$(@)$(NO_COLOR)\n";
	$(PRECOMPILE)
	$(CXX) $(DEPFLAGS) $(CFLAGS) -c -o $@ $<
	@$(POSTCOMPILE)

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.cc $(OBJS)
	@printf "%b" "$(COM_COLOR)$(COM_MESSAGE) $(OBJ_COLOR)$(@)$(NO_COLOR)\n";
	@$(PRECOMPILE)
	$(CXX) $(DEPFLAGS) $(CFLAGS) -c -o $@ $<
	@$(POSTCOMPILE)

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(basename $(HEADERS))))

# main programmes
# --- QM ---
$(MAIN_QM): $(BUILD_DIR)/driver_qm.o $(OBJS) $(SOURCE_DIR)/config.h
	@printf "%b" "$(COM_COLOR)$(LINK_MESSAGE) $(OBJ_COLOR)$(@)$(NO_COLOR)\n";
	@$(CXX) $(LFLAGS) -o $(BUILD_DIR)/$(MAIN_QM) $(OBJS) $(BUILD_DIR)/driver_qm.o $(LLIBS)

# --- QFT (Schwinger- and non-linear sigma- model) ---
$(MAIN_QFT): $(BUILD_DIR)/driver_qft.o $(OBJS) $(SOURCE_DIR)/config.h
	@printf "%b" "$(COM_COLOR)$(LINK_MESSAGE) $(OBJ_COLOR)$(@)$(NO_COLOR)\n";
	@$(CXX) $(LFLAGS) -o $(BUILD_DIR)/$(MAIN_QFT) $(OBJS) $(BUILD_DIR)/driver_qft.o $(LLIBS)

# --- test distributions ---
$(MAIN_TESTDIST): $(BUILD_DIR)/test_distribution.o $(OBJS) $(SOURCE_DIR)/config.h
	@printf "%b" "$(COM_COLOR)$(LINK_MESSAGE) $(OBJ_COLOR)$(@)$(NO_COLOR)\n";
	@$(CXX) $(LFLAGS) -o $(BUILD_DIR)/test_distribution $(OBJS) $(BUILD_DIR)/test_distribution.o $(LLIBS)

# --- test Schwinger fillin ---
$(MAIN_TESTSCHWINGERFILLIN): $(BUILD_DIR)/test_schwinger_fillin_distribution.o $(OBJS) $(SOURCE_DIR)/config.h
	@printf "%b" "$(COM_COLOR)$(LINK_MESSAGE) $(OBJ_COLOR)$(@)$(NO_COLOR)\n";
	@$(CXX) $(LFLAGS) -o $(BUILD_DIR)/test_schwinger_fillin_distribution $(OBJS) $(BUILD_DIR)/test_schwinger_fillin_distribution.o $(LLIBS)

# --- test fast Bessel evaluation ---
$(MAIN_FASTBESSEL): $(BUILD_DIR)/test_fastbessel.o $(OBJS) $(SOURCE_DIR)/config.h
	@printf "%b" "$(COM_COLOR)$(LINK_MESSAGE) $(OBJ_COLOR)$(@)$(NO_COLOR)\n";
	@$(CXX) $(LFLAGS) -o $(BUILD_DIR)/test_fastbessel $(OBJS) $(BUILD_DIR)/test_fastbessel.o $(LLIBS)

doc: $(SOURCES) ./src/config.h
	doxygen doxygen.cfg

.phony: clean
clean:
	rm -rf *~ $(BUILD_DIR)/*

.phony: cleandoc
cleandoc:
	rm -rf doc/
