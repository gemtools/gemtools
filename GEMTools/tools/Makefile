#==================================================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 02/10/2012
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Builds the GT-tools
#==================================================================================================

# Definitions
ROOT_PATH=..
include ../Makefile.mk

GEM_TOOLS=gt.construct gt.stats gt.filter gt.mapset gt.map2sam align_stats sel_reads gt.gtfcount gt.region

GEM_TOOLS_SRC=$(addsuffix .c, $(GEM_TOOLS))
GEM_TOOLS_BIN=$(addprefix $(FOLDER_BIN)/, $(GEM_TOOLS))
LIBS:=-lgemtools -lpthread -lm
ifeq ($(HAVE_OPENMP),1)
LIBS:=$(LIBS) -fopenmp
endif
ifeq ($(HAVE_ZLIB),1)
LIBS:=$(LIBS) -lz
endif
ifeq ($(HAVE_BZLIB),1)
LIBS:=$(LIBS) -lbz2
endif

all: GEM_TOOLS_FLAGS=-O4 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(SUPPRESS_CHECKS) $(OPTIMIZTION_FLAGS) $(ARCH_FLAGS_OPTIMIZTION_FLAGS)
all: $(GEM_TOOLS_BIN)

static: GEM_TOOLS_FLAGS=-O4 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(SUPPRESS_CHECKS) $(OPTIMIZTION_FLAGS) $(ARCH_FLAGS_OPTIMIZTION_FLAGS) -static
static: $(GEM_TOOLS_BIN)

debug: GEM_TOOLS_FLAGS=-O0 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(DEBUG_FLAGS)
debug: $(GEM_TOOLS_BIN)

$(FOLDER_BIN)/gt.stats: EXTRA_LIBS=-ljson
$(GEM_TOOLS_BIN): $(FOLDER_LIB)/libgemtools.a $(GEM_TOOLS_SRC)
	$(CC) $(GEM_TOOLS_FLAGS) -o $@ $(notdir $@).c $(LIB_PATH_FLAGS) $(INCLUDE_FLAGS) $(LIBS) $(EXTRA_LIBS)
