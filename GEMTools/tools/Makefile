#==================================================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 02/10/2012
# DESCRIPTION: Builds the GT-tools
#==================================================================================================

# Definitions
ROOT_PATH=..
include ../Makefile.mk

GEM_TOOLS=gem-tools-examples gt.stats 
LIBS=-lgemtools -lpthread -fopenmp

all: GEM_TOOLS_FLAGS=-O4 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(SUPPRESS_CHECKS) $(OPTIMIZTION_FLAGS) $(ARCH_FLAGS_OPTIMIZTION_FLAGS)
all: $(GEM_TOOLS)

debug: GEM_TOOLS_FLAGS=-O0 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(DEBUG_FLAGS)
debug: $(GEM_TOOLS)

$(GEM_TOOLS):
	$(CC) $(GEM_TOOLS_FLAGS) -o $(FOLDER_BIN)/$@ $@.c $(LIB_PATH_FLAGS) $(INCLUDE_FLAGS) $(LIBS) 
