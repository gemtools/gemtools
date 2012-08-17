#============================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 04/06/2012
# DESCRIPTION: // TODO
#============================================================================

GEM_TOOLS=gem-source #gem-tools-examples
SRCS=gem_tools.c \
     gt_input_file.c gt_buffered_input_file.c \
       gt_input_map_parser.c gt_input_sam_parser.c gt_input_fastq_parser.c \
     gt_buffered_output_file.c gt_output_buffer.c \
       gt_output_map.c gt_generic_printer.c \
     gt_template.c gt_alignment.c gt_map.c gt_misms.c \
     gt_template_handling.c gt_alignment_handling.c \
     gt_hash.c gt_vector.c gt_error.c gt_commons.c 
CC=gcc
AR=ar
LIBS=-lpthread

# Compilation flags
GENERAL_FLAGS=-fPIC
ARCH_FLAGS=-D__LINUX__
SUPPRESS_CHECKS=-DNDEBUG -DGT_NO_CONSISTENCY_CHECKS
OPTIMIZTION_FLAGS=-O4 $(SUPPRESS_CHECKS) # -fomit-frame-pointer -ftree-vectorize
ARCH_FLAGS_OPTIMIZTION_FLAGS=-msse3 -mssse3 -msse4.2
DEBUG_FLAGS=-g -ggdb3

OBJS=$(SRCS:.c=.o)

all: GEM_TOOLS_FLAGS=-O4 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(SUPPRESS_CHECKS) $(ARCH_FLAGS_OPTIMIZTION_FLAGS)
all: gt_library $(GEM_TOOLS)

debug: GEM_TOOLS_FLAGS=-O0 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(DEBUG_FLAGS) 
debug: gt_library $(GEM_TOOLS)

gt_library: $(OBJS)
	$(AR) -rcs gem_tools.a gem_tools.o

$(GEM_TOOLS): $(OBJS)
	$(CC) $(GEM_TOOLS_FLAGS) -o $@ $@.c $(OBJS) $(LIBS)

.c.o:
	$(CC) $(GEM_TOOLS_FLAGS) -c $<

clean:
	@rm -f $(GEM_TOOLS) *.o *.a core
