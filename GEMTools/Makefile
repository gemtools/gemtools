#============================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 04/06/2012
# DESCRIPTION: // TODO
#============================================================================

GEM_TOOLS=gem-tools-examples gem-map-2-fastq
SRCS=gem_tools.c gt_input_file.c gt_buffered_map_input.c  \
     gt_template.c gt_alignment.c gt_map.c gt_misms.c \
     gt_vector.c gt_error.c gt_commons.c \
     gt_buffered_output_file.c gt_output_map.c gt_generic_printer.c gt_output_buffer.c
     # gt_output_sam.c 
CC=gcc
AR=ar
ARCH_FLAGS=-D__LINUX__ -msse3 -mssse3 -msse4.2 -fPIC
SUPPRESS_CHECKS=-DNDEBUG -DGT_NO_CONSISTENCY_CHECKS
OPTIMIZTION_FLAGS=-O4 -fomit-frame-pointer -ftree-vectorize $(SUPPRESS_CHECKS)
DEBUG_FLAGS=-g -ggdb3
LIBS=-lpthread

OBJS=$(SRCS:.c=.o)

all: gt_library $(GEM_TOOLS)

gt_library: $(OBJS)
	$(AR) -rcs gem_tools.a gem_tools.o

$(GEM_TOOLS): $(OBJS)
	$(CC) $(ARCH_FLAGS) $(OPTIMIZTION_FLAGS) -o $@ $@.c $(OBJS) $(LIBS)

.c.o:
	$(CC) $(ARCH_FLAGS) $(OPTIMIZTION_FLAGS) -c $<

clean:
	@rm -f $(GEM_TOOLS) *.o *.a core
