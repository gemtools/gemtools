#==================================================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 02/10/2012
# DESCRIPTION: Builds all the sources of the GT-library, and the library itself (.a)
#==================================================================================================

# Definitions
ROOT_PATH=..
include ../Makefile.mk

SRCS=gem_tools.c \
     gt_ihash.c gt_shash.c gt_vector.c gt_string.c gt_error.c gt_commons.c \
     gt_dna_string.c gt_dna_read.c \
     gt_template.c gt_alignment.c gt_map.c gt_misms.c \
     gt_template_utils.c gt_alignment_utils.c \
     gt_input_file.c gt_buffered_input_file.c \
     gt_input_map_parser.c gt_input_sam_parser.c gt_input_fastq_parser.c \
     gt_generic_printer.c gt_output_buffer.c gt_output_map.c \
     gt_buffered_output_file.c gt_output_file.c
#       \
#      gt_map_align.c \
#       gt_compact_dna_string.c gt_sequence_archive.c 
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))
GT_LIB=$(FOLDER_LIB)/libgemtools.a

all: GEM_TOOLS_FLAGS=$(GENERAL_FLAGS) $(ARCH_FLAGS) $(SUPPRESS_CHECKS) $(OPTIMIZTION_FLAGS) $(ARCH_FLAGS_OPTIMIZTION_FLAGS) -Wall
all: $(GT_LIB)

debug: GEM_TOOLS_FLAGS=-O0 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(DEBUG_FLAGS) -Wall
debug: $(GT_LIB)

$(GT_LIB): $(OBJS)
	$(AR) -rcs $(GT_LIB) $(FOLDER_BUILD)/*.o

$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(GEM_TOOLS_FLAGS) $(INCLUDE_FLAGS) -c $< -o $@
