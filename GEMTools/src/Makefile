#==================================================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 02/10/2012
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Builds all the sources of the GT-library, and the library itself (.a)
#==================================================================================================

# Definitions
ROOT_PATH=..
include ../Makefile.mk

SRCS=gem_tools.c \
     gt_commons.c gt_error.c gt_mm.c gt_fm.c gt_profiler.c \
     gt_ihash.c gt_shash.c gt_vector.c gt_string.c \
     gt_attributes.c gt_dna_string.c gt_dna_read.c gt_compact_dna_string.c \
     gt_template.c gt_alignment.c gt_map.c gt_misms.c \
     gt_template_utils.c gt_alignment_utils.c gt_counters_utils.c \
     gt_map_metrics.c gt_map_align.c gt_map_score.c \
     gt_sequence_archive.c gt_segmented_sequence.c \
     gt_input_file.c gt_buffered_input_file.c \
     gt_input_parser.c gt_input_map_parser.c gt_input_fasta_parser.c gt_input_generic_parser.c \
     gt_input_sam_parser.c gt_sam_attributes.c \
     gt_buffered_output_file.c gt_output_file.c gt_generic_printer.c gt_output_buffer.c \
     gt_output_printer.c gt_output_map.c gt_output_fasta.c gt_output_sam.c gt_output_generic_printer.c \
     gt_stats.c gt_gemIdx_loader.c gt_gtf.c
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.c=.o))
GT_LIB=$(FOLDER_LIB)/libgemtools.a

all: GEM_TOOLS_FLAGS=$(GENERAL_FLAGS) $(ARCH_FLAGS) $(SUPPRESS_CHECKS) $(OPTIMIZTION_FLAGS) $(ARCH_FLAGS_OPTIMIZTION_FLAGS)
all: $(GT_LIB)

debug: GEM_TOOLS_FLAGS=-O0 $(GENERAL_FLAGS) $(ARCH_FLAGS) $(DEBUG_FLAGS)
debug: $(GT_LIB)

$(GT_LIB): $(OBJS)
	$(AR) -rcs $(GT_LIB) $(FOLDER_BUILD)/*.o

$(FOLDER_BUILD)/gt_fm.o : gt_fm.c
	$(CC) $(GEM_TOOLS_FLAGS) $(INCLUDE_FLAGS) -c $< -o $@ -fopenmp
$(FOLDER_BUILD)/gt_mm.o : gt_mm.c
	$(CC) $(GEM_TOOLS_FLAGS) $(INCLUDE_FLAGS) -c $< -o $@ -fopenmp

$(FOLDER_BUILD)/%.o : %.c
	$(CC) $(GEM_TOOLS_FLAGS) $(INCLUDE_FLAGS) -c $< -o $@
