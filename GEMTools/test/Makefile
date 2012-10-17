#==================================================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 02/10/2012
# DESCRIPTION: Builds and runs all the tests
#==================================================================================================

# Definitions
ROOT_PATH=..
include ../Makefile.mk

GT_TESTS=gt_test_commons gt_test_core_structures gt_test_parsers
GT_TESTS_FLAGS=$(ARCH_FLAGS) $(DEBUG_FLAGS)
GT_COVERAGE_FLAGS=-g -Wall -fprofile-arcs -ftest-coverage $(GT_TESTS_FLAGS)

LIBS=-lpthread -lgemtools -lcheck

all: $(GT_TESTS) end_banner

coverage: $(GT_COVERAGE) end_banner

setup:
	@mkdir -p build

$(GT_TESTS): clean setup
	$(CC) $(GT_TESTS_FLAGS) $(INCLUDE_FLAGS) $(LIB_PATH_FLAGS) -o ./build/$@ $@.c $(LIBS)
	@echo "=======================================================================>>"
	@echo "==>> Testing " $@
	@echo "=======================================================================>>"
	./build/$@
	
$(GT_COVERAGE): clean setup
	$(CC) $(GT_COVERAGE_FLAGS) $(INCLUDE_FLAGS) $(LIB_PATH_FLAGS) -o ./build/$@ $@.c $(LIBS)
	./build/$@
	
end_banner: 	
	@echo "=======================================================================>>"
	
clean:
	@rm -rf $(FOLDER_TEST_BUILD)/*
	@rm -rf build
