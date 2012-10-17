#============================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 02/10/2012
# DESCRIPTION: Top level makefile
#============================================================================

# Definitions
ROOT_PATH=$(CURDIR)
include Makefile.mk

all: release
	
release:
	@mkdir -p build bin lib
	$(MAKE) --directory=src
	$(MAKE) --directory=tools

debug:
	$(MAKE) --directory=src   debug
	$(MAKE) --directory=tools debug

check:
	$(MAKE) --directory=test

clean:
	@rm -f $(FOLDER_BIN)/* $(FOLDER_BUILD)/*.o $(FOLDER_LIB)/*.a $(FOLDER_TEST_BUILD)/*
	@rm -Rf build bin lib
