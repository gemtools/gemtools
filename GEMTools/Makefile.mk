#==================================================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile.mk
# DATE: 02/10/2012
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Makefile definitions' file
#==================================================================================================

# Utilities
CC=gcc
AR=ar

# Folders
FOLDER_BIN=$(ROOT_PATH)/bin
FOLDER_BUILD=$(ROOT_PATH)/build
FOLDER_DATASETS=$(ROOT_PATH)/datasets
FOLDER_INCLUDE=$(ROOT_PATH)/include
FOLDER_LIB=$(ROOT_PATH)/lib
FOLDER_RESOURCES=$(ROOT_PATH)/resources
FOLDER_RESOURCES_LIB=$(ROOT_PATH)/resources/lib
FOLDER_RESOURCES_INCLUDE=$(ROOT_PATH)/resources/include
FOLDER_SOURCE=$(ROOT_PATH)/src
FOLDER_TEST=$(ROOT_PATH)/test
FOLDER_TOOLS=$(ROOT_PATH)/tools

# Flags
ARCH_FLAGS = -D__LINUX__
HAVE_ZLIB = 1
HAVE_BZLIB = 0
HAVE_OPENMP = 1

GENERAL_FLAGS=-std=gnu89 -fPIC -Wall
ifeq ($(HAVE_ZLIB),1)
GENERAL_FLAGS:=$(GENERAL_FLAGS) -DHAVE_ZLIB
endif
ifeq ($(HAVE_BZLIB),1)
GENERAL_FLAGS:=$(GENERAL_FLAGS) -DHAVE_BZLIB
endif
ifeq ($(HAVE_OPENMP),1)
GENERAL_FLAGS:=$(GENERAL_FLAGS) -DHAVE_OPENMP
endif

OPTIMIZTION_FLAGS=-O4 # -fomit-frame-pointer -ftree-vectorize
ARCH_FLAGS_OPTIMIZTION_FLAGS= # -msse3 -mssse3 -msse4.2

INCLUDE_FLAGS=-I$(FOLDER_INCLUDE) -I$(FOLDER_RESOURCES_INCLUDE)
LIB_PATH_FLAGS=-L$(FOLDER_LIB) -L$(FOLDER_RESOURCES_LIB)

SUPPRESS_CHECKS=-DNDEBUG -DGT_NO_CONSISTENCY_CHECKS
DEBUG_FLAGS=-g -ggdb3 -rdynamic -DGT_DEBUG

PLATFORM=$(shell uname)
ifeq ($(PLATFORM),Darwin)
FOLDER_RESOURCES_LIB=$(ROOT_PATH)/resources/lib/mac
endif

