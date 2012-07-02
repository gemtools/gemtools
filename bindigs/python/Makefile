#============================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 04/06/2012
# DESCRIPTION: // TODO
#============================================================================

GEM_TOOLS=../../GEMTools
PY_HOME=/apps/PYTHON/2.7.2
CC=gcc
FLAGS=-fPIC
PY_INCLUDE=-I/apps/PYTHON/2.7.2/include/python2.7 -I/apps/PYTHON/2.7.2/lib/python2.7/config
INCLUDES=$(PY_INCLUDE) -I$(GEM_TOOLS)
ARCH_FLAGS=-D__LINUX__ -msse3 -mssse3 -msse4.2 

SRC=py_iterator.c py_mismatch.c py_template_iterator.c gempymodule.c py_alignment.c py_map.c \
	py_template.c

OBJS=$(SRCS:.c=.o)

all: gempy

gempy: $(OBJS)
	ld -shared $(GEM_TOOLS)/*.o $(OBJS) -o gempy.so

.c.o:
	$(CC) $(FLAGS) $(ARCH_FLAGS) $(INCLUDES) -c $<

test: gempy
	python -m unittest discover  tests '*_tests.py'

clean:
	@rm -f *.o *.a *.so gemtools.py *.pyc gem/*.pyc tests/*.pyc
	@rm -Rf testresults
