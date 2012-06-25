#============================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 04/06/2012
# DESCRIPTION: // TODO
#============================================================================

GEM_TOOLS=../../GEMTools
PY_HOME=/apps/PYTHON/2.7.2
CC=gcc
FLAGS=-lpthread -lc -fPIC
PY_INCLUDE=-I/apps/PYTHON/2.7.2/include/python2.7 -I/apps/PYTHON/2.7.2/lib/python2.7/config
INCLUDES=$(PY_INCLUDE) -I$(GEM_TOOLS)
all: gemtools.py

gemtools.py: gem_tools_wrap.o
	ld -shared $(GEM_TOOLS)/*.o gem_tools_wrap.o -o _gemtools.so

gem_tools_wrap.c: gem_tools.i
	swig -Wall -I$(GEM_TOOLS) -python gem_tools.i

gem_tools_wrap.o: gem_tools_wrap.c
	$(CC) $(FLAGS) $(INCLUDES) -c gem_tools_wrap.c

test: gemtools.py
	python -m unittest discover  tests '*_tests.py'

clean:
	@rm -f *.o *.a *.c *.so gemtools.py *.pyc gem/*.pyc
