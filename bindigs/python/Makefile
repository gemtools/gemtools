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
all: gemtools.py

gemtools.py: gem_tools_wrap.o gempymodule.o
	ld -shared $(GEM_TOOLS)/*.o gem_tools_wrap.o -o _gemtools.so
	ld -shared $(GEM_TOOLS)/*.o gempymodule.o -o gempy.so

gem_tools_wrap.c: gem_tools.i
	swig -Wall -I$(GEM_TOOLS) -builtin -nocontract -nodefaultctor -nodefaultdtor -fastdispatch -python gem_tools.i

gem_tools_wrap.o: gem_tools_wrap.c
	$(CC) $(FLAGS) $(INCLUDES) -c gem_tools_wrap.c

gempymodule.o:
	$(CC) $(FLAGS) $(ARCH_FLAGS) $(INCLUDES) -c gempymodule.c

test: gemtools.py
	python -m unittest discover  tests '*_tests.py'

clean:
	@rm -f *.o *.a *_wrap.c *.so gemtools.py *.pyc gem/*.pyc tests/*.pyc
	@rm -Rf testresults
