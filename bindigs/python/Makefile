#============================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 04/06/2012
# DESCRIPTION: // TODO
#============================================================================


all: gempy

gempy: $(OBJS)
	rm -Rf build/
	python setup.py build

test: gempy
	python -m unittest discover  tests '*_tests.py'

clean:
	@rm -f *.o *.a *.so gemtools.py *.pyc gem/*.pyc tests/*.pyc
	@rm -Rf testresults
