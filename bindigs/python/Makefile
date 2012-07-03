#============================================================================
# PROJECT: GEM-Tools library
# FILE: Makefile
# DATE: 04/06/2012
# DESCRIPTION: // TODO
#============================================================================


all: gempy

gempy: $(OBJS)
	$(MAKE) --directory=../../GEMTools
	rm -Rf build/
	python setup.py build

test: gempy
	python setup.py nosetests 

clean:
	@rm -f *.o *.a *.so gemtools.py *.pyc gem/*.pyc tests/*.pyc
	@rm -Rf testresults build Gem.egg-info
	@rm -f gem/*.pyc gem/*.so
	@rm -f test/*.pyc 
