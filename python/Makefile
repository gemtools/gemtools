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
	cp build/lib.*/gem/gemtools.so gem/
	nosetests 

install: gempy
	python setup.py install

install-user: gempy
	python setup.py install --user

install-devel-cnag: gempy
	PYTHONPATH=/project/devel/apps/gemtools/Devel/lib64/python/ python setup.py install --home /project/devel/apps/gemtools/Devel

clean:
	@rm -f *.o *.a *.so gemtools.py *.pyc gem/*.pyc tests/*.pyc
	@rm -Rf testresults build Gem.egg-info
	@rm -f gem/*.pyc gem/*.so
	@rm -f test/*.pyc 
	@rm -Rf dist 
	@rm -Rf gemtools-env 
	@rm -f nosetests.xml 
