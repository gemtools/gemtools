
all:
	$(MAKE) -C GEMTools
	python setup.py build

install: all
	python setup.py install

install-user: all
	python setup.py install-user

test: all
	-$(MAKE) -C GEMTools check 
	python setup.py nosetests

clean:
	$(MAKE) -C GEMTools clean
	@rm -Rf build dist
	@rm -Rf python/Gem.egg-info


