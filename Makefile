VERSION			   = "1.7"
VIRTUALENV         = dist-utils/virtualenv.py
VIRTUALENV_NAME    = environment

all:
	$(MAKE) -C GEMTools
	python setup.py build

gemtools:
	$(MAKE) -C GEMTools release

install: all
	python setup.py install

install-user: all
	python setup.py install --user

test: all
	-$(MAKE) -C GEMTools check
	python setup.py nosetests

test-c: all
	-$(MAKE) -C GEMTools check

test-python: all
	python setup.py nosetests --with-xunit

package:
	python setup.py package

# pyinstaller distribution
prep-dist: all
	python setup.py build_ext -i

fetch-pyinstaller:
	test -s dist-utils/pyinstaller-2.0.tar.bz2 || wget -O dist-utils/pyinstaller-2.0.tar.bz2  http://sourceforge.net/projects/pyinstaller/files/2.0/pyinstaller-2.0.tar.bz2/download

install-pyinstaller: fetch-pyinstaller
	test -s dist-utils/pyinstaller-2.0 || tar -C dist-utils -xjvf dist-utils/pyinstaller-2.0.tar.bz2

dist: prep-dist install-pyinstaller
	@mkdir -p dist
	./dist-utils/pyinstaller-2.0/pyinstaller.py -F -p python/ --hidden-import=matplotlib --hidden-import=pylab --hidden-import=numpy --hidden-import=gem --hidden-import=gem.gemtools -o dist/ -y python/gem/commands.py -n gemtools
	python setup.py package_static

clean:
	$(MAKE) -C GEMTools clean
	@rm -Rf build dist
	@rm -Rf python/Gem.egg-info
	@rm -Rf python/Gemtools.egg-info
	@rm -Rf python/gem/*.so
	@rm -Rf python/gem/*.pyc
	@rm -Rf python/src/gemtools.c


devel: venv
venv: $(VIRTUALENV_NAME)/bin/activate
$(VIRTUALENV_NAME)/bin/activate: requirements.txt
	test -d $(VIRTUALENV_NAME)/bin/activate || python $(VIRTUALENV) $(VIRTUALENV_NAME)
	. $(VIRTUALENV_NAME)/bin/activate; pip install -Ur requirements.txt
	. $(VIRTUALENV_NAME)/bin/activate; python setup.py develop
	@echo "----------------------------------------------------------"
	@echo "---  Virtualenv set up in $(VIRTUALENV_NAME)"
	@echo "---  and develop install performced. The command in"
	@echo "---  your path now refer to the current development"
	@echo "---  version"
	@echo "---"
	@echo "---  Activate the environment with:"
	@echo "---"
	@echo "---  . $(VIRTUALENV_NAME)/bin/activate"
	@echo "----------------------------------------------------------"


