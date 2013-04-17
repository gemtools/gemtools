GEM-Tools
===================
GEM-Tools is a C API and a Python module to support and simplify usage of the
[GEM Mapper](http://algorithms.cnag.cat/wiki/The_GEM_library).

The C API currently supports fast and exhaustive support for parsing GEM .map
files and other mapping formats.

The Python *gem* module allows to integrate the GEM Mapper in python scripts
and simplifies development of mapping pipelines.

In addition to the library functionality, GEM-Tools also provides a command
line tool `gemtools`, that you can use to start the RNASeq pipeline and various
other tools, like the indexer and the statistics module.

Licensing
---------
The Python modules also distributes current GEM binaries that are compatible
with the library code. Note that the GEM binaries are distributed under the
[*GEM Non commercial binary license*](http://algorithms.cnag.cat/wiki/GEM:Non_commercial_binary_license)
while the GEM-Tools is licensed under GPL.

Installation
==================
The GEM-Tools library is distributed in three different flavors. You can get a
statically compiled binary bundle if you are just interested in teh command
line tools. If you would like to explore the library functionality, you can
install the latest release from pypi. If you want to get the latest and
greatest, you can clone the git repository and install the library from source. 

Static distribution
-------------------

The GEM-Tools command line tools are bundled and distributed as a static binary
bundle that includes all dependencies. You do not have to install any other
software to use the binary bundle, but you will not be able to write your own
python scripts. You are restricted to the available command line tools.

*NOTE* that the GEM binaries are compiled with optimizations for i3 processors.
We also provide a bundle for older CPU's that lack certain features. If you are
not sure which package to download, i3 or core2, run the following command on
the machine you want to install GEM-Tools:
	
	python -c 'print ("i3 compatible" if set(["popcnt", "ssse3", "sse4_1", "sse4_2"]).issubset(set([f for sl in map(lambda x:x.split(), filter(lambda x:x.startswith("flags"), open("/proc/cpuinfo").readlines())) for f in sl])) else "No i3 support detected. Please use the core2 bundle")'

This command checks your ``/proc/cpuinfo`` file for the flags *popcnt*, *ssse3*,
*sse4_1* and *sse4_2* are present. If that is the case, the CPU is supported by the
i3 bundle. Otherwise use core2.

The latest release can be downloaded here:

* [GEM-Tools static binary bundle 1.6 for i3](http://barnaserver.com/gemtools/releases/GEMTools-static-i3-1.6.tar.gz)
* [GEM-Tools static binary bundle 1.6 for core2](http://barnaserver.com/gemtools/releases/GEMTools-static-core2-1.6.tar.gz)

Library distribution
--------------------

In order to install the GEM-Tool library on your system and get access to the
command line tools as well as the library functionality, you need to have a
couple of dependencies installed. Unfortunately we can currently not install
those dependencies for you. After all dependencies are installed, you have all
the options to build and install the library in your local machine. These
include:

* Install a released version from pypi
* Install form the git repository
* Build the library bundle
* Build the static binary distribution

Dependencies
------------

The C API needs to have gzlib and bzlib installed with header files.  Both
libraries are used to transparently open compressed files. On a Debian/Ubuntu
system the packages are libbz2-dev and zlib-dev.

Here is an example of how you can install the necessary dependencies to build
the C-library on a Debian/Ubuntu system: 

	sudo apt-get install make gcc libbz2-dev

For the python part, the library will work both with Python 2.6 and Python 2.7.
In order to compile the the C binding, you need to have
[Cython](http://www.cython.org/) installed, as well as the python header files. 

Here is an example of how you can install the necessary dependencies to build
the Python library on a Debian/Ubuntu system: 

	sudo apt-get install make gcc libbz2-dev python-dev


Install from pypi 
-----------------

Gemtools is distributed through pypi and you can install the latest released
version using `pip` or `easy_install`. But please make sure you have all the
dependencies installed first.

As root:

    pip install gemtools

As non-root user:
    
    pip install gemtools --user


Install from github
-------------------

You can install GEM-Tool from the github repository using the distributed 
setup.py script. 

Clone the repository
	
	git clone https://github.com/gemtools/gemtool

Change into the gemtools folder and install the python library

	cd gemtools

As root:
    
    python setup.py install

As non-root user:
    
	python setup.py install --user 

This will install into a $HOME/.local folder and you might want to add
$HOME/.local/bin to your path. For example, add

    export PATH=$PATH:$HOME/.local/bin

to your ~/.bashrc to make the change permanent.


Build the library package
-------------------------

From the repository, you can build the python library package. This you can use
to moved to a dedicated folder and can be managed by your prefered module
manager easily. 

In order to build the library package, clone the the repository and call

	make package

This will create a *dist/* folder that contains two tarballs, one for i3, and
one for core2. In addition you will see the unpacked folders. They are
structured as follows:

	bin/ -- all the executables
	lib/ -- symlink to lib64
	lib64/ -- contains the gemtools c library as limgemtools.a
	lib64/python<version>/site-packages -- the python libraries for your python version (i.e. 2.6 or 2.7)
	include/ -- the C API header files

Say you moved the package to `/opt/gemtools` and you run python2.6. You can
activate the package by exporting the following variables to your environment:

	export PATH=$PATH:/opt/gemtools/bin
	export PYTHONPATH=$PYTHONPATH:/opt/gemtools/lib64/python2.6/site-packages

This will make all the executables available in your path and put the python
library into the python search path so you can leverage it from your scripts.


Build the static binary bundle
------------------------------

If you are not interested in using any python library functions from your 
script, you can build portable static binary package by calling the `dist` 
target in the Makefile.

	make dist

This will create a *dist/* folder that contains two tarballs, one for i3, and
one for core2. In addition you will see the unpacked folders. They are
structured as follows:

	bin/ -- all the executables
	lib/ -- symlink to lib64
	lib64/ -- contains the gemtools c library as limgemtools.a
	include/ -- the C API header files

Say you moved the package to `/opt/gemtools`. You can activate the package by
exporting the following variables to your environment:

	export PATH=$PATH:/opt/gemtools/bin

	
Bugs and feature requests
=====================
Please feel free to use the Github bug tracker to report issues and feature
requests that you find in the GEM-Tools library. 

If you run into problems with GEM and the included GEM binaries, please use 
the [GEM bugtracker](http://algorithms.cnag.cat/mantis).



Change log
=====================
	1.6.1
    - GEMTools SAM conversion doesnow include the header again
    - Pair detection fixed for reads with space in the id and not casava
    - Pipeline -o paramter does create the output folder now automatically
    - samtools sorter max memory is now properly set for the non-threads version
    - Added option to create sam-compact format in the pipeline

    1.6
	- Parallel file conversion for pipeline
	- Cython C bindings for the parser library
	- JSON export of the statistics
	- Index and transcript index support from gemtools
	- Restartable pipeline options
	- gem-2-sam now writes the XS field for cufflinks support
	- Added the new rna seq pipeline
	- Added stats and report tools to the gemtools cli tool
    - Added new bundle system and detection for i3 vs core2
    - Added optional appending NH field to SAM entries
    - Added *extra* options to mapper/splitmapper/pairalign to enable 
      passing custom command line options to the tools.
    - Modify C::counters_functions to use zero-based indexes

    1.5
    - Fixed issues with gem-rna-mapper output and validation
    - Moved validation to be a default step after split-mapping
    - Added stats filter
    - Fixed issue with gem-2-sam binary

