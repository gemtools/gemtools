GEM-Tools
===================
GEM-Tools is a C API and a Python module to support and simplify usage of the
[GEM Mapper](http://algorithms.cnag.cat/wiki/The_GEM_library).

The C API curtly supports fast and exhaustive support for parsing GEM .map
files and other mapping formats.

The Python *gem* module allows to integrate the GEM Mapper in python scripts
and simplifies development of mapping pipelines.

The Python modules also distributes current GEM binaries that are compatible
with the library code. Note that the GEM binaries are distributed under the
[*GEM Non commercial binary license*](http://algorithms.cnag.cat/wiki/GEM:Non_commercial_binary_license)
while the GEM-Tools is licensed under GPL.

Installation
==================
Installing the gemtools python module can be done in different ways. The easiest
way is to install using *pip*. If you do not have *pip* installed or you
do not want to use it, clone the githup repository and insall from source.

Dependencies
----------------------------

For the C API you need to have gzlib and bzlib installed with header files.
Both libraries are used to transparently open compressed files.

For python, the library will work both with Python 2.6 and Python 2.7, but for
2.6 you have ti install the argparse library. An easy way to install gemtools
is to work with virtualenv to set up a custom environment. See the virtualenv
section of the installation instructions. 

Install the latest release
----------------------------
Gemtools is distributed through pypi and you can install the lates released
vesion using pip or easy install

As root:

    pip install gemtools

As non-root user:
    
    pip install gemtools --user

Install from github
---------------------------
You can use pip to install directly from github:

As root:
    
    pip install git+http://github.com/gemtools/gemtools

As non-root user:
    
    pip install git+http://github.com/gemtools/gemtools --user

Install without pip
-------------------
If you do not have *pip* or you cloned the github repository already,
you can simply run:

As root:
    
    python setup.py install

As non-root user:
    
    python setup.py install --user
    
This will install into a $HOME/.local folder and you might want to add
$HOME/.local/bin to your path. For example, add

    export PATH=$PATH:$HOME/.local/bin

to your ~/.bashrc to make the change permanent.

You can also install into a custom folder. For example, to install into a
$HOME/usr/gemtools folder, do the following:

    PYTHONPATH=$PYTHONPATH:$HOME/usr/gemtools/lib python setup.py install --home $HOME/usr/gemtools    

This will also create a $HOME/usr/gemtools/bin folder that you might want to
add to your path. In addtition you have to adopt your global PYTHONPATH and
include the new location. Put this into your .bashrc:

    export PATH=$PATH:$HOME/usr/gemtools/bin
    export PYTHONPATH=$PYTHONPATH:$HOME/usr/gemtools/lib

Install as non-root user
--------------------------
There are a couple of different ways to install the gemtools module as non
root user, but the simplest is to append the --user option to the installtion
method of choice. This will install gemtools (and all other python modules you
install that way) into the $HOME/.local folder. This folder is already in your
default python paath and you do not have to modify any environment variables.

Install into a custom environment with virtualenv
-------------------------------------------------
Virtualenv is an excellent tool to set up independend python environment and
install dependencies without root previleges and without introducing any
conflicts.

If you do not have virtualenv installed, simply download the
[https://raw.github.com/pypa/virtualenv/master/virtualenv.py](virtualenv.py])
script, that is all you need to get started. Then create the folder where you
want to install the library to and create a virtual environment there. Here we
use $HOME/usr/gemtools as an example:

    mkdir -p $HOME/usr/gemtools
    python virtualenv.py $HOME/usr/gemtools

Verify the installation
-------------------------
To quickly check if gemtools is installed, run:
    
    python -m gem.__main__

If gemtools is installed properly, it will print the paths to the bundled binaries.

Examples
======================
We host [https://github.com/gemtools/gemtools-examples](another GIT reposirory) with 
pipeline examples that are self contained and can be run immediately to explore
the basic functionallity of the library. 

Bugs and feature requests
=====================
Please feel free to use the Github bug tracker to report issues and feature
requests that you find in the GEM-Tools library. If you run into problems with
GEM and the included GEM binaries, please use the [GEM
bugtracker](http://algorithms.cnag.cat/mantis).

Development
=====================
In order to develop on the GEM-Tools library, node that the repository contains
a [Vagrant](http://vagrantup.com/) files to be able to quickly setup a virtual 
machine to test and run GEM on you local non-linux-64bit machine.

Changelog
=====================

    1.6
    - Added new bundle system and detection for i3 vs core2
    - Added optional appending NH field to SAM entries
    - Added *extra* options to mapper/splitmapper/pairalign to enable 
      passing custom command line options to the tools. 

    1.5
    - Fixed issues with gem-rna-mapper output and validation
    - Moved validation to be a default step after split-mapping
    - Added stats filter
    - Fixed issue with gem-2-sam binary

