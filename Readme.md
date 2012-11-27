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

If you do not have *pip* or you cloned the github repository already,
you can simply run:

As root:
    
    python setup.py install

As non-root user:
    
    python setup.py install --user

Install as non-root user
--------------------------
There are a couple of different ways to install the gemtools module as non
root user, but the simplest is to append the --user option to the installtion
method of choice. This will install gemtools (and all other python modules you
install that way) into the $HOME/.local folder. This folder is already in your
default python paath and you do not have to modify any environment variables.

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

    1.5
    - Fixed issues with gem-rna-mapper output and validation
    - Moved validation to be a default step after split-mapping
    - Added stats filter
    - Fixed issue with gem-2-sam binary
