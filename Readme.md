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

Note that the python module is also distributed through pypi and can easily be
installed with either pip or easy_install:

    pip install gemtools
    easy_install gemtools

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
bugtracker](algorithms.cnag.cat/mantis).

Deveopment
=====================

In order to develop on the GEM-Tools library, node that the repository contains
a [Vagrant](http://vagrantup.com/) files to be able to quickly setup a virtual 
machine to test and run GEM on you local non-linux-64bit machine.


Changes
=====================

Version 1.4

    [GT-22] - Add the updated gem-2-sam converter
    [GT-23] - Initial mapping error
    [GT-28] - Splitmapper does not contain the trim parameter to do inline trimming
    [GT-19] - Update gem binaries
    [GT-29] - Create a filter for concatenating sets of files

