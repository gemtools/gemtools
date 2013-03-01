Welcome to GEMTools's documentation!
====================================

GEMTools is a library and a set of command line tools around the `GEM mapper <http://algorithms.cnag.cat/>`_.
The packages consists of a :ref:`python_api`, a :ref:`c_api` and a set of :ref:`command_line_tools`.

The package also bundles the GEMTools :ref:`RNASeq pipeline<rna_pipeline>` that allows
you to map raw reads to both a genome and a transcriptome in a quick and easy fashion.

Download and Installation
-------------------------

There are three different ways you can install GEMTools.

First, we provide a binary package that aims to be used from the command line.
It bundles the *GEM* binaries and the *gemtools* command and allows you to run the RNASeq pipeline directly
from your bash shell. The package is self contained and no further dependencies are needed.
The *gemtools* command offers an easy way to start the pipeline and various other
tools.

Secondly, you can install GEMTools as a python library. The feature set is exactly
the same, the *GEM* binaries come bundled with the library and the *gemtools* command
is exposed, but in addition, you will be able to use the Python API to write your
own scripts to parse GEM output and to build customized mapping pipelines that
integrate more tightly with your own workflow.
In order to install the GEMTools library, you have to have some dependencies installed on you local machine.

The third way to get GEMTools is to get the source distribution directly from our
`GitHub Repository <https://github.com/gemtools/gemtools>`_. Installing from source
will not only allow you to easily stay up to date with the latest development, but
you will be able to use out :ref:`C API <c_api>` and write very fast, parallelized
programs to analyze the results created with GEM.

Both, the binary bundle as well as the library bundle come in two version. One for
the *core2* and one for *i3*. The right package basically depends on the capabilities
of your CPU. The *core2* package will also work with older CPUs, but some of the
optimizations are not available. If you don't know if your CPU supports the *i3*
package, run the following one liner on the machine you want to run GEM on::

    python -c 'print ("i3 compatible" if set(["popcnt", "ssse3", "sse4_1", "sse4_2"]).issubset(set([f for sl in map(lambda x:x.split(), filter(lambda x:x.startswith("flags"), open("/proc/cpuinfo").readlines())) for f in sl])) else "No i3 support detected. Please use the core2 bundle")'

This command checks your ``/proc/cpuinfo`` file for the flags *popcnt*, *ssse3*,
*sse4_1* and *sse4_2* are present. If that is the case, the CPU is supported by the
i3 bundle. Otherwise use core2.


Downloads binary bundle
***********************

* `GEMTools binary bundle 1.6 for i3 <http://bla>`_
* `GEMTools binary bundle 1.6 for core2 <http://bla>`_

Downloads library bundle
************************

* `GEMTools binary bundle 1.6 for i3 <http://bla>`_
* `GEMTools binary bundle 1.6 for core2 <http://bla>`_

Souce Code
**********



Contents:

.. toctree::
   :maxdepth: 2

   command_line_tools
   rna_pipeline
   python_api
   c_api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

