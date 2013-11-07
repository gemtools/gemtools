from distribute_setup import use_setuptools
use_setuptools()

#
# make sur ewe use gcc on mac os x
# for openmp support
#
import platform
if platform.system() == "Darwin":
    from distutils import sysconfig
    sysconfig._config_vars['CC'] = "gcc"
    sysconfig._config_vars['LDSHARED'] = "gcc -Wl,-F. -bundle -undefined dynamic_lookup"

import os
import sys
from setuptools import setup, Command
#from setuptools.extension import Extension
from distutils.extension import Extension
from setuptools.command.install import install as _install
from setuptools.command.build_py import build_py as _build_py
#from setuptools.command.build_ext import build_ext as _build_ext
#from Cython.Distutils import build_ext as _build_ext
from distutils.command.build_ext import build_ext as _build_ext

import subprocess
import urllib
import tempfile
import shutil


__VERSION_MAJOR = "1"
__VERSION_MINOR = "7.1"
__VERSION__ = "%s.%s" % (__VERSION_MAJOR, __VERSION_MINOR)
__DOWNLOAD_FILE_TEMPLATE__ = "GEM-gemtools-%s-%s.tar.gz"
__DOWNLOAD_URL__ = "http://barnaserver.com/gemtools/"


def _is_i3_compliant():
    """Reads lines from /proc/cpuinfo and scans "flags" lines.
    Returns true if the flags are compatible with the GEM
    i3 bundle.
    """
    if not os.path.exists("/proc/cpuinfo"):
        return False
    stream = open("/proc/cpuinfo", 'r')
    i3_flags = set(["popcnt", "ssse3", "sse4_1", "sse4_2"])
    cpu_flags = set([])
    for line in iter(stream.readline, ''):
        line = line.rstrip()
        if line.startswith("flags"):
            for e in line.split(":")[1].strip().split(" "):
                cpu_flags.add(e)
    stream.close()
    return i3_flags.issubset(cpu_flags)


def download(type, target_dir=None):
    file_name = __DOWNLOAD_FILE_TEMPLATE__ % (__VERSION__, type)
    base_url = "%s/%s" % (__DOWNLOAD_URL__, file_name)

    if target_dir is None:
        parent_dir = os.path.split(os.path.abspath(__file__))[0]
        dirpath = "%s/%s" % (parent_dir, "downloads")
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)
        # support jenkins builds with latest artefacts
        target = "%s/%s" % (dirpath, file_name)
    else:
        target = "%s/%s" % (target_dir, file_name)

    parent_dir = os.path.split(os.path.abspath(target))[0]
    jenkins_src = "%s/current-%s.tar.gz" % (parent_dir, type)

    if os.path.exists(jenkins_src):
        print "copy jenkins build for", type
        shutil.copyfile(jenkins_src, target)
    else:
        if not os.path.exists(target):
            print "Downloading %s bundle from %s to %s" % (type, base_url, target)
            urllib.urlretrieve(base_url, target)


class fetch(Command):
    """Fetch binaries  package"""
    description = "Fetch binaries"
    user_options = []

    def run(self):
        print "Fetching binaries"
        download("i3")
        download("core2")

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass


class package(Command):
    """Package distribution"""
    description = "Package distribution"
    user_options = []

    def run(self):
        print "Fetching binaries"
        download("i3")
        download("core2")
        subprocess.Popen(["dist-utils/create_distribution.sh", __VERSION__, "i3"]).wait()
        subprocess.Popen(["dist-utils/create_distribution.sh", __VERSION__, "core2"]).wait()

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

class package_static(Command):
    """Package static distribution"""
    description = "Package static distribution"
    user_options = []

    def run(self):
        print "Fetching binaries"
        download("i3")
        download("core2")
        subprocess.Popen(["dist-utils/create_static_distribution.sh", __VERSION__, "i3"]).wait()
        subprocess.Popen(["dist-utils/create_static_distribution.sh", __VERSION__, "core2"]).wait()

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass


def _install_bundle(install_dir, base=None):
    """Download GEM bundle and move
    the bundled executables into the given target directory
    """
    if install_dir is None:
        print "Unable to determine installation directory !"
        exit(1)

    # pick the type
    type = "core2"
    if _is_i3_compliant():
        type = "i3"

    if not os.path.exists(install_dir):
        os.mkdir(install_dir)

    file_name = __DOWNLOAD_FILE_TEMPLATE__ % (__VERSION__, type)
    base_url = "%s/%s" % (__DOWNLOAD_URL__, file_name)

    if base is not None:
        dirpath = base
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)
    else:
        dirpath = tempfile.mkdtemp()

    target = "%s/%s" % (dirpath, file_name)
    keep = False

    # check jenksin currents
    jenkins_src = "%s/current-%s.tar.gz" % (dirpath, type)
    if os.path.exists(jenkins_src):
        target = jenkins_src
        keep = True
    else:
        if os.path.exists(file_name):
            target = file_name
            keep = True
        else:
            if not os.path.exists(target):
                download(type, dirpath)

    tar = subprocess.Popen("tar xzf %s --exclude \"._*\"" % (os.path.abspath(target)), shell=True, cwd=dirpath)
    if tar.wait() != 0:
        print "Error while extracting gem bundle"
        exit(1)
    if base is None and not keep:
        os.remove(target)

    bins = [x for x in os.listdir(dirpath)]
    for file in bins:
        if not file.endswith("gz"):
            print "Move bundle library: %s to %s" % (file, install_dir)
            result_file = "%s/%s" % (install_dir, file)
            if os.path.exists(result_file):
                os.remove(result_file)
            shutil.move("%s/%s" % (dirpath, file), install_dir)
            os.chmod(result_file, 0755)
    # copy gt.* binaries
    bins = [x for x in os.listdir("GEMTools/bin")]
    for file in bins:
        if not file.endswith("gz") and file != "gt.construct":
            print "Copy binary library: %s to %s" % (file, install_dir)
            result_file = "%s/%s" % (install_dir, file)
            if os.path.exists(result_file):
                os.remove(result_file)
            shutil.copy("%s/%s" % ("GEMTools/bin", file), install_dir)
            os.chmod(result_file, 0755)


    # remove temp directory
    if base is None:
        os.removedirs(dirpath)


## install bundle command
class install_bundle(Command):
    def run(self):
        parent_dir = os.path.split(os.path.abspath(__file__))[0]
        target_dir = "%s/%s" % (parent_dir, "python/gem/gembinaries")
        _install_bundle(target_dir, base=parent_dir + "/downloads")


# hack the setup tools installation
# to make sure bundled binaries are
# executable after install
class install(_install):

    def run(self):
        _install.run(self)
        if os.getenv("GEM_NO_BUNDLE", None) is None:
            # find target folder
            install_dir = None
            for file in self.get_outputs():
                if file.endswith("/gem/__init__.py"):
                    install_dir = "%s/gembinaries" % os.path.split(file)[0]
                    break

            _install_bundle(install_dir)


def compile_gemtools():
        process = subprocess.Popen(['make'], shell=True, cwd='GEMTools')
        if process.wait() != 0:
            print >> sys.stderr, """

Error while compiling GEMTools. That is very unfortunate.

A possible reason might be a missing dependency. Please take a look at the lines
before this one. You need the following programs and libraries installed to compile
the GEMTools library.

Programms needed:
    * make
    * gcc

Libraris needed:
    * python-dev (the python headers and include files)
    * libbz2-dev (for bz compression support)

On a Debian/Ubuntu system you should be able to get all needed dependencies with:

sudo apt-get install make gcc python-dev libbz2-dev

"""
            exit(1)

class build_ext(_build_ext):
    """Custom implementation of the extension builder to make sure that
    libgemtools is build and to trigger the cython build AFTER cython dependency
    was installed."""
    def run(self):
        compile_gemtools()
        try:
            """Fix the extensions, the .pyx file extensions are changed
            to .c somewhere along the line. This is a DIRTY hack and
            at some point we have to figure a better way"""
            for ext in self.extensions:
                new_sources = []
                for source in ext.sources:
                    output_dir = os.path.dirname(source)
                    (base, ex) = os.path.splitext(os.path.basename(source))
                    if not os.path.exists(source) and ex == ".c":
                        ex = ".pyx"

                    new_sources.append(os.path.join(output_dir, base + ex))
                ext.sources = new_sources

            # import cython, create a new instance of cythons build_ext
            # initialize it and run it
            from Cython.Distutils import build_ext as c_build_ext
            cb = c_build_ext(self.distribution)
            cb.finalize_options()
            cb.inplace = self.inplace
            cb.run()
        except ImportError:
            sys.stderr.write("\nERROR: Unable to load Cython builder. Please make sure Cython is installed on your system.\n\n")
            exit(1)

    def build_extensions(self):
        pass


class build_py(_build_py):
    def run(self):
        compile_gemtools()
        parent_dir = os.path.split(os.path.abspath(__file__))[0]
        target_dir = "%s/%s" % (parent_dir, "python/gem/gembinaries")
        _install_bundle(target_dir, base=parent_dir + "/downloads")
        _build_py.run(self)


_commands = {'install': install, 'build_ext': build_ext, 'fetch': fetch, 'package': package, 'package_static': package_static, 'build_py': build_py}

# extend nosetests command to
# ensure we have the bundle installed and
# locally
try:
    from nose.commands import nosetests as _nosetests

    class nosetests(_nosetests):
        def run(self):
            parent_dir = os.path.split(os.path.abspath(__file__))[0]
            target_dir = "%s/%s" % (parent_dir, "python/gem/gembinaries")
            _install_bundle(target_dir, base=parent_dir + "/downloads")
            _nosetests.run(self)
    _commands['nosetests'] = nosetests
except:
    pass

gemtools = Extension("gem.gemtools", sources=["python/src/gemtools_binding.c", "python/src/gemtools.pyx", "python/src/gemapi.pxd"],
                    include_dirs=['GEMTools/include', 'GEMTools/resources/include/'],
                    library_dirs=['GEMTools/lib'],
                    libraries=['z', 'bz2', 'gemtools'],
                    extra_compile_args=['-fopenmp'],
                    extra_link_args=["-fopenmp"]
)


setup(
        cmdclass=_commands,
        name='Gemtools',
        version=__VERSION__,
        description='Python support library for the GEM mapper and the gemtools library',
        author='Thasso Griebel, Santiago Marco Sola',
        author_email='thasso.griebel@gmail.com',
        url='https://github.com/gemtools/gemtools',
        license="GNU General Public License (GPL)",
        long_description='''This is the python binding and wrapper library around the GEM mapper.
The module allows you to run teh GEM mapper and simplifies building mapping
pipeline in python. In addition, we provide a fast C based parsing library that
is used to parse GEM results and extract mapping information.

For more information about the GEM see

http://algorithms.cnag.cat/wiki

The code for this project can be found on github:

https://github.com/gemtools/gemtools
''',
        package_dir={'': 'python'},
        packages=['gem'],
#        package_data={"": ["%s/%s" % ("gem/gembinaries",x) for x in os.listdir("python/gem/gembinaries")]},
        package_data={"": ["%s/%s" % ("python/gem/gembinaries", x) for x in ["gem-2-sam",
                                                                     "gem-indexer_bwt-dna",
                                                                     "gem-indexer_generate",
                                                                     "gem-2-gem",
                                                                     "gem-retriever",
                                                                     "gtf-2-junctions",
                                                                     "gem-indexer",
                                                                     "gem-indexer_fasta2meta+cont",
                                                                     "gem-info",
                                                                     "gem-mapper",
                                                                     "gem-rna-tools"
                                                                     ]]},
        ext_modules=[gemtools],
        test_suite='nose.collector',
        zip_safe=False,
        include_package_data=False,
        platforms=['lx64'],
        classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Developers',
          'License :: OSI Approved :: GNU General Public License (GPL)',
          'Operating System :: POSIX :: Linux',
          'Programming Language :: Python',
          'Programming Language :: C',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        setup_requires=["cython==0.18"],
        install_requires=[
                "argparse",
#                "numpy==1.7.0",
#                "matplotlib==1.2.0"
                ],
        entry_points={
            'console_scripts': [
                'gemtools = gem.commands:gemtools'
            ]
        },
)
