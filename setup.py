import os
from setuptools import setup, Command
from distutils.core import Extension
from setuptools.command.install import install as _install
from setuptools.command.build_ext import build_ext as _build_ext
import subprocess
import urllib
import tempfile
import shutil


__VERSION_MAJOR="1"
__VERSION_MINOR="6"
__VERSION__="%s.%s" %(__VERSION_MAJOR, __VERSION_MINOR)
__DOWNLOAD_FILE_TEMPLATE__="GEM-gemtools-%s-%s.tgz"
__DOWNLOAD_URL__="https://github.com/downloads/gemtools/gemtools"


def _is_i3_compliant():
    """Reads lines from /proc/cpuinfo and scans "flags" lines.
    Returns true if the flags are compatible with the GEM
    i3 bundle.
    """
    if not os.path.exists("/proc/cpuinfo"): return False
    stream = open("/proc/cpuinfo", 'r')
    i3_flags = set(["popcnt", "ssse3", "sse4_1", "sse4_2"])
    cpu_flags = set([])
    for line in iter(stream.readline,''):
        line = line.rstrip()
        if line.startswith("flags"):
            for e in line.split(":")[1].strip().split(" "):
                cpu_flags.add(e)
    stream.close()
    return i3_flags.issubset(cpu_flags)

def _install_bundle(install_dir):
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
    base_url = "%s/%s" % (__DOWNLOAD_URL__,file_name)
    dirpath = tempfile.mkdtemp()
    target = "%s/%s" %(dirpath, file_name)
    print "Downloading %s bundle from %s to %s" % (type, base_url, target)
    urllib.urlretrieve (base_url, target)
    tar = subprocess.Popen("tar xzvf %s" % (target), shell=True, cwd=dirpath)
    if tar.wait() != 0:
        print "Error while extracting gem bundle"
        exit(1)
    os.remove(target)

    bins=[x for x in os.listdir(dirpath)]

    for file in bins:
        print "Move bundle library: %s to %s" % (file, install_dir)
        result_file = "%s/%s" % (install_dir, file)
        if os.path.exists(result_file):
            os.remove(result_file)
        shutil.move("%s/%s" % (dirpath, file), install_dir)
        os.chmod(result_file, 0755)

    # remove temp directory
    os.removedirs(dirpath)



# hack the setup tools installation
# to make sure bundled binaries are
# executable after install
class install(_install):
    def run(self):
        _install.run(self)

        # find target folder
        install_dir = None
        for file in self.get_outputs():
            if file.endswith("gem/__init__.py"):
                install_dir = "%s/gembinaries" % os.path.split(file)[0]
                break

        _install_bundle(install_dir)

class build_ext(_build_ext):
    def run(self):
        process = subprocess.Popen(['make'], shell=True, cwd='GEMTools')
        if process.wait() != 0:
            raise ValueError("Error while compiling GEMTools")
        _build_ext.run(self)


gemtools = Extension('gem.gemtools',
                    define_macros=[('MAJOR_VERSION', __VERSION_MAJOR),
                                   ('MINOR_VERSION', __VERSION_MINOR)],
                    include_dirs=['GEMTools/include', 'GEMTools/resources/include/'],
                    library_dirs=['GEMTools/lib'],
                    libraries=['gemtools'],
                    sources=['python/src/py_iterator.c', 'python/src/py_template_iterator.c',
                               'python/src/py_mismatch.c', 'python/src/py_map.c', 'python/src/py_alignment.c',
                               'python/src/py_template.c', 'python/src/gemtoolsmodule.c', 'python/src/py_mappings_iterator.c'])

setup(
        cmdclass={'install': install, 'build_ext': build_ext},
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
        package_dir={'':'python'},
        packages=['gem'],
#        package_data={"": ["%s/%s" % ("gem/gembinaries",x) for x in os.listdir("python/gem/gembinaries")]},
        package_data={"": ["%s/%s" % ("gem/gembinaries",x) for x in ["gem-2-sam",
                                                                     "gem-indexer_bwt-dna",
                                                                     "gem-indexer_generate",
                                                                     "gem-map-2-map",
                                                                     "gem-retriever",
                                                                     "gtf-2-junctions",
                                                                     "gem-indexer",
                                                                     "gem-indexer_fasta2meta+cont",
                                                                     "gem-info",
                                                                     "gem-mapper",
                                                                     "gem-rna-mapper",
                                                                     "splits-2-junctions"]]},
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
)
