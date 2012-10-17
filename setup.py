import os
from setuptools import setup
from distutils.core import Extension
from setuptools.command.install import install as _install
from setuptools.command.build_ext import build_ext as _build_ext
import sys
import subprocess

# hack the setup tools installation
# to make sure bundled binaries are 
# executable after install
class install(_install):
    def run(self):
        _install.run(self)
        bins=[x for x in os.listdir("python/gem/gembinaries")]
        for file in self.get_outputs():
            if os.path.basename(file) in bins:
                print "Making %s executable" % file
                os.chmod(file, 0755)


class build_ext(_build_ext):
    def run(self):
        process = subprocess.Popen(['make'], shell=True, cwd='GEMTools')
        if process.wait() != 0:
            raise ValueError("Error while compiling GEMTools")
        _build_ext.run(self)


gemtools = Extension('gem.gemtools',
                    define_macros=[('MAJOR_VERSION', '1'),
                                   ('MINOR_VERSION', '3')],
                    include_dirs=['GEMTools/include', 'GEMTools/resources/include/'],
                    library_dirs=['GEMTools/lib'],
                    libraries=['gemtools'],
                    sources=['python/src/py_iterator.c', 'python/src/py_template_iterator.c',
                               'python/src/py_mismatch.c', 'python/src/py_map.c', 'python/src/py_alignment.c',
                               'python/src/py_template.c', 'python/src/gemtoolsmodule.c', 'python/src/py_mappings_iterator.c'])

setup(
        cmdclass={'install': install, 'build_ext': build_ext},
        name='Gemtools',
        version='1.3',
        description='Python support library for the GEM mapper and the gemtools library',
        author='Thasso Griebel',
        author_email='thasso.griebel@gmail.com',
        url='http://algorithms.cnag.cat/',
        license="GNU Library or Lesser General Public License (LGPL)",
        long_description='''This is the python binding and wrapper library around the GEM mapper.
The module allows you to run teh GEM mapper and simplifies building mapping
pipeline in python. In addition, we provide a fast C based parsing library that
is used to parse GEM results and extract mapping information.

For more information about the GEM see

http://algorithms.cnag.cat/

The code for this project can be found on github:


''',
        package_dir={'':'python'},
        packages=['gem'],
        package_data={"": ["%s/%s" % ("gem/gembinaries",x) for x in os.listdir("python/gem/gembinaries")]},
        ext_modules=[gemtools],
        test_suite='nose.collector',
        zip_safe=False,
        include_package_data=True,
        platforms=['lx64'],
        classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Developers',
          'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)',
          'Operating System :: POSIX :: Linux',
          'Programming Language :: Python',
          'Programming Language :: C',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
)
