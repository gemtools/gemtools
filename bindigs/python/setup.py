from ez_setup import use_setuptools

use_setuptools()
from setuptools import setup
from distutils.core import Extension


gemtools = Extension('gem.gemtools',
                    define_macros=[('MAJOR_VERSION', '1'),
                                   ('MINOR_VERSION', '0')],
                    include_dirs=['../../GEMTools'],
                    extra_objects=['../../GEMTools/gem_tools.o', '../../GEMTools/gt_error.o',
                                   '../../GEMTools/gt_output_map.o', '../../GEMTools/gt_alignment.o',
                                   '../../GEMTools/gt_input_file.o', '../../GEMTools/gt_template.o',
                                   '../../GEMTools/gt_buffered_map_input.o', '../../GEMTools/gt_map.o',
                                   '../../GEMTools/gt_vector.o', '../../GEMTools/gt_commons.o',
                                   '../../GEMTools/gt_misms.o'],
                    sources=['gem/py_iterator.c', 'gem/py_template_iterator.c',
                               'gem/py_mismatch.c', 'gem/py_map.c', 'gem/py_alignment.c',
                               'gem/py_template.c', 'gem/gemtoolsmodule.c'])


setup(
        name='Gem',
        version='1.0',
        description='Python support library for the GEM mapper and gemtools',
        author='Thasso Griebel',
        author_email='thasso.griebel@gmail.com',
        url='http://barnaserver.com/gemtools',
        long_description='''
        This is the python binding to the gemtools library.
        ''',
        packages=['gem'],
        ext_modules=[gemtools],
        setup_requires=['nose>=1.0'],
)
