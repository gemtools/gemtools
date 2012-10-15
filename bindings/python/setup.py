import os
from setuptools import setup
from distutils.core import Extension

gemtools_dir = "../../GEMTools/"
objs = []
for f in os.listdir(gemtools_dir):
    if f.endswith(".o"):
        objs.append(gemtools_dir+f)



gemtools = Extension('gem.gemtools',
                    define_macros=[('MAJOR_VERSION', '1'),
                                   ('MINOR_VERSION', '3')],
                    include_dirs=['../../GEMTools'],
                    extra_objects=objs,
                    sources=['src/py_iterator.c', 'src/py_template_iterator.c',
                               'src/py_mismatch.c', 'src/py_map.c', 'src/py_alignment.c',
                               'src/py_template.c', 'src/gemtoolsmodule.c', 'src/py_mappings_iterator.c'])

setup(
        name='Gem',
        version='1.3',
        description='Python support library for the GEM mapper and gemtools',
        author='Thasso Griebel',
        author_email='thasso.griebel@gmail.com',
        url='http://barnaserver.com/gemtools',
        long_description='''
        This is the python binding to the gemtools library.
        ''',
        packages=['gem'],
        data_files=[("gem/gembinaries/", ["%s/%s" % ("gem/gembinaries",x) for x in os.listdir("gem/gembinaries")])],
        ext_modules=[gemtools],
        setup_requires=['nose'],
        test_suite = 'nose.collector',
)
