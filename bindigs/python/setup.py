from distutils.core import setup, Extension

module1 = Extension('gempy',
                    define_macros=[('MAJOR_VERSION', '1'),
                                   ('MINOR_VERSION', '0')],
                    include_dirs=['../../GEMTools'],
                    extra_objects=['../../GEMTools/gem_tools.o', '../../GEMTools/gt_error.o',
                                   '../../GEMTools/gt_output_map.o', '../../GEMTools/gt_alignment.o',
                                   '../../GEMTools/gt_input_file.o', '../../GEMTools/gt_template.o',
                                   '../../GEMTools/gt_buffered_map_input.o', '../../GEMTools/gt_map.o',
                                   '../../GEMTools/gt_vector.o', '../../GEMTools/gt_commons.o',
                                   '../../GEMTools/gt_misms.o'],
                    sources=['py_iterator.c', 'py_template_iterator.c',
                               'py_mismatch.c', 'py_map.c', 'py_alignment.c',
                               'py_template.c', 'gempymodule.c'])

setup(name='GemPy',
       version='1.0',
       description='Pythobn binding to gemtools',
       author='Thasso Griebel',
       author_email='thasso.griebel@gmail.com',
       url='http://barnaserver.com/gemtools',
       long_description='''
This is the python binding to the gemtools library.
''',
       ext_modules=[module1])
