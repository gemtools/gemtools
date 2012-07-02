from distutils.core import setup, Extension

module1 = Extension('gempy',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0')],
                    include_dirs = ['../../GEMTools'],
                    library_dirs = ['../../GEMTools'],
                    libraries = ['gemtools'],
                    sources = ['py_iterator.c', 'py_template_iterator.c',
                               'py_mismatch.c', 'py_map.c', 'py_alignment.c',
                               'py_template.c', 'gempymodule.c'])

setup (name = 'GemPy',
       version = '1.0',
       description = 'Pythobn binding to gemtools',
       author = 'Thasso Griebel',
       author_email = 'thasso.griebel@gmail.com',
       url = 'http://barnaserver.com/gemtools',
       long_description = '''
This is the python binding to the gemtools library.
''',
       ext_modules = [module1])
