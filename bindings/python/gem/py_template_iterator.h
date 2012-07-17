#ifndef GEMPY_PY_TEMPLATE_ITERATOR
#define GEMPY_PY_TEMPLATE_ITERATOR

#include "Python.h"
#include "gem_tools.h"

typedef struct {
  PyObject_HEAD
  gt_template* template;
  gt_buffered_map_input* map_input;
  gt_input_file* input_file;
  PyObject* tmpl;
} gempy_template_iterator;

gempy_template_iterator* create_template_stream_iterator(FILE* file);
gempy_template_iterator* create_template_file_iterator(char* filename, bool memorymap);

PyObject* gempy_template_iterator_iter(PyObject *self);
PyObject* gempy_template_iterator_iternext(PyObject *self);
void gempy_template_iterator_dealloc(gempy_template_iterator* self);

#endif
