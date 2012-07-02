#ifndef GEMPY_PY_TEMPLATE_ITERATOR
#define GEMPY_PY_TEMPLATE_ITERATOR

#include "Python.h"
#include "gem_tools.h"
#include "py_template.h"

typedef struct {
  PyObject_HEAD
  gt_template* template;
  gt_buffered_map_input* map_input;
  PyObject* tmpl;
} gempy_template_iterator;

gempy_template_iterator* create_template_stream_iterator(FILE* file);
gempy_template_iterator* create_template_file_iterator(char* filename, bool memorymap);

PyObject* gempy_template_iterator_iter(PyObject *self);
PyObject* gempy_template_iterator_iternext(PyObject *self);

static PyTypeObject gempy_template_iteratorType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "gempy._template_iterator",            /*tp_name*/
    sizeof(gempy_template_iterator),       /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    0,                         /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER,
      /* tp_flags: Py_TPFLAGS_HAVE_ITER tells python to
         use tp_iter and tp_iternext fields. */
    "Internal Template iterator",           /* tp_doc */
    0,  /* tp_traverse */
    0,  /* tp_clear */
    0,  /* tp_richcompare */
    0,  /* tp_weaklistoffset */
    gempy_template_iterator_iter,  /* tp_iter: __iter__() method */
    gempy_template_iterator_iternext  /* tp_iternext: next() method */
};



#endif
