#ifndef GEMPY_PY_ITERATOR
#define GEMPY_PY_ITERATOR

#include "Python.h"

typedef PyObject* (*gempy_converter_function)(void*,...);
typedef void* (*gempy_getter_function)(void*, uint64_t);

typedef struct {
    PyObject_HEAD
    uint64_t start;
    uint64_t length;
    uint64_t pos;
    gempy_getter_function getter;
    void* arg;
    gempy_converter_function converter;
    int parent;
} gempy_iterator;

PyObject* gempy_iterator_iter(PyObject *self);
PyObject* gempy_iterator_iternext(PyObject *self);
static PyObject* create_gempy_iterator(uint64_t start, uint64_t length, void* getter, void* arg, void* converter, int parent);

static PyTypeObject gempy_iteratorType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "gempy._iterator",            /*tp_name*/
    sizeof(gempy_iterator),       /*tp_basicsize*/
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
    "Internal iterator object.",           /* tp_doc */
    0,  /* tp_traverse */
    0,  /* tp_clear */
    0,  /* tp_richcompare */
    0,  /* tp_weaklistoffset */
    gempy_iterator_iter,  /* tp_iter: __iter__() method */
    gempy_iterator_iternext  /* tp_iternext: next() method */
};


#endif
