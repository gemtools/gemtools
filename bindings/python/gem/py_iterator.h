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

void gempy_iterator_dealloc(PyObject* self);
PyObject* gempy_iterator_iter(PyObject *self);
PyObject* gempy_iterator_iternext(PyObject *self);
Py_ssize_t gempy_iterator_len(PyObject* self);
PyObject* create_gempy_iterator(uint64_t start, uint64_t length, void* getter, void* arg, void* converter, int parent);




#endif
