#include "py_iterator.h"

PyObject* gempy_iterator_iter(PyObject *self){
  Py_INCREF(self);
  return self;
}

PyObject* gempy_iterator_iternext(PyObject *self){
    gempy_iterator *p = (gempy_iterator *) self;
    if(p->pos >= p->length){
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
    }
    gempy_getter_function getter = p->getter;
    gempy_converter_function converter = p->converter;
    uint64_t s = p->pos + p->start;
    p->pos++;
    void* result = getter(p->arg, s);
    if(!result){
        Py_RETURN_NONE;
    }else{
        if(p->parent){
            return converter(result, p->arg);
        }else{
            return converter(result);
        }
    }
}


