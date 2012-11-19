#include "py_iterator.h"

//PyObject* gempy_iterator_iter(PyObject *self){
//  Py_INCREF(self);
//  return self;
//}
//
//void gempy_iterator_dealloc(PyObject* self){
//    self->ob_type->tp_free(self);
//}
//
//PyObject* gempy_iterator_iternext(PyObject *self){
//    gempy_iterator *p = (gempy_iterator *) self;
//    if(p->pos >= p->length){
//        PyErr_SetNone(PyExc_StopIteration);
//        return NULL;
//    }
//    gempy_getter_function getter = p->getter;
//    gempy_converter_function converter = p->converter;
//    uint64_t s = p->pos + p->start;
//    p->pos++;
//    void* result = getter(p->arg, s);
//    PyObject* ret = NULL;
//    if(p->parent){
//        ret = converter(result, p->arg, (p->length > 1) ? s+1 : 0);
//    }else{
//        ret = converter(result);
//    }
//    return ret;
//}
//Py_ssize_t gempy_iterator_len(PyObject *self){
//    gempy_iterator* i = (gempy_iterator*) self;
//    return i->length;
//}
//
//
