#include "py_mismatch.h"

int Mismatch_init(Mismatch *self, PyObject *args, PyObject *kwds){
    gt_misms* ms = gt_misms_new();
    self->misms = ms;
    return 0;
}

PyObject* Mismatch_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
    Mismatch *self;
    self = (Mismatch *)type->tp_alloc(type, 0);
    Mismatch_init(self, args, kwds);
    return (PyObject *)self;
}

void Mismatch_dealloc(PyObject* self){
    self->ob_type->tp_free(self);
}

PyObject* Mismatch_getposition(Mismatch *self, void *closure){
    PyObject* ret = PyLong_FromUnsignedLongLong(gt_misms_get_position(self->misms));
    return ret;
}

int Mismatch_setposition(Mismatch *self, PyObject *value, void *closure){
    //gt_misms_set_position(self->misms, PyLong_AsUnsignedLongLong(value));
    return 0;
}

PyObject* Mismatch_getsize(Mismatch *self, void *closure){
    if(gt_misms_get_type(self->misms) != MISMS){
      PyObject* ret = PyLong_FromUnsignedLongLong(gt_misms_get_size(self->misms));
      return ret;
    }
    return PyLong_FromLongLong(-1);
}

int Mismatch_setsize(Mismatch *self, PyObject *value, void *closure){
    //gt_misms_set_size(self->misms, PyLong_AsUnsignedLongLong(value));
    return 0;
}

PyObject* Mismatch_gettype(Mismatch *self, void *closure){
    PyObject* ret = PyInt_FromLong(gt_misms_get_type(self->misms));
    return ret;
}

int Mismatch_settype(Mismatch *self, PyObject *value, void *closure){
    //gt_misms_set_type(self->misms, PyInt_AsLong(value));
    return 0;
}

PyObject* Mismatch_getbase(Mismatch *self, void *closure){
    PyObject* ret = Py_BuildValue("c", gt_misms_get_base(self->misms));
    return ret;
}

int Mismatch_setbase(Mismatch *self, PyObject *value, void *closure){
    //gt_misms_set_base(self->misms, PyInt_AsInt(value));
    return 0;
}
