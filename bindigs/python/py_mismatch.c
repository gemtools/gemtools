#include "py_mismatch.h"

static int Mismatch_init(Mismatch *self, PyObject *args, PyObject *kwds){
    gt_misms* ms = malloc(sizeof(gt_misms));
    gt_cond_fatal_error(!ms,MEM_HANDLER);
    ms->misms_type = INS;
    ms->position = 0;
    ms->size = 0;
    self->misms = ms;
    return 0;
}

static PyObject* Mismatch_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
    Mismatch *self;
    self = (Mismatch *)type->tp_alloc(type, 0);
    Mismatch_init(self, args, kwds);
    return (PyObject *)self;
}

static PyObject* Mismatch_getposition(Mismatch *self, void *closure){
    return PyLong_FromUnsignedLongLong(gt_misms_get_position(self->misms));
}

static int Mismatch_setposition(Mismatch *self, PyObject *value, void *closure){
    gt_misms_set_position(self->misms, PyLong_AsUnsignedLongLong(value));
    return 0;
}

static PyObject* Mismatch_getsize(Mismatch *self, void *closure){
    return PyLong_FromUnsignedLongLong(gt_misms_get_size(self->misms));
}

static int Mismatch_setsize(Mismatch *self, PyObject *value, void *closure){
    gt_misms_set_size(self->misms, PyLong_AsUnsignedLongLong(value));
    return 0;
}

static PyObject* Mismatch_gettype(Mismatch *self, void *closure){
    return PyInt_FromLong(gt_misms_get_type(self->misms));
}

static int Mismatch_settype(Mismatch *self, PyObject *value, void *closure){
    gt_misms_set_type(self->misms, PyInt_AsLong(value));
    return 0;
}

static PyObject* Mismatch_getbase(Mismatch *self, void *closure){
    return Py_BuildValue("c", gt_misms_get_base(self->misms));
}

static int Mismatch_setbase(Mismatch *self, PyObject *value, void *closure){
    gt_misms_set_base(self->misms, PyInt_AsInt(value));
    return 0;
}


Mismatch* create_mismatch(gt_misms* map){
    Mismatch* a = PyObject_New(Mismatch, &MismatchType);
    a->misms = map;
    return a;
}
