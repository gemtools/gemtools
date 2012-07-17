#include "py_template.h"
#include "py_iterator.h"
#include "py_alignment.h"
#include "gem_tools.h"

int Template_init(Template *self, PyObject *args, PyObject *kwds){
    self->template = gt_template_new();
    return 0;
}

void Template_dealloc(PyObject* self){
    self->ob_type->tp_free(self);
}


PyObject* Template_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
    Template *self = (Template *)type->tp_alloc(type, 0);
    Template_init(self, args, kwds);
    return (PyObject *)self;
}

PyObject* Template_gettag(Template *self, void *closure){
    PyObject* ret = PyString_FromString(gt_template_get_tag(self->template));
    //Py_DECREF(ret);
    return ret;
}

int Template_settag(Template *self, PyObject *value, void *closure){
    gt_template_set_tag(self->template, PyString_AsString(value));
    return 0;
}

PyObject* Template_getmax_complete_strata(Template *self, void *closure){
    PyObject* ret = PyLong_FromUnsignedLongLong(gt_template_get_mcs(self->template));
    //Py_DECREF(ret);
    return ret;
}

int Template_setmax_complete_strata(Template *self, PyObject *value, void *closure){
    gt_template_set_mcs(self->template, PyLong_AsUnsignedLongLong(value));
    return 0;
}

PyObject* Template_getblocks(PyObject *self, PyObject *closure){
    Template* tmpl = (Template*) self;
    PyObject* ret = create_gempy_iterator(0, gt_template_get_num_blocks(tmpl->template), gt_template_get_block, tmpl->template, create_alignment, 1);
    //Py_DECREF(ret);
    return ret;
}

PyObject* Template_getcounters(Template *self, void *closure){
    PyObject* ret = create_gempy_iterator(1, gt_template_get_num_counters(self->template), gt_template_get_counter, self->template, PyLong_FromUnsignedLongLong, 0);
    //Py_DECREF(ret);
    return ret;
}

int Template_setcounters(Template *self, PyObject *value, void *closure){
    PyErr_SetString(PyExc_TypeError, "Setting blocks is currently not supported");
    return -1;
}


