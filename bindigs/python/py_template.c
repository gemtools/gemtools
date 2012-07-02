#include "py_template.h"

static int Template_init(Template *self, PyObject *args, PyObject *kwds){
    self->template = gt_template_new();
    return 0;
}

static PyObject* Template_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
    Template *self = (Template *)type->tp_alloc(type, 0);
    Template_init(self, args, kwds);
    return (PyObject *)self;
}

static PyObject* Template_gettag(Template *self, void *closure){
    return PyString_FromString(gt_template_get_tag(self->template));
}

static int Template_settag(Template *self, PyObject *value, void *closure){
    gt_template_set_tag(self->template, PyString_AsString(value));
    return 0;
}

static PyObject* Template_getmax_complete_strata(Template *self, void *closure){
    return PyLong_FromUnsignedLongLong(gt_template_get_mcs(self->template));
}

static int Template_setmax_complete_strata(Template *self, PyObject *value, void *closure){
    gt_template_set_mcs(self->template, PyLong_AsUnsignedLongLong(value));
    return 0;
}

static PyObject* Template_getblocks(Template *self, void *closure){
    return create_gempy_iterator(0, gt_template_get_num_blocks(self->template), gt_template_get_block, self->template, create_alignment, 1);
}

static int Template_setblocks(Template *self, PyObject *value, void *closure){
    PyErr_SetString(PyExc_TypeError, "Setting blocks is currently not supported");
    return -1;
}

static PyObject* Template_getcounters(Template *self, void *closure){
    return create_gempy_iterator(1, gt_template_get_num_counters(self->template), gt_template_get_counter, self->template, PyLong_FromUnsignedLongLong, 0);
}

static int Template_setcounters(Template *self, PyObject *value, void *closure){
    PyErr_SetString(PyExc_TypeError, "Setting blocks is currently not supported");
    return -1;
}

Template* create_template(gt_template* template){
    Template* tmpl = PyObject_New(Template, &TemplateType);
    tmpl->template = template;
    return tmpl;
}
