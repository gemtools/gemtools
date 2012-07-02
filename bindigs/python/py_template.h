#ifndef GEMPY_PY_TEMPLATE
#define GEMPY_PY_TEMPLATE

#include "Python.h"
#include "structmember.h"
#include "gem_tools.h"

typedef struct {
    PyObject_HEAD
    gt_template* template;
} Template;


static int Template_init(Template *self, PyObject *args, PyObject *kwds);
static PyObject* Template_new(PyTypeObject *type, PyObject *args, PyObject *kwds);

static PyObject* Template_gettag(Template *self, void *closure);
static int Template_settag(Template *self, PyObject *value, void *closure);
static PyObject* Template_getmax_complete_strata(Template *self, void *closure);
static int Template_setmax_complete_strata(Template *self, PyObject *value, void *closure);
static PyObject* Template_getblocks(Template *self, void *closure);
static int Template_setblocks(Template *self, PyObject *value, void *closure);
static PyObject* Template_getcounters(Template *self, void *closure);
static int Template_setcounters(Template *self, PyObject *value, void *closure);
Template* create_template(gt_template* template);

static PyGetSetDef Template_getseters[] = {
    {"tag", (getter) Template_gettag, (setter) Template_settag, "Template Tag", NULL},
    {"max_complete_strata", (getter) Template_getmax_complete_strata, 
        (setter) Template_setmax_complete_strata, "Max Complete Strata", NULL},
    {"blocks", (getter) Template_getblocks, 
        (setter) Template_setblocks, "Alignment blocks", NULL},
    {"counters", (getter) Template_getcounters, 
        (setter) Template_setcounters, "Counters", NULL},
    {NULL}  /* Sentinel */
};

static PyTypeObject TemplateType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "gempy.Template",             /*tp_name*/
    sizeof(Template),             /*tp_basicsize*/
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
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "Template objects",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    0,             /* tp_methods */
    0,             /* tp_members */
    Template_getseters,        /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Template_init,      /* tp_init */
    0,                         /* tp_alloc */
    Template_new,                 /* tp_new */
};
#endif
