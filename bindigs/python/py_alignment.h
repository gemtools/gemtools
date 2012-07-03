#ifndef GEMPY_PY_ALIGNMENT
#define GEMPY_PY_ALIGNMENT
/*
 * Alignment class
 */
#include "Python.h"
#include "structmember.h"
#include "gem_tools.h"

typedef struct {
    PyObject_HEAD
    gt_alignment* alignment;
    gt_template* template;
} Alignment;

Alignment* create_alignment(gt_alignment* alignment, gt_template* parent);

int Alignment_init(Alignment *self, PyObject *args, PyObject *kwds);
PyObject* Alignment_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* Alignment_gettag(Alignment *self, void *closure);
int Alignment_settag(Alignment *self, PyObject *value, void *closure);
PyObject* Alignment_getread(Alignment *self, void *closure);
int Alignment_setread(Alignment *self, PyObject *value, void *closure);
PyObject* Alignment_getqualities(Alignment *self, void *closure);
int Alignment_setqualities(Alignment *self, PyObject *value, void *closure);
PyObject* Alignment_getmax_complete_strata(Alignment *self, void *closure);
int Alignment_setmax_complete_strata(Alignment *self, PyObject *value, void *closure);
PyObject* Alignment_getcounters(Alignment *self, void *closure);
int Alignment_setcounters(Alignment *self, PyObject *value, void *closure);



static PyGetSetDef Alignment_getseters[] = {
    {"tag", (getter) Alignment_gettag, (setter) Alignment_settag, "Alignment Tag", NULL},
    {"read", (getter) Alignment_getread, (setter) Alignment_setread, "Alignment Read", NULL},
    {"qualities", (getter) Alignment_getqualities, (setter) Alignment_setqualities, "Alignment Qualities", NULL},
    {"max_complete_strata", (getter) Alignment_getmax_complete_strata, 
        (setter) Alignment_setmax_complete_strata, "Max Complete Strata", NULL},
    {"counters", (getter) Alignment_getcounters, 
        (setter) Alignment_setcounters, "Counters", NULL},
    {NULL}  /* Sentinel */
};

static PyTypeObject AlignmentType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "gempy.Alignment",             /*tp_name*/
    sizeof(Alignment),             /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    0, /*tp_dealloc*/
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
    "Alignment objects",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    0,             /* tp_methods */
    0,             /* tp_members */
    Alignment_getseters,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Alignment_init,      /* tp_init */
    0,                         /* tp_alloc */
    Alignment_new,                 /* tp_new */
};
#endif
