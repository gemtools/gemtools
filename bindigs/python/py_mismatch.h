#ifndef GEMPY_PY_MISMATCH
#define GEMPY_PY_MISMATCH

#include "Python.h"
#include "structmember.h"
#include "gem_tools.h"

/*
 * Python mismatch structure
 */
typedef struct {
    PyObject_HEAD
    gt_misms* misms;
} Mismatch;
/*
 * Create a fully initialized new mismatch instance from the
 * given gr_misms struct
 */
Mismatch* create_mismatch(gt_misms* map);
static int Mismatch_init(Mismatch *self, PyObject *args, PyObject *kwds);
static PyObject* Mismatch_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
static PyObject* Mismatch_getposition(Mismatch *self, void *closure);
static int Mismatch_setposition(Mismatch *self, PyObject *value, void *closure);
static PyObject* Mismatch_getsize(Mismatch *self, void *closure);
static int Mismatch_setsize(Mismatch *self, PyObject *value, void *closure);
static PyObject* Mismatch_gettype(Mismatch *self, void *closure);
static int Mismatch_settype(Mismatch *self, PyObject *value, void *closure);
static PyObject* Mismatch_getbase(Mismatch *self, void *closure);
static int Mismatch_setbase(Mismatch *self, PyObject *value, void *closure);
static PyGetSetDef Mismatch_getseters[] = {
    {"type", (getter) Mismatch_gettype, (setter) Mismatch_settype, "Set type", NULL},
    {"position", (getter) Mismatch_getposition, (setter) Mismatch_setposition, "Set position", NULL},
    {"base", (getter) Mismatch_getbase, (setter) Mismatch_setbase, "Set base", NULL},
    {"size", (getter) Mismatch_getsize, (setter) Mismatch_setsize, "Set size", NULL},
    {NULL}  /* Sentinel */
};
static PyTypeObject MismatchType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "gempy.Mismatch",          /*tp_name*/
    sizeof(Mismatch),          /*tp_basicsize*/
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
    "Mismatch objects",       /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    0,                     /* tp_methods */
    0,                     /* tp_members */
    Mismatch_getseters,    /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Mismatch_init,      /* tp_init */
    0,                         /* tp_alloc */
    Mismatch_new,                 /* tp_new */
};
#endif
