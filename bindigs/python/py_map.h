#ifndef GEMPY_PY_MAP
#define GEMPY_PY_MAP

#include "Python.h"
#include "structmember.h"
#include "gem_tools.h"
#include "py_mismatch.h"

typedef struct {
    PyObject_HEAD
    gt_map* map;
} Map;

Map* create_map(gt_map* map);

int Map_init(Map *self, PyObject *args, PyObject *kwds);
PyObject* Map_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* Map_getseq_name(Map *self, void *closure);
int Map_setseq_name(Map *self, PyObject *value, void *closure);
PyObject* Map_getposition(Map *self, void *closure);
int Map_setposition(Map *self, PyObject *value, void *closure);
PyObject* Map_getscore(Map *self, void *closure);
int Map_setscore(Map *self, PyObject *value, void *closure);
PyObject* Map_getbase_lengt(Map *self, void *closure);
int Map_setbase_length(Map *self, PyObject *value, void *closure);
int Map_setdirection(Map *self, PyObject *value, void *closure);
PyObject* Map_getmismatches(Map *self, void *closure);
int Map_setmismatches(Map *self, PyObject *value, void *closure);
PyObject* Map_getdirection(Map *self, void *closure);
int Map_setdirection(Map *self, PyObject *value, void *closure);

static PyGetSetDef Map_getseters[] = {
    {"seq_name", (getter) Map_getseq_name, (setter) Map_setseq_name, "Genomic sequence name", NULL},
    {"position", (getter) Map_getposition, (setter) Map_setposition, "Genomic Position", NULL},
    {"base_length", (getter) Map_getbase_lengt, (setter) Map_setbase_length, "Base length", NULL},
    {"direction", (getter) Map_getdirection, (setter) Map_setdirection, "Strand", NULL},
    {"score", (getter) Map_getscore, (setter) Map_setscore, "Score", NULL},
    {"mismatches", (getter) Map_getmismatches, (setter) Map_setmismatches, "Mismatches", NULL},
    {NULL}  /* Sentinel */
};

static PyTypeObject MapType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "gempy.Map",             /*tp_name*/
    sizeof(Map),             /*tp_basicsize*/
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
    "Map objects",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    0,             /* tp_methods */
    0,             /* tp_members */
    Map_getseters,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Map_init,      /* tp_init */
    0,                         /* tp_alloc */
    Map_new,                 /* tp_new */
};
#endif
