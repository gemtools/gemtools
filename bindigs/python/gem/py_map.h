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


#endif
