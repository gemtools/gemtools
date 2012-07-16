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


/*Construction*/
int Map_init(Map *self, PyObject *args, PyObject *kwds);
PyObject* Map_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
// internal construction
Map* create_map(gt_map* map);
void Map_dealloc(Map* map);

// seq name
PyObject* Map_getseq_name(Map *self, void *closure);
int Map_setseq_name(Map *self, PyObject *value, void *closure);

// position
PyObject* Map_getposition(Map *self, void *closure);
int Map_setposition(Map *self, PyObject *value, void *closure);

// strand
int Map_setdirection(Map *self, PyObject *value, void *closure);
PyObject* Map_getdirection(Map *self, void *closure);

// map score
PyObject* Map_getscore(Map *self, void *closure);
int Map_setscore(Map *self, PyObject *value, void *closure);

// distance
PyObject* Map_getdistance(Map *self, void *closure);

// levenshtein
PyObject* Map_getlevenshtein(Map *self, void *closure);

// length of the read
PyObject* Map_getbase_length(Map *self, void *closure);
int Map_setbase_length(Map *self, PyObject *value, void *closure);

// length incl indels
PyObject* Map_getlength(Map *self, void *closure);

// global length 
PyObject* Map_getglobal_length(Map *self, void *closure);

// global distance
PyObject* Map_getglobal_distance(Map *self, void *closure);

// global score
PyObject* Map_getglobal_score(Map *self, void *closure);

// global score
PyObject* Map_getglobal_levenshtein(Map *self, void *closure);

// iterate the mismatches
PyObject* Map_getmismatches(Map *self, void *closure);
int Map_setmismatches(Map *self, PyObject *value, void *closure);
PyObject* Map_getnum_mismatches(Map *self, void *closure);



#endif
