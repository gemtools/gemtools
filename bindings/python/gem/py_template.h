#ifndef GEMPY_PY_TEMPLATE
#define GEMPY_PY_TEMPLATE

#include "Python.h"
#include "structmember.h"
#include "gem_tools.h"

typedef struct {
    PyObject_HEAD
    gt_template* template;
} Template;

Template* create_template(gt_template* template);

int Template_init(Template *self, PyObject *args, PyObject *kwds);
PyObject* Template_new(PyTypeObject *type, PyObject *args, PyObject *kwds);

PyObject* Template_gettag(Template *self, void *closure);
int Template_settag(Template *self, PyObject *value, void *closure);
PyObject* Template_getmax_complete_strata(Template *self, void *closure);
int Template_setmax_complete_strata(Template *self, PyObject *value, void *closure);
PyObject* Template_getblocks(Template *self, void *closure);
int Template_setblocks(Template *self, PyObject *value, void *closure);
PyObject* Template_getcounters(Template *self, void *closure);
int Template_setcounters(Template *self, PyObject *value, void *closure);

//PyObject* Tempalte_iterate_mappings(PyObject* self, PyObject* args);


#endif
