#ifndef GEMPY_PY_TEMPLATE
#define GEMPY_PY_TEMPLATE

#include "Python.h"
#include "structmember.h"
#include "gem_tools.h"

typedef struct {
    PyObject_HEAD
    gt_template* template;
} Template;

// allocation and creation
Template* create_template(gt_template* template);
int Template_init(Template *self, PyObject *args, PyObject *kwds);
void Template_dealloc(PyObject* self);
PyObject* Template_new(PyTypeObject *type, PyObject *args, PyObject *kwds);

/*
Fill the template from the input line
*/
PyObject* Template_fill(Template *self, PyObject *args);

//// set/get tag
//PyObject* Template_gettag(Template *self, void *closure);
//int Template_settag(Template *self, PyObject *value, void *closure);
//

//// get mcs
PyObject* Template_get_max_complete_strata(Template *self, void *closure);


//
//// num blocks
//PyObject* Template_get_num_blocks(Template *self, void *closure);
//// num counters
//PyObject* Template_get_num_counters(Template *self, void *closure);
//
//// get counters
//PyObject* Template_get_counters(PyObject *self, PyObject *closure);
//// get blocks
//PyObject* Template_get_blocks(PyObject* self, PyObject *closure);


#endif
