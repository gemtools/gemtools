/*
 * Template class
 */
#include "Python.h"
#include "structmember.h"


typedef struct {
    PyObject_HEAD
    PyObject *tag; /* The read ID */
    PyObject *alignments;  /* List of alignments*/
    PyObject *counters; /*List of counters*/
    PyObject *multimaps; /*List of tuples for multimaps*/
} Template;

static void
Template_dealloc(Template* self)
{
    register uint64_t i = 0;
    register uint64_t size = 0;

    Py_XDECREF(self->tag);
    if(self->alignments){
        size = PyList_GET_SIZE(self->alignments);
        for(i=0; i<size;i++){
            Py_XDECREF(PyList_GET_ITEM(self->alignments, i));
        }
    }
    Py_XDECREF(self->alignments);
    /*
    if(self->counters){
        size = PyList_GET_SIZE(self->counters);
        for(i=0; i<size;i++){
            Py_XDECREF(PyList_GET_ITEM(self->counters, i));
        }
    }
    */
    Py_XDECREF(self->counters);

    if(self->multimaps){
        size = PyList_GET_SIZE(self->multimaps);
        for(i=0; i<size;i++){
            Py_XDECREF(PyList_GET_ITEM(self->multimaps, i));
        }
    }
    Py_XDECREF(self->multimaps);

    self->ob_type->tp_free((PyObject*)self);
}



static PyObject *
Template_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    Template *self;

    self = (Template *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->tag = PyString_FromString("");
        if (self->tag == NULL)
          {
            Py_DECREF(self);
            return NULL;
          }
        self->alignments = NULL;
        self->counters= NULL;
        self->multimaps= NULL;

        /*
        self->alignments = PyList_New(0);
        if (self->alignments == NULL)
          {
            Py_DECREF(self);
            return NULL;
          }
        self->counters = PyList_New(0);
        if (self->counters == NULL)
          {
            Py_DECREF(self);
            return NULL;
          }
        self->multimaps= PyList_New(0);
        if (self->multimaps == NULL)
          {
            Py_DECREF(self);
            return NULL;
          }
        */

    }
    return (PyObject *)self;
}

static int
Template_init(Template *self, PyObject *args, PyObject *kwds)
{
    PyObject *tag=NULL, *alignments=NULL, *counters=NULL, *multimaps=NULL, *tmp;

    static char *kwlist[] = {"tag", "alignments", "counters","multimaps", NULL};

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "|sOOO", kwlist,
                                      &tag, &alignments, &counters, &multimaps
                                      ))
        return -1;

    if (tag) {
        tmp = self->tag;
        Py_INCREF(tag);
        self->tag = tag;
        Py_XDECREF(tag);
    }

    if (alignments) {
        tmp = self->alignments;
        Py_INCREF(alignments);
        self->alignments = alignments;
        Py_XDECREF(tmp);
    }

    if (counters) {
        tmp = self->counters;
        Py_INCREF(counters);
        self->counters = counters;
        Py_XDECREF(tmp);
    }
    if (multimaps) {
        tmp = self->multimaps;
        Py_INCREF(multimaps);
        self->multimaps = multimaps;
        Py_XDECREF(tmp);
    }

    return 0;
}


static PyMemberDef Template_members[] = {
    {"tag", T_OBJECT_EX, offsetof(Template, tag), 0,
     "Sequence Tag"},
    {"alignments", T_OBJECT_EX, offsetof(Template, alignments), 0,
     "Alignments"},
    {"counters", T_OBJECT_EX, offsetof(Template, counters), 0,
     "Counters"},
    {"multimaps", T_OBJECT_EX, offsetof(Template, multimaps), 0,
     "Multimaps"},
    {NULL}  /* Sentinel */
};


static PyMethodDef Template_methods[] = {
    {NULL}  /* Sentinel */
};

static PyTypeObject TemplateType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "gempy.Template",             /*tp_name*/
    sizeof(Template),             /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)Template_dealloc, /*tp_dealloc*/
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
    Template_methods,             /* tp_methods */
    Template_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Template_init,      /* tp_init */
    0,                         /* tp_alloc */
    Template_new,                 /* tp_new */
};

Template* create_template(){
    Template* tmpl = PyObject_New(Template, &TemplateType);
    tmpl->tag = NULL;
    tmpl->alignments = NULL;
    tmpl->counters = NULL;
    tmpl->multimaps = NULL;
    Py_INCREF(tmpl);
    return tmpl;
}
