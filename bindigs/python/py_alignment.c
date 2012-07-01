
/*
 * Alignment class
 */
#include "Python.h"
#include "structmember.h"


typedef struct {
    PyObject_HEAD
    PyObject *tag; /* The read ID */
    PyObject *qualities;  /* List of alignments*/
    PyObject *read;  /* List of alignments*/
    PyObject *counters; /*List of counters*/
    PyObject *maps; /*List of tuples for maps*/
} Alignment;

static void
Alignment_dealloc(Alignment* self)
{
    register int i = 0;
    register int size = 0;
    Py_XDECREF(self->tag);
    Py_XDECREF(self->read);
    Py_XDECREF(self->qualities);
    /*
    if(self->counters){
        size = PyList_GET_SIZE(self->counters);
        for(i=0; i<size;i++){
            Py_XDECREF(PyList_GET_ITEM(self->counters, i));
        }
    }
    */
    Py_XDECREF(self->counters);
    if(self->maps){
        size = PyList_GET_SIZE(self->maps);
        for(i=0; i<size;i++){
            Py_XDECREF(PyList_GET_ITEM(self->maps, i));
        }
    }
    Py_XDECREF(self->maps);
    self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
Alignment_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    Alignment *self;

    self = (Alignment *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->tag = PyString_FromString("");
        if (self->tag == NULL)
          {
            Py_DECREF(self);
            return NULL;
          }
        self->qualities = PyString_FromString("");
        if (self->qualities == NULL)
          {
            Py_DECREF(self);
            return NULL;
          }
        self->read = PyList_New(0);
        if (self->read == NULL)
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
        self->maps= PyList_New(0);
        if (self->maps == NULL)
          {
            Py_DECREF(self);
            return NULL;
          }

    }
    return (PyObject *)self;
}

static int
Alignment_init(Alignment *self, PyObject *args, PyObject *kwds)
{
    PyObject *tag=NULL, *read=NULL, *qualities=NULL, *counters=NULL, *maps=NULL, *tmp;

    static char *kwlist[] = {"tag", "read", "qualities", "counters","maps", NULL};

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "|sOOO", kwlist,
                                      &tag, &read, &qualities, &counters, &maps
                                      ))
        return -1;

    if (tag) {
        tmp = self->tag;
        Py_INCREF(tag);
        self->tag = tag;
        Py_XDECREF(tag);
    }

    if (read) {
        tmp = self->read;
        Py_INCREF(read);
        self->read = read;
        Py_XDECREF(tmp);
    }

    if (qualities) {
        tmp = self->qualities;
        Py_INCREF(qualities);
        self->qualities = qualities;
        Py_XDECREF(tmp);
    }

    if (counters) {
        tmp = self->counters;
        Py_INCREF(counters);
        self->counters = counters;
        Py_XDECREF(tmp);
    }
    if (maps) {
        tmp = self->maps;
        Py_INCREF(maps);
        self->maps = maps;
        Py_XDECREF(tmp);
    }

    return 0;
}


static PyMemberDef Alignment_members[] = {
    {"tag", T_OBJECT_EX, offsetof(Alignment, tag), 0,
     "Sequence Tag"},
    {"read", T_OBJECT_EX, offsetof(Alignment, read), 0,
     "read"},
    {"qualities", T_OBJECT_EX, offsetof(Alignment, qualities), 0,
     "qualities"},
    {"counters", T_OBJECT_EX, offsetof(Alignment, counters), 0,
     "Counters"},
    {"maps", T_OBJECT_EX, offsetof(Alignment, maps), 0,
     "maps"},
    {NULL}  /* Sentinel */
};


static PyMethodDef Alignment_methods[] = {
    {NULL}  /* Sentinel */
};

static PyTypeObject AlignmentType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "gempy.Alignment",             /*tp_name*/
    sizeof(Alignment),             /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)Alignment_dealloc, /*tp_dealloc*/
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
    Alignment_methods,             /* tp_methods */
    Alignment_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Alignment_init,      /* tp_init */
    0,                         /* tp_alloc */
    Alignment_new,                 /* tp_new */
};


Alignment* create_alignment(){
    Alignment* a = PyObject_New(Alignment, &AlignmentType);
    a->tag = NULL;
    a->qualities = NULL;
    a->read = NULL;
    a->counters = NULL;
    a->maps = NULL;
    Py_INCREF(a);
    return a;
}
