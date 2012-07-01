#include "Python.h"
#include "gem_tools.h"

#include "py_template.c"
#include "py_alignment.c"

typedef struct {
  PyObject_HEAD
  gt_template* template;
  gt_buffered_map_input* map_input;
  PyObject* tmpl;
} spam_MyIter;


PyObject* spam_MyIter_iter(PyObject *self)
{
  Py_INCREF(self);
  return self;
}


PyObject* spam_MyIter_iternext(PyObject *self)
{
    spam_MyIter *p = (spam_MyIter *)self;
    gt_status error_code;
    Template* tmpl = NULL;
    char* tag;
    char* ali_tag;
    char* quals;
    char* read;
    char* qualities;
    gt_template* template = p->template;
    uint64_t alignment_counter = 0;
    uint64_t num_counters = 0;
    uint64_t j = 0;
    uint64_t count = 0;
    PyObject* cvalue = NULL;
    while ((error_code=gt_buffered_map_input_get_template(p->map_input,p->template))) {
        if (error_code==GT_BMI_FAIL) continue;
        GT_TEMPLATE_ITERATE(template, map_array) {
            tmpl = create_template();
            // fill python template
            tag = gt_template_get_tag(template);
            tmpl->tag = PyString_FromString(tag);
            // createa counter list
            num_counters = gt_template_get_num_counters(template);
            tmpl->counters = PyList_New(num_counters);
            for(j=0;j<num_counters;j++){
                count = gt_template_get_counter(template, j+1);
                cvalue = PyLong_FromUnsignedLongLong(count);
                PyList_SET_ITEM(tmpl->counters, j, cvalue);
            }

            // createa alignment list
            register const uint64_t template_num_blocks = gt_template_get_num_blocks(template);
            tmpl->alignments = PyList_New(template_num_blocks);
            // fill arrays
            GT_ALIGNMENT_ITERATE(template, alignment){
                register Alignment* ali = create_alignment();
                PyList_SET_ITEM(tmpl->alignments, alignment_counter++, (PyObject*)ali);
                ali_tag = gt_alignment_get_tag(alignment);
                ali->tag = PyString_FromString(ali_tag ? ali_tag : tag);
                ali->read = PyString_FromString(gt_alignment_get_read(alignment));
                qualities= gt_alignment_get_qualities(alignment);
                if(qualities){
                    ali->qualities = PyString_FromString(qualities);
                }
                // fill counters
                num_counters = gt_alignment_get_num_counters(alignment);
                ali->counters = PyList_New(num_counters);
                for(j=0;j<num_counters;j++){
                    count = gt_alignment_get_counter(alignment, j);
                    cvalue = PyLong_FromUnsignedLongLong(count);
                    PyList_SET_ITEM(ali->counters, j, cvalue);
                }
                printf("Counter list applied %s %d\n", tag, num_counters);
            }
            if(p->tmpl){
                Py_CLEAR(p->tmpl);
            }
            p->tmpl = (PyObject*) tmpl;
            return (PyObject*) tmpl;
        }
  }
  /* Raising of standard StopIteration exception with empty value. */
  PyErr_SetNone(PyExc_StopIteration);
  return NULL;
}


static PyTypeObject spam_MyIterType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "spam._MyIter",            /*tp_name*/
    sizeof(spam_MyIter),       /*tp_basicsize*/
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
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER,
      /* tp_flags: Py_TPFLAGS_HAVE_ITER tells python to
         use tp_iter and tp_iternext fields. */
    "Internal myiter iterator object.",           /* tp_doc */
    0,  /* tp_traverse */
    0,  /* tp_clear */
    0,  /* tp_richcompare */
    0,  /* tp_weaklistoffset */
    spam_MyIter_iter,  /* tp_iter: __iter__() method */
    spam_MyIter_iternext  /* tp_iternext: next() method */
};


static PyObject * spam_myiter(PyObject *self, PyObject *args){
    void* m;
    spam_MyIter *p;

    if (!PyArg_ParseTuple(args, "O", &m))  return NULL;

    /* I don't need python callable __init__() method for this iterator,
         so I'll simply allocate it as PyObject and initialize it by hand. */

    p = PyObject_New(spam_MyIter, &spam_MyIterType);
    if (!p) return NULL;

    /* I'm not sure if it's strictly necessary. */
    if (!PyObject_Init((PyObject *)p, &spam_MyIterType)) {
        Py_DECREF(p);
        return NULL;
    }
    FILE* file = PyFile_AsFile(m);
    p->map_input = gt_buffered_map_input_new(gt_input_stream_open(file));
    p->template = gt_template_new();
    p->tmpl = NULL;
    return (PyObject *)p;
}


static PyMethodDef GempyMethods[] = {
    {"myiter",  spam_myiter, METH_VARARGS, "Iterate from i=0 while i<m."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initgempy(void)
{
  PyObject* m;

  spam_MyIterType.tp_new = PyType_GenericNew;
  if (PyType_Ready(&spam_MyIterType) < 0)  return;
  if (PyType_Ready(&TemplateType) < 0)  return;
  if (PyType_Ready(&AlignmentType) < 0)  return;

  m = Py_InitModule("gempy", GempyMethods);

  Py_INCREF(&spam_MyIterType);
  Py_INCREF(&TemplateType);
  Py_INCREF(&AlignmentType);
  PyModule_AddObject(m, "_MyIter", (PyObject *)&spam_MyIterType);
  PyModule_AddObject(m, "Template", (PyObject *)&TemplateType);
  PyModule_AddObject(m, "Alignment", (PyObject *)&AlignmentType);
}



/*
static PyObject * gempy_testme(PyObject* self, PyObject* args){
    return Py_BuildValue("s", "lala lulu");
}

static PyMethodDef GempyMethods[] = {
    {"testme",  gempy_testme, METH_VARARGS,"Test me method."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initgempy(void){
    (void) Py_InitModule("gempy", GempyMethods);
}
*/
/*
PyObject* gtpy_fill_template(PyObject *py_tmpl, gt_template* gt_tmpl, gt_buffered_map_input* map_input){
    gt_status error_code;
    int py_err;
    char* tag;
    while ((error_code=gt_buffered_map_input_get_template(map_input,gt_tmpl))) {
        if (error_code==GT_BMI_FAIL) continue;
        tag = gt_template_get_tag(gt_tmpl);
        py_err = PyObject_SetAttrString(py_tmpl, "tag", PyString_FromString(tag));
        printf("%d %s", py_err, tag);
        //register const uint64_t num_blocks = gt_template_get_num_blocks(gt_tmpl);
        //register uint64_t i;
        //for (i=0;i<num_blocks;++i) {
        //    register gt_alignment* alignment = gt_template_get_block(gt_tmpl,0);

        //    //return py_tmpl;
        //    //printf("@%s\n%s\n+\n%s\n",gt_template_get_tag(gt_tmpl),
        //    //    gt_alignment_get_read(alignment),gt_alignment_get_qualities(alignment));
        //}
        break;
    }
    if(error_code == 0){
        return Py_None;
    }else{
        return py_tmpl;
    }
}
*/
