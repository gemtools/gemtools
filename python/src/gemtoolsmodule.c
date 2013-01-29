#include "Python.h"
#include "gem_tools.h"
//#include "py_mismatch.h"
//#include "py_map.h"
//#include "py_alignment.h"
#include "py_template.h"
//#include "py_template_iterator.h"
//#include "py_mappings_iterator.h"
//#include "py_iterator.h"

/******* TEMPLATE ********/
//static PyObject* Template_iterate_mappings(PyObject* self, PyObject* closure){
//    Template* t = (Template*)self;
//    return (PyObject*) create_template_mappings_iterator(t->template);
//}
//
static PyGetSetDef Template_getseters[] = {
//    {"tag", (getter) Template_gettag, (setter) Template_settag, "Template Tag", NULL},
    {"max_complete_strata", (getter) Template_get_max_complete_strata, NULL, "Max Complete Strata", NULL},
//    {"num_blocks", (getter) Template_get_num_blocks, NULL, "Get the number of blocks in the template", NULL},
//    {"num_counters", (getter) Template_get_num_blocks, NULL, "Get the number of counters for the template", NULL},
    {NULL}  /* Sentinel */
};
//
static PyMethodDef Template_Methods[] = {
    {"fill", Template_fill, METH_VARARGS, "Fill the template"},
//    {"blocks", Template_get_blocks, METH_VARARGS, "Iterator over the alignment blocks in the template"},
//    {"counters", Template_get_counters, METH_VARARGS, "Iterator over the counter list"},
//    {"mappings", Template_iterate_mappings, METH_VARARGS, "Iterate over the mappings"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


static PyTypeObject TemplateType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "gempy.Template",             /*tp_name*/
    sizeof(Template),             /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    Template_dealloc,                         /*tp_dealloc*/
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
    "Represents a gemtools gt_template",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    Template_Methods,             /* tp_methods */
    0,             /* tp_members */
    Template_getseters,        /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Template_init,      /* tp_init */
    0,                         /* tp_alloc */
    Template_new,                 /* tp_new */
};

Template* create_template(gt_template* template){
    Template* tmpl = PyObject_New(Template, &TemplateType);
    tmpl->template = template;
    return tmpl;
}
///******* END TEMPLATE ********/
//
//
///******* TEMPLATE ITERATOR ********/
//static PyTypeObject gempy_template_iteratorType = {
//    PyObject_HEAD_INIT(NULL)
//    0,                         /*ob_size*/
//    "gem.gemtools._template_iterator",            /*tp_name*/
//    sizeof(gempy_template_iterator),       /*tp_basicsize*/
//    0,                         /*tp_itemsize*/
//    gempy_template_iterator_dealloc, /*tp_dealloc*/
//    0,                         /*tp_print*/
//    0,                         /*tp_getattr*/
//    0,                         /*tp_setattr*/
//    0,                         /*tp_compare*/
//    0,                         /*tp_repr*/
//    0,                         /*tp_as_number*/
//    0,                         /*tp_as_sequence*/
//    0,                         /*tp_as_mapping*/
//    0,                         /*tp_hash */
//    0,                         /*tp_call*/
//    0,                         /*tp_str*/
//    0,                         /*tp_getattro*/
//    0,                         /*tp_setattro*/
//    0,                         /*tp_as_buffer*/
//    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER,
//      /* tp_flags: Py_TPFLAGS_HAVE_ITER tells python to
//         use tp_iter and tp_iternext fields. */
//    "Internal Template iterator",           /* tp_doc */
//    0,  /* tp_traverse */
//    0,  /* tp_clear */
//    0,  /* tp_richcompare */
//    0,  /* tp_weaklistoffset */
//    gempy_template_iterator_iter,  /* tp_iter: __iter__() method */
//    gempy_template_iterator_iternext  /* tp_iternext: next() method */
//};
//
//gempy_template_iterator* create_template_stream_iterator(FILE* file){
//    gempy_template_iterator *p;
//    p = PyObject_New(gempy_template_iterator, &gempy_template_iteratorType);
//    if (!p) return NULL;
//    p->input_file = gt_input_stream_open(file);
//    p->map_input = gt_buffered_input_file_new(p->input_file);
//    p->template = gt_template_new();
//    p->tmpl = NULL;
//    return (gempy_template_iterator *)p;
//}
//
//gempy_template_iterator* create_template_file_iterator(char* filename, bool memorymap){
//    gempy_template_iterator *p;
//    p = PyObject_New(gempy_template_iterator, &gempy_template_iteratorType);
//    if (!p) return NULL;
//    /* I'm not sure if it's strictly necessary. */
//    if (!PyObject_Init((PyObject *)p, &gempy_template_iteratorType)) {
//        printf("ERROR template iterator init failed!\n");
//        Py_DECREF(p);
//        return NULL;
//    }
//    p->input_file = gt_input_file_open(filename, memorymap);
//    p->map_input = gt_buffered_input_file_new(p->input_file); // false disable memory map
//    p->template = gt_template_new();
//    p->tmpl = NULL;
//    return (gempy_template_iterator *)p;
//}
///******* END TEMPLATE ITERATOR ********/
//
//
//
///******* ALIGNMENT  ********/
//
//static PyObject* Alignment_iterate_mappings(PyObject* self, PyObject* closure){
//    Alignment* t = (Alignment*)self;
//    gempy_alignment_mappings_iterator* ret =  create_alignment_mappings_iterator(t->alignment);
//    ret->template = t->template;
//    return (PyObject*) ret;
//}
//
//static PyMethodDef Alignmnt_methods[] = {
//    {"to_sequence", Alignment_to_sequence, METH_VARARGS, "Convert alignment Alignment_to_sequence"},
//    {"mappings", Alignment_iterate_mappings, METH_VARARGS, "Iterator over all mappings fir this alignment"},
//    {"counters", Alignment_get_counters, METH_VARARGS, "Counters"},
//    {NULL, NULL, 0, NULL}        /* Sentinel */
//};
//
//static PyGetSetDef Alignment_getseters[] = {
//    {"tag", (getter) Alignment_gettag, (setter) Alignment_settag, "Alignment Tag", NULL},
//    {"read", (getter) Alignment_getread, (setter) Alignment_setread, "Alignment Read", NULL},
//    {"qualities", (getter) Alignment_getqualities, (setter) Alignment_setqualities, "Alignment Qualities", NULL},
//    {"max_complete_strata", (getter) Alignment_getmax_complete_strata,
//        (setter) Alignment_setmax_complete_strata, "Max Complete Strata", NULL},
//    {NULL}  /* Sentinel */
//};
//
//static PyTypeObject AlignmentType = {
//    PyObject_HEAD_INIT(NULL)
//    0,                         /*ob_size*/
//    "gempy.Alignment",             /*tp_name*/
//    sizeof(Alignment),             /*tp_basicsize*/
//    0,                         /*tp_itemsize*/
//    Alignment_dealloc, /*tp_dealloc*/
//    0,                         /*tp_print*/
//    0,                         /*tp_getattr*/
//    0,                         /*tp_setattr*/
//    0,                         /*tp_compare*/
//    0,                         /*tp_repr*/
//    0,                         /*tp_as_number*/
//    0,                         /*tp_as_sequence*/
//    0,                         /*tp_as_mapping*/
//    0,                         /*tp_hash */
//    0,                         /*tp_call*/
//    0,                         /*tp_str*/
//    0,                         /*tp_getattro*/
//    0,                         /*tp_setattro*/
//    0,                         /*tp_as_buffer*/
//    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
//    "Alignment objects",           /* tp_doc */
//    0,		               /* tp_traverse */
//    0,		               /* tp_clear */
//    0,		               /* tp_richcompare */
//    0,		               /* tp_weaklistoffset */
//    0,		               /* tp_iter */
//    0,		               /* tp_iternext */
//    Alignmnt_methods,             /* tp_methods */
//    0,             /* tp_members */
//    Alignment_getseters,                         /* tp_getset */
//    0,                         /* tp_base */
//    0,                         /* tp_dict */
//    0,                         /* tp_descr_get */
//    0,                         /* tp_descr_set */
//    0,                         /* tp_dictoffset */
//    (initproc)Alignment_init,      /* tp_init */
//    0,                         /* tp_alloc */
//    Alignment_new,                 /* tp_new */
//};
//
//
//Alignment* create_alignment(gt_alignment* alignment, gt_template* parent, uint64_t index){
//    Alignment* a = PyObject_New(Alignment, &AlignmentType);
//    a->alignment = alignment;
//    a->template = parent;
//    a->index = index;
//    return a;
//}
///******* END ALIGNMENT  ********/
//
//
//
///******* MAP  ********/
//static PyGetSetDef Map_getseters[] = {
//    {"seq_name", (getter) Map_getseq_name, (setter) Map_setseq_name, "Genomic sequence name", NULL},
//    {"position", (getter) Map_getposition, (setter) Map_setposition, "Genomic Position", NULL},
//    {"base_length", (getter) Map_getbase_length, (setter) Map_setbase_length, "Base length without indels. This is the length of the read.", NULL},
//    {"length", (getter) Map_getlength, NULL, "Length of the mapping including indels", NULL},
//    {"direction", (getter) Map_getdirection, (setter) Map_setdirection, "Strand", NULL},
//    {"score", (getter) Map_getscore, (setter) Map_setscore, "Score", NULL},
//    {"distance", (getter) Map_getdistance, NULL, "Map distance", NULL},
//    {"levenshtein", (getter) Map_getlevenshtein, NULL, "Map levenshtein", NULL},
//    {"global_length", (getter) Map_getglobal_length, NULL, "Map global_length", NULL},
//    {"global_score", (getter) Map_getglobal_score, NULL, "Map global_score", NULL},
//    {"global_distance", (getter) Map_getglobal_distance, NULL, "Map global_distance", NULL},
//    {"global_levenshtein", (getter) Map_getglobal_levenshtein, NULL, "Map global_levenshtein", NULL},
//    {"num_mismatches", (getter) Map_getnum_mismatches, NULL, "Mismatches", NULL},
//    {NULL}  /* Sentinel */
//};
//
//static PyMethodDef Map_methods[] = {
//    {"mismatches", Map_get_mismatches, METH_VARARGS, "Iterator over the mismatches"},
//    {NULL, NULL, 0, NULL}        /* Sentinel */
//};
//
//static PyTypeObject MapType = {
//    PyObject_HEAD_INIT(NULL)
//    0,                         /*ob_size*/
//    "gempy.Map",             /*tp_name*/
//    sizeof(Map),             /*tp_basicsize*/
//    0,                         /*tp_itemsize*/
//    Map_dealloc, /*tp_dealloc*/
//    0,                         /*tp_print*/
//    0,                         /*tp_getattr*/
//    0,                         /*tp_setattr*/
//    0,                         /*tp_compare*/
//    0,                         /*tp_repr*/
//    0,                         /*tp_as_number*/
//    0,                         /*tp_as_sequence*/
//    0,                         /*tp_as_mapping*/
//    0,                         /*tp_hash */
//    0,                         /*tp_call*/
//    0,                         /*tp_str*/
//    0,                         /*tp_getattro*/
//    0,                         /*tp_setattro*/
//    0,                         /*tp_as_buffer*/
//    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
//    "Map objects",           /* tp_doc */
//    0,		               /* tp_traverse */
//    0,		               /* tp_clear */
//    0,		               /* tp_richcompare */
//    0,		               /* tp_weaklistoffset */
//    0,		               /* tp_iter */
//    0,		               /* tp_iternext */
//    Map_methods,             /* tp_methods */
//    0,             /* tp_members */
//    Map_getseters,                         /* tp_getset */
//    0,                         /* tp_base */
//    0,                         /* tp_dict */
//    0,                         /* tp_descr_get */
//    0,                         /* tp_descr_set */
//    0,                         /* tp_dictoffset */
//    (initproc)Map_init,      /* tp_init */
//    0,                         /* tp_alloc */
//    Map_new,                 /* tp_new */
//};
//
//Map* create_map(gt_map* map){
//    Map* a = PyObject_New(Map, &MapType);
//    a->map = map;
//    return a;
//}
///******* END MAP  ********/
//
//
///******* MISMATCH  ********/
//static PyGetSetDef Mismatch_getseters[] = {
//    {"type", (getter) Mismatch_gettype, (setter) Mismatch_settype, "Set type", NULL},
//    {"position", (getter) Mismatch_getposition, (setter) Mismatch_setposition, "Set position", NULL},
//    {"base", (getter) Mismatch_getbase, (setter) Mismatch_setbase, "Set base", NULL},
//    {"size", (getter) Mismatch_getsize, (setter) Mismatch_setsize, "Set size", NULL},
//    {NULL}  /* Sentinel */
//};
//static PyTypeObject MismatchType = {
//    PyObject_HEAD_INIT(NULL)
//    0,                         /*ob_size*/
//    "gempy.Mismatch",          /*tp_name*/
//    sizeof(Mismatch),          /*tp_basicsize*/
//    0,                         /*tp_itemsize*/
//    Mismatch_dealloc,                         /*tp_dealloc*/
//    0,                         /*tp_print*/
//    0,                         /*tp_getattr*/
//    0,                         /*tp_setattr*/
//    0,                         /*tp_compare*/
//    0,                         /*tp_repr*/
//    0,                         /*tp_as_number*/
//    0,                         /*tp_as_sequence*/
//    0,                         /*tp_as_mapping*/
//    0,                         /*tp_hash */
//    0,                         /*tp_call*/
//    0,                         /*tp_str*/
//    0,                         /*tp_getattro*/
//    0,                         /*tp_setattro*/
//    0,                         /*tp_as_buffer*/
//    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
//    "Mismatch objects",       /* tp_doc */
//    0,		               /* tp_traverse */
//    0,		               /* tp_clear */
//    0,		               /* tp_richcompare */
//    0,		               /* tp_weaklistoffset */
//    0,		               /* tp_iter */
//    0,		               /* tp_iternext */
//    0,                     /* tp_methods */
//    0,                     /* tp_members */
//    Mismatch_getseters,    /* tp_getset */
//    0,                         /* tp_base */
//    0,                         /* tp_dict */
//    0,                         /* tp_descr_get */
//    0,                         /* tp_descr_set */
//    0,                         /* tp_dictoffset */
//    (initproc)Mismatch_init,      /* tp_init */
//    0,                         /* tp_alloc */
//    Mismatch_new,                 /* tp_new */
//};
//
//Mismatch* create_mismatch(gt_misms* map){
//    Mismatch* a = PyObject_New(Mismatch, &MismatchType);
//    a->misms = map;
//    return a;
//}
///******* END MISMATCH  ********/
//
//
///******* GENERIC ITERATOR ********/
//
//
//static PySequenceMethods gempy_iterator_sequence_methods = {
//    gempy_iterator_len, /* sq_length */
//};
//
//static PyTypeObject gempy_iteratorType = {
//    PyObject_HEAD_INIT(NULL)
//    0,                         /*ob_size*/
//    "gempy._iterator",            /*tp_name*/
//    sizeof(gempy_iterator),       /*tp_basicsize*/
//    0,                         /*tp_itemsize*/
//    gempy_iterator_dealloc,                         /*tp_dealloc*/
//    0,                         /*tp_print*/
//    0,                         /*tp_getattr*/
//    0,                         /*tp_setattr*/
//    0,                         /*tp_compare*/
//    0,                         /*tp_repr*/
//    0,                         /*tp_as_number*/
//    &gempy_iterator_sequence_methods, /*tp_as_sequence*/
//    0,                         /*tp_as_mapping*/
//    0,                         /*tp_hash */
//    0,                         /*tp_call*/
//    0,                         /*tp_str*/
//    0,                         /*tp_getattro*/
//    0,                         /*tp_setattro*/
//    0,                         /*tp_as_buffer*/
//    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER,
//      /* tp_flags: Py_TPFLAGS_HAVE_ITER tells python to
//         use tp_iter and tp_iternext fields. */
//    "Internal iterator object.",           /* tp_doc */
//    0,  /* tp_traverse */
//    0,  /* tp_clear */
//    0,  /* tp_richcompare */
//    0,  /* tp_weaklistoffset */
//    gempy_iterator_iter,  /* tp_iter: __iter__() method */
//    gempy_iterator_iternext  /* tp_iternext: next() method */
//};
//
//
//PyObject* create_gempy_iterator(uint64_t start, uint64_t length, void* getter, void* arg, void* converter, int parent){
//    gempy_iterator *p;
//    p = PyObject_New(gempy_iterator, &gempy_iteratorType);
//    if (!p) return NULL;
//    p->start = start;
//    p->length = length;
//    p->pos = 0;
//    p->getter = getter;
//    p->arg = arg;
//    p->converter = converter;
//    p->parent = parent;
//    return (PyObject *)p;
//}
///******* END GENERIC ITERATOR  ********/
//
///******* MAPPINGS ITERATOR ********/
//static PyMemberDef gempy_mappings_iterator_members[] = {
//    {"distance", T_ULONGLONG, offsetof(gempy_mappings_iterator, distance), 1, "Edit distance of the mapping"},
//    {"score", T_ULONGLONG, offsetof(gempy_mappings_iterator, score), 1, "Score of the mapping"},
//    {NULL}
//};
//
//static PyTypeObject gempy_mappings_iteratorType = {
//    PyObject_HEAD_INIT(NULL)
//    0,                         /*ob_size*/
//    "gem.gemtools._mappings_iterator",            /*tp_name*/
//    sizeof(gempy_mappings_iterator),       /*tp_basicsize*/
//    0,                         /*tp_itemsize*/
//    gempy_mappings_iterator_dealloc, /*tp_dealloc*/
//    0,                         /*tp_print*/
//    0,                         /*tp_getattr*/
//    0,                         /*tp_setattr*/
//    0,                         /*tp_compare*/
//    0,                         /*tp_repr*/
//    0,                         /*tp_as_number*/
//    0,                         /*tp_as_sequence*/
//    0,                         /*tp_as_mapping*/
//    0,                         /*tp_hash */
//    0,                         /*tp_call*/
//    0,                         /*tp_str*/
//    0,                         /*tp_getattro*/
//    0,                         /*tp_setattro*/
//    0,                         /*tp_as_buffer*/
//    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER,
//      /* tp_flags: Py_TPFLAGS_HAVE_ITER tells python to
//         use tp_iter and tp_iternext fields. */
//    "Internal Mappings iterator",           /* tp_doc */
//    0,  /* tp_traverse */
//    0,  /* tp_clear */
//    0,  /* tp_richcompare */
//    0,  /* tp_weaklistoffset */
//    gempy_mappings_iterator_iter,  /* tp_iter: __iter__() method */
//    gempy_mappings_iterator_iternext,  /* tp_iternext: next() method */
//    0,             /* tp_methods */
//    gempy_mappings_iterator_members,             /* tp_members */
//    0,        /* tp_getset */
//};
//
//gempy_mappings_iterator* create_mappings_iterator(gt_map* map_block){
//    gempy_mappings_iterator *p;
//    p = PyObject_New(gempy_mappings_iterator, &gempy_mappings_iteratorType);
//    if (!p) return NULL;
//    if (!PyObject_Init((PyObject *)p, &gempy_mappings_iteratorType)) {
//        printf("ERROR mappings iterator init failed!\n");
//        Py_DECREF(p);
//        return NULL;
//    }
//    p->map_block = map_block;
//    p->alignment_iterator = NULL;
//    p->distance = 0;
//    p->score = 0;
//    //printf("Created mappings_iterator %p\n", p);
//    return p;
//}
//
///******* END MAPPINGS ITERATOR ********/
//
///******* ALIGNMENTS ITERATOR ********/
//static PyTypeObject gempy_alignment_iteratorType = {
//    PyObject_HEAD_INIT(NULL)
//    0,                         /*ob_size*/
//    "gem.gemtools._alignment_iterator",            /*tp_name*/
//    sizeof(gempy_alignment_iterator),       /*tp_basicsize*/
//    0,                         /*tp_itemsize*/
//    gempy_alignment_iterator_dealloc, /*tp_dealloc*/
//    0,                         /*tp_print*/
//    0,                         /*tp_getattr*/
//    0,                         /*tp_setattr*/
//    0,                         /*tp_compare*/
//    0,                         /*tp_repr*/
//    0,                         /*tp_as_number*/
//    0,                         /*tp_as_sequence*/
//    0,                         /*tp_as_mapping*/
//    0,                         /*tp_hash */
//    0,                         /*tp_call*/
//    0,                         /*tp_str*/
//    0,                         /*tp_getattro*/
//    0,                         /*tp_setattro*/
//    0,                         /*tp_as_buffer*/
//    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER,
//      /* tp_flags: Py_TPFLAGS_HAVE_ITER tells python to
//         use tp_iter and tp_iternext fields. */
//    "Internal Mappings iterator",           /* tp_doc */
//    0,  /* tp_traverse */
//    0,  /* tp_clear */
//    0,  /* tp_richcompare */
//    0,  /* tp_weaklistoffset */
//    gempy_alignment_iterator_iter,  /* tp_iter: __iter__() method */
//    gempy_alignment_iterator_iternext  /* tp_iternext: next() method */
//};
//
//
//static PyTypeObject gempy_alignment_mappings_iteratorType = {
//    PyObject_HEAD_INIT(NULL)
//    0,                         /*ob_size*/
//    "gem.gemtools._alignment_mappings_iterator",            /*tp_name*/
//    sizeof(gempy_alignment_mappings_iterator),       /*tp_basicsize*/
//    0,                         /*tp_itemsize*/
//    gempy_alignment_mappings_iterator_dealloc, /*tp_dealloc*/
//    0,                         /*tp_print*/
//    0,                         /*tp_getattr*/
//    0,                         /*tp_setattr*/
//    0,                         /*tp_compare*/
//    0,                         /*tp_repr*/
//    0,                         /*tp_as_number*/
//    0,                         /*tp_as_sequence*/
//    0,                         /*tp_as_mapping*/
//    0,                         /*tp_hash */
//    0,                         /*tp_call*/
//    0,                         /*tp_str*/
//    0,                         /*tp_getattro*/
//    0,                         /*tp_setattro*/
//    0,                         /*tp_as_buffer*/
//    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER,
//      /* tp_flags: Py_TPFLAGS_HAVE_ITER tells python to
//         use tp_iter and tp_iternext fields. */
//    "Internal Alignment Mappings iterator",           /* tp_doc */
//    0,  /* tp_traverse */
//    0,  /* tp_clear */
//    0,  /* tp_richcompare */
//    0,  /* tp_weaklistoffset */
//    gempy_alignment_mappings_iterator_iter,  /* tp_iter: __iter__() method */
//    gempy_alignment_mappings_iterator_iternext  /* tp_iternext: next() method */
//};
//
//gempy_alignment_iterator* create_template_mappings_iterator(gt_template* template){
//    gempy_alignment_iterator *p;
//    p = PyObject_New(gempy_alignment_iterator, &gempy_alignment_iteratorType);
//    if (!p) return NULL;
//    if (!PyObject_Init((PyObject *)p, &gempy_alignment_iteratorType)) {
//        printf("ERROR alignment iterator init failed!\n");
//        Py_DECREF(p);
//        return NULL;
//    }
//    p->template = template;
//    p->map_array = NULL;
//    p->jump=-1;
//    p->num_blocks = gt_template_get_num_blocks(template);
//    p->position = 0;
//    //gt_template_new_maps_iterator(template,&(p->maps_iterator));
//    return p;
//}
//
//gempy_alignment_mappings_iterator* create_alignment_mappings_iterator(gt_alignment* alignment){
//    gempy_alignment_mappings_iterator *p;
//    p = PyObject_New(gempy_alignment_mappings_iterator, &gempy_alignment_mappings_iteratorType);
//    if (!p) return NULL;
//    if (!PyObject_Init((PyObject *)p, &gempy_alignment_mappings_iteratorType)) {
//        printf("ERROR alignment mappings iterator init failed!\n");
//        Py_DECREF(p);
//        return NULL;
//    }
//    p->alignment = alignment;
//    p->total = gt_alignment_get_num_maps(alignment);
//    p->current = 0;
//    p->template = NULL;
//    return p;
//}
//
//
///******* END MAPPINGS ITERATOR ********/
//
//
//
///**** MODULE METHODS AND DEFINITIONS ******/
//
//static PyObject* gempy_open_stream(PyObject *self, PyObject *args){
//    void* m;
//    gempy_template_iterator *p;
//    if (!PyArg_ParseTuple(args, "O", &m))  return NULL;
//    p = create_template_stream_iterator(PyFile_AsFile(m));
//    if (!p) return NULL;
//    Py_DECREF(m);
//    return (PyObject *)p;
//};
//
//static PyObject* gempy_open_file(PyObject *self, PyObject *args){
//    char* filename;
//    PyObject* pp;
//    gempy_template_iterator *p;
//    if (!PyArg_ParseTuple(args, "s", &filename))  return NULL;
//    p = create_template_file_iterator(filename, false);
//    if (!p) return NULL;
//    pp = (PyObject *)p;
//    return pp;
//};

/**
Pass two Template instances and merge the templates. Returns the
merged Template
*/
static PyObject* gemtools_merge_templates(PyObject *self, PyObject *args){
    PyObject *t1 = NULL;
    PyObject *t2 = NULL;
    if (!PyArg_UnpackTuple(args, "ref", 2, 2, &t1, &t2)) {
        PyErr_SetString(PyExc_Exception, "Error while parsing arguments!");
        return NULL;
    }

    Template* t_1 = (Template*)t1;
    Template* t_2 = (Template*)t2;
    register gt_template *ptemplate;

    //gt_template_recalculate_counters(t_1->template);
    //gt_template* tmpl = gt_template_copy(t_1->template, true, false);
    //gt_template_merge_template_mmaps(t_1->template, t_2->template);

    //gt_output_map_fprint_template(stdout, t_1->template, GT_ALL, false);
    //gt_output_map_fprint_template(stdout, t_2->template, GT_ALL, false);

    ptemplate = gt_template_union_template_mmaps_fx(gt_mmap_cmp,gt_map_cmp,t_1->template,t_2->template);


    gt_string* string = gt_string_new(1024);
    gt_output_map_sprint_template(string, ptemplate, GT_ALL, true);
    char * line = gt_string_get_string(string);
    PyObject* returnvalue = PyString_FromString(line);
    gt_string_delete(string);
    gt_template_delete(ptemplate);
    return returnvalue;
};



GT_INLINE gt_status gt_mapset_read_template_sync(
    gt_buffered_input_file* const buffered_input_master,gt_buffered_input_file* const buffered_input_slave,
    gt_template* const template_master,gt_template* const template_slave, bool paired_end) {
  // Read master
  register gt_status error_code_master, error_code_slave;
  gt_generic_parser_attr generic_parser_attr = GENERIC_PARSER_ATTR_DEFAULT(paired_end);
  if ((error_code_master=gt_input_generic_parser_get_template(
      buffered_input_master,template_master,&generic_parser_attr))==GT_IMP_FAIL) {
    gt_fatal_error_msg("Fatal error parsing file <<Master>>");
  }
  // Read slave
  if ((error_code_slave=gt_input_generic_parser_get_template(
      buffered_input_slave,template_slave,&generic_parser_attr))==GT_IMP_FAIL) {
    gt_fatal_error_msg("Fatal error parsing file <<Slave>>");
  }
  // Check EOF conditions
  if (error_code_master==GT_IMP_EOF) {
    if (error_code_slave!=GT_IMP_EOF) {
      gt_fatal_error_msg("<<Slave>> contains more/different reads from <<Master>>");
    }
    return GT_IMP_EOF;
  } else if (error_code_slave==GT_IMP_EOF) { // Slave exhausted. Dump master & return EOF
    do {
      if (error_code_master==GT_IMP_FAIL) gt_fatal_error_msg("Fatal error parsing file <<Master>>");
      gt_output_map_fprint_gem_template(stdout,template_master,GT_ALL,true);
    } while ((error_code_master=gt_input_generic_parser_get_template(
                buffered_input_master,template_master,&generic_parser_attr)));
    return GT_IMP_EOF;
  }
  // Synch loop
  while (!gt_streq(gt_template_get_tag(template_master),gt_template_get_tag(template_slave))) {
    // Print non correlative master's template
    gt_output_map_fprint_gem_template(stdout,template_master,GT_ALL,true);
    // Fetch next master's template
    if ((error_code_master=gt_input_generic_parser_get_template(
        buffered_input_master,template_master,&generic_parser_attr))!=GT_IMP_OK) {
      gt_fatal_error_msg("<<Slave>> contains more/different reads from <<Master>>");
    }
  }
  return GT_IMP_OK;
}

static PyObject* gemtools_merge_files(PyObject *self, PyObject *args){
    PyObject *py_instream1 = NULL;
    PyObject *py_instream2 = NULL;
    PyObject *py_outstream = NULL;
    PyObject *py_paired = NULL;
    if (!PyArg_UnpackTuple(args, "ref", 4, 4, &py_instream1, &py_instream2, &py_outstream, &py_paired)) {
        PyErr_SetString(PyExc_Exception, "Error while parsing arguments!");
        return NULL;
    }

    bool paired = PyObject_IsTrue(py_paired);
    // Open file IN/OUT
    gt_input_file* input_file_1 = gt_input_stream_open(PyFile_AsFile(py_instream1));
    gt_input_file* input_file_2 = gt_input_stream_open(PyFile_AsFile(py_instream2));
    //gt_output_file* output_file = gt_output_stream_new(PyFile_AsFile(py_outstream),SORTED_FILE);
    FILE* output_file = PyFile_AsFile(py_outstream);
    // Parallel reading+process
    gt_buffered_input_file* buffered_input_1 = gt_buffered_input_file_new(input_file_1);
    gt_buffered_input_file* buffered_input_2 = gt_buffered_input_file_new(input_file_2);

    gt_template *template_1 = gt_template_new();
    gt_template *template_2 = gt_template_new();
    while (gt_mapset_read_template_sync(buffered_input_1,buffered_input_2,template_1,template_2,paired)) {
        // Record current read length
        //current_read_length = gt_template_get_total_length(template_1);
        // Apply operation
        register gt_template *ptemplate;
        ptemplate=gt_template_union_template_mmaps(template_1,template_2);
        // Print template
        gt_output_map_fprint_gem_template(output_file, ptemplate, GT_ALL, true);
        // Delete template
        gt_template_delete(ptemplate);
    }
    // Clean
    gt_template_delete(template_1);
    gt_template_delete(template_2);
    gt_buffered_input_file_close(buffered_input_1);
    gt_buffered_input_file_close(buffered_input_2);

    // Clean
    //gt_input_file_close(input_file_1);
    //gt_input_file_close(input_file_2);
    //gt_output_file_close(output_file);
    //fclose(output_file);
    //Py_RETURN_NONE;
    return Py_BuildValue("");
}



static PyMethodDef GempyMethods[] = {
    {"merge_templates", gemtools_merge_templates, METH_VARARGS, "Merge two Templates and returns the merged result Template"},
    {"merge_files", gemtools_merge_files, METH_VARARGS, "Merge two streams of mappings into one output stream"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC initgemtools(void){
  PyObject* m;

//  gempy_iteratorType.tp_new = PyType_GenericNew;
//  gempy_template_iteratorType.tp_new = PyType_GenericNew;
//  gempy_mappings_iteratorType.tp_new = PyType_GenericNew;
//  gempy_alignment_iteratorType.tp_new = PyType_GenericNew;
//  gempy_alignment_mappings_iteratorType.tp_new = PyType_GenericNew;
//  if (PyType_Ready(&gempy_iteratorType) < 0)  return;
//  if (PyType_Ready(&gempy_template_iteratorType) < 0)  return;
//  if (PyType_Ready(&MismatchType) < 0)  return;
//  if (PyType_Ready(&MapType) < 0)  return;
//  if (PyType_Ready(&AlignmentType) < 0)  return;
  if (PyType_Ready(&TemplateType) < 0)  return;
//  if (PyType_Ready(&gempy_mappings_iteratorType) < 0)  return;
//  if (PyType_Ready(&gempy_alignment_iteratorType) < 0)  return;
//  if (PyType_Ready(&gempy_alignment_mappings_iteratorType) < 0)  return;

  m = Py_InitModule("gemtools", GempyMethods);

//  Py_INCREF(&gempy_iteratorType);
//  Py_INCREF(&gempy_template_iteratorType);
//  Py_INCREF(&gempy_mappings_iteratorType);
  Py_INCREF(&TemplateType);
//  Py_INCREF(&AlignmentType);
//  Py_INCREF(&MapType);
//  Py_INCREF(&MismatchType);
//  Py_INCREF(&gempy_alignment_iteratorType);
//  Py_INCREF(&gempy_alignment_mappings_iteratorType);

//  PyModule_AddObject(m, "_template_iterator", (PyObject *)&gempy_template_iteratorType);
//  PyModule_AddObject(m, "_iterator", (PyObject *)&gempy_iteratorType);
  PyModule_AddObject(m, "Template", (PyObject *)&TemplateType);
//  PyModule_AddObject(m, "Alignment", (PyObject *)&AlignmentType);
//  PyModule_AddObject(m, "Map", (PyObject *)&MapType);
//  PyModule_AddObject(m, "Mismatch", (PyObject *)&MismatchType);
//  PyModule_AddObject(m, "_mappings_iterator", (PyObject *)&gempy_mappings_iteratorType);
//  PyModule_AddObject(m, "_alignment_iterator", (PyObject *)&gempy_alignment_iteratorType);
//  PyModule_AddObject(m, "_alignment_mappings_iterator", (PyObject *)&gempy_alignment_mappings_iteratorType);
//
//
//  // define some constants in the gemtools module
//  PyModule_AddObject(m, "FORWARD", PyInt_FromLong(FORWARD));
//  PyModule_AddObject(m, "REVERSE", PyInt_FromLong(REVERSE));
//  PyModule_AddObject(m, "SKIP_NO_JUNCTION", PyInt_FromLong(NO_JUNCTION));
//  PyModule_AddObject(m, "SKIP_SPLICE", PyInt_FromLong(SPLICE));
//  PyModule_AddObject(m, "SKIP_POSITIVE", PyInt_FromLong(POSITIVE_SKIP));
//  PyModule_AddObject(m, "SKIP_NEGATIVE", PyInt_FromLong(NEGATIVE_SKIP));
//  PyModule_AddObject(m, "SKIP_INSERT", PyInt_FromLong(INSERT));
//
//  PyModule_AddObject(m, "MISMS", PyInt_FromLong(MISMS));
//  PyModule_AddObject(m, "INS", PyInt_FromLong(INS));
//  PyModule_AddObject(m, "DEL", PyInt_FromLong(DEL));

}
