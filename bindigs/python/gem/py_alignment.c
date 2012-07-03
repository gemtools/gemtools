#include "py_alignment.h"
#include "py_iterator.h"

int Alignment_init(Alignment *self, PyObject *args, PyObject *kwds){
    self->alignment = gt_alignment_new();
    self->template = NULL;
    return 0;
}

PyObject* Alignment_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
    Alignment *self;
    self = (Alignment *)type->tp_alloc(type, 0);
    Alignment_init(self, args, kwds);
    return (PyObject *)self;
}

PyObject* Alignment_gettag(Alignment *self, void *closure){
    char* tag = gempy_alignment_get_tag(self);
    if(tag){
        return PyString_FromString(tag);
    }else if(self->template){
        Py_RETURN_NONE;
    }
}

char* gempy_alignment_get_tag(Alignment* self){
    char* tag = gt_alignment_get_tag(self->alignment);
    if(tag){
        return tag;
    }else if(self->template){
        return gt_template_get_tag(self->template);
    }
}

int Alignment_settag(Alignment *self, PyObject *value, void *closure){
    gt_alignment_set_tag(self->alignment, PyString_AsString(value));
    return 0;
}

PyObject* Alignment_getread(Alignment *self, void *closure){
    char* read = gt_alignment_get_read(self->alignment);
    if(read){
        return PyString_FromString(read);
    }else{
        Py_RETURN_NONE;
    }
}

int Alignment_setread(Alignment *self, PyObject *value, void *closure){
    gt_alignment_set_read(self->alignment, PyString_AsString(value));
    return 0;
}

PyObject* Alignment_getqualities(Alignment *self, void *closure){
    char* qual = gt_alignment_get_qualities(self->alignment);
    if(qual){
        return PyString_FromString(qual);
    }else{
        Py_RETURN_NONE;
    }
}

int Alignment_setqualities(Alignment *self, PyObject *value, void *closure){
    gt_alignment_set_qualities(self->alignment, PyString_AsString(value));
    return 0;
}

PyObject* Alignment_getmax_complete_strata(Alignment *self, void *closure){
    return PyLong_FromUnsignedLongLong(gt_alignment_get_mcs(self->alignment));
}

int Alignment_setmax_complete_strata(Alignment *self, PyObject *value, void *closure){
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the mcs attribute");
        return -1;
    }
    if (! PyLong_Check(value)) {
        PyErr_SetString(PyExc_TypeError,
                        "The first attribute value must be a long");
                        return -1;
    }
    gt_alignment_set_mcs(self->alignment, PyLong_AsUnsignedLongLong(value));
    return 0;
}

PyObject* Alignment_getcounters(Alignment *self, void *closure){
    return create_gempy_iterator(1, gt_alignment_get_num_counters(self->alignment), gt_alignment_get_counter, self->alignment, PyLong_FromUnsignedLongLong, 0);
}

int Alignment_setcounters(Alignment *self, PyObject *value, void *closure){
    PyErr_SetString(PyExc_TypeError, "Setting blocks is currently not supported");
    return -1;
}

PyObject* Alignment_to_sequence(PyObject *self, PyObject *args){

    Alignment* a = (Alignment*) self;
    char* tag = gempy_alignment_get_tag(a);
    char* read = gt_alignment_get_read(a->alignment);
    char* qualities = gt_alignment_get_qualities(a->alignment);
    char* fastq = NULL;
    PyObject* result = NULL;

    if(qualities){
        fastq = malloc((strlen(tag)+strlen(read)+strlen(qualities)+5)* sizeof(char));
        sprintf(fastq, "@%s\n%s\n+\n%s", tag, read, qualities);
    }else{
        fastq = malloc((strlen(tag)+strlen(read)+2)* sizeof(char));
        sprintf(fastq, ">%s\n%s", tag, read);
    }
    result = PyString_FromString(fastq);
    free(fastq);
    return result;
}


