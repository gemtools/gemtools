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

void Alignment_dealloc(PyObject* self){
    self->ob_type->tp_free(self);
}


PyObject* Alignment_gettag(Alignment *self, void *closure){
    char* tag = gempy_alignment_get_tag(self);
    if(tag){
        PyObject* ret = PyString_FromString(tag);
        // GT-13 workaround we need to free the tag if we
        // appended /1 or /2
        if(self->index){
            free(tag);
        }
        return ret;
    }
    Py_RETURN_NONE;
}

char* gempy_alignment_get_tag(Alignment* self){
    char* tag = gt_alignment_get_tag(self->alignment);
    if(tag){
        return tag;
    }else if(self->template){
        char* org = gt_template_get_tag(self->template);
        // GT-13 workaround
        if(self->index){
            char* new = malloc((strlen(org)+2)* sizeof(char));
            sprintf(new, "%s/%"PRIx64, org, self->index);
            return new;
        }
        return org;
    }
    return NULL;
}

int Alignment_settag(Alignment *self, PyObject *value, void *closure){
    //gt_alignment_set_tag(self->alignment, PyString_AsString(value));
    return 0;
}

PyObject* Alignment_getread(Alignment *self, void *closure){
    char* read = gt_alignment_get_read(self->alignment);
    if(read){
        PyObject* ret = PyString_FromString(read);
        //Py_DECREF(ret);
        return ret;
    }else{
        Py_RETURN_NONE;
    }
}

int Alignment_setread(Alignment *self, PyObject *value, void *closure){
    //gt_alignment_set_read(self->alignment, PyString_AsString(value));
    return 0;
}

PyObject* Alignment_getqualities(Alignment *self, void *closure){
    char* qual = gt_alignment_get_qualities(self->alignment);
    if(qual){
        PyObject* ret = PyString_FromString(qual);
        //Py_DECREF(ret);
        return ret;
    }else{
        Py_RETURN_NONE;
    }
}

int Alignment_setqualities(Alignment *self, PyObject *value, void *closure){
    //gt_alignment_set_qualities(self->alignment, PyString_AsString(value));
    return 0;
}

PyObject* Alignment_getmax_complete_strata(Alignment *self, void *closure){
    uint64_t mcs = gt_alignment_get_mcs(self->alignment);
    PyObject* ret;
    if( mcs != UINT64_MAX){
        ret = PyLong_FromUnsignedLongLong(mcs);
    }else{
        ret = PyLong_FromLongLong(-1);
    }
    //Py_DECREF(ret);
    return ret;
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
    //gt_alignment_set_mcs(self->alignment, PyLong_AsUnsignedLongLong(value));
    return 0;
}

PyObject* Alignment_get_counters(PyObject *self, PyObject *closure){
    Alignment* ali = (Alignment*) self;
    PyObject* ret = create_gempy_iterator(1, gt_alignment_get_num_counters(ali->alignment), gt_alignment_get_counter, ali->alignment, PyLong_FromUnsignedLongLong, 0);
    return ret;
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
    Py_DECREF(result);
    free(fastq);  
    // GT-13 workaround we need to free the tag if we
    // appended /1 or /2
    if(a->index){
       free(tag);
    }
    return result;
}
