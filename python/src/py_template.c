#include "py_template.h"
#include "py_iterator.h"
#include "py_alignment.h"
#include "gem_tools.h"

int Template_init(Template *self, PyObject *args, PyObject *kwds){
    self->template = gt_template_new();
    return 0;
}

void Template_dealloc(PyObject* self){
    self->ob_type->tp_free(self);
}


PyObject* Template_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
    Template *self = (Template *)type->tp_alloc(type, 0);
    Template_init(self, args, kwds);
    return (PyObject *)self;
}

PyObject* Template_fill(Template *self, PyObject *args){
    PyObject *line;
    if (PyArg_UnpackTuple(args, "ref", 1, 1, &line)) {
        if(gt_input_map_parse_template(PyString_AsString(line), self->template) != 0){
            PyErr_SetString(PyExc_Exception, "Error while parsing mapping!");
            return NULL;
        }
    }
    Py_RETURN_TRUE;
}
//
//PyObject* Template_gettag(Template *self, void *closure){
//    PyObject* ret = PyString_FromString(gt_template_get_tag(self->template));
//    //Py_DECREF(ret);
//    return ret;
//}
//
//int Template_settag(Template *self, PyObject *value, void *closure){
//    //gt_template_set_tag(self->template, PyString_AsString(value));
//    return 0;
//}
//

PyObject* Template_get_max_complete_strata(Template *self, void *closure){
    PyObject* ret = PyLong_FromUnsignedLongLong(gt_template_get_mcs(self->template));
    //Py_DECREF(ret);
    return ret;
}
//
//PyObject* Template_get_blocks(PyObject *self, PyObject* closure){
//    Template* tmpl = (Template*) self;
//    PyObject* ret = create_gempy_iterator(0, gt_template_get_num_blocks(tmpl->template), gt_template_get_block, tmpl->template, create_alignment, 1);
//    return ret;
//}
//
//PyObject* Template_get_counters(PyObject* self, PyObject *closure){
//    Template* tmpl = (Template*) self;
//    PyObject* ret = create_gempy_iterator(1, gt_template_get_num_counters(tmpl->template), gt_template_get_counter, tmpl->template, PyLong_FromUnsignedLongLong, 0);
//    return ret;
//}
//
//PyObject* Template_get_num_blocks(Template *self, void *closure){
//    PyObject* ret = PyLong_FromUnsignedLongLong(gt_template_get_num_blocks(self->template));
//    return ret;
//}
//PyObject* Template_get_num_counters(Template *self, void *closure){
//    PyObject* ret = PyLong_FromUnsignedLongLong(gt_template_get_num_counters(self->template));
//    return ret;
//}

