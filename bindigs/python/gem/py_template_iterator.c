#include "py_template_iterator.h"
#include "py_template.h"

PyObject* gempy_template_iterator_iter(PyObject *self){
  Py_INCREF(self);
  return self;
}

PyObject* gempy_template_iterator_iternext(PyObject *self){
    gempy_template_iterator *p = (gempy_template_iterator *)self;
    gt_status error_code;
    Template* tmpl = NULL;
    gt_template* template = p->template;
    while ((error_code=gt_buffered_map_input_get_template(p->map_input,p->template))) {
        if (error_code==GT_BMI_FAIL) continue;
        tmpl = create_template(template);
        return (PyObject*) tmpl;
    }
    /* Raising of standard StopIteration exception with empty value. */
    PyErr_SetNone(PyExc_StopIteration);
    return (PyObject*) NULL;
}

void gempy_template_iterator_dealloc(gempy_template_iterator* self){
    // close the stream and map
    gt_buffered_map_input_close(self->map_input);
    // Close files
    gt_input_file_close(self->input_file);
    self->ob_type->tp_free((PyObject*)self);
}

