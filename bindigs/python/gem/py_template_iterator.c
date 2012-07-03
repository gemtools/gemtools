#include "py_template_iterator.h"
#include "py_template.h"

PyObject* gempy_template_iterator_iter(PyObject *self){
  Py_INCREF(self);
  return self;
}

PyObject* gempy_template_iterator_iternext(PyObject *self){
/*    gempy_template_iterator *p = (gempy_template_iterator *)self;
    gt_status error_code;
    Template* tmpl = NULL;
    gt_template* template = p->template;
    while ((error_code=gt_buffered_map_input_get_template(p->map_input,p->template))) {
        if (error_code==GT_BMI_FAIL) continue;
        GT_TEMPLATE_ITERATE(template, map_array) {
            tmpl = create_template(template);
            return (PyObject*) tmpl;
        }
    }
*/
    printf("Nothing to iterate on \n");
    /* Raising of standard StopIteration exception with empty value. */
    PyErr_SetNone(PyExc_StopIteration);
    return (PyObject*) NULL;
}

gempy_template_iterator* create_template_stream_iterator(FILE* file){
    gempy_template_iterator *p;
    p = PyObject_New(gempy_template_iterator, &gempy_template_iteratorType);
    if (!p) return NULL;
    p->map_input = gt_buffered_map_input_new(gt_input_stream_open(file));
    p->template = gt_template_new();
    return (gempy_template_iterator *)p;
}

gempy_template_iterator* create_template_file_iterator(char* filename, bool memorymap){
    gempy_template_iterator *p;
    p = PyObject_New(gempy_template_iterator, &gempy_template_iteratorType);
    if (!p) return NULL;
    /* I'm not sure if it's strictly necessary. */
    if (!PyObject_Init((PyObject *)p, &gempy_template_iteratorType)) {
        printf("ERROR template iterator init failed!\n");
        Py_DECREF(p);
        return NULL;
    }
    p->map_input = gt_buffered_map_input_new(gt_input_file_open(filename, memorymap)); // false disable memory map
    p->template = gt_template_new();
    return (gempy_template_iterator *)p;
}
