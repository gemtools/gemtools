#include "Python.h"
#include "gem_tools.h"


#include "py_mismatch.h"
#include "py_map.h"
#include "py_alignment.h"
#include "py_template.h"
#include "py_template_iterator.h"
#include "py_iterator.h"


static PyObject* gempy_open_stream(PyObject *self, PyObject *args){
    void* m;
    gempy_template_iterator *p;
    if (!PyArg_ParseTuple(args, "O", &m))  return NULL;
    p = create_template_stream_iterator(PyFile_AsFile(m));
    if (!p) return NULL;
    return (PyObject *)p;
}

static PyObject* gempy_open_file(PyObject *self, PyObject *args, PyObject *keywds){
    char* filename;
    gempy_template_iterator *p;
    if (!PyArg_ParseTuple(args, "s", &m))  return NULL;

    p = create_template_file_iterator(filename, false);
    if (!p) return NULL;
    return (PyObject *)p;
}

static PyMethodDef GempyMethods[] = {
    {"open_file", gempy_open_file, METH_VARARGS, "Iterator over a .map file"},
    {"open_stream", gempy_open_stream, METH_VARARGS, "Iterator over a .map stream"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC initgempy(void){
  PyObject* m;

  gempy_iterator.tp_new = PyType_GenericNew;
  gempy_template_iterator.tp_new = PyType_GenericNew;
  if (PyType_Ready(&gempy_iteratorType) < 0)  return;
  if (PyType_Ready(&gempy_template_iteratorType) < 0)  return;
  if (PyType_Ready(&MismatchType) < 0)  return;
  if (PyType_Ready(&MapType) < 0)  return;
  if (PyType_Ready(&AlignmentType) < 0)  return;
  if (PyType_Ready(&TemplateType) < 0)  return;

  m = Py_InitModule("gempy", GempyMethods);

  Py_INCREF(&gempy_iteratorType);
  Py_INCREF(&gempy_template_iteratorType);
  Py_INCREF(&TemplateType);
  Py_INCREF(&AlignmentType);
  Py_INCREF(&MapType);
  Py_INCREF(&MismatchType);

  PyModule_AddObject(m, "_template_iterator", (PyObject *)&gempy_template_iteratorType);
  PyModule_AddObject(m, "_iterator", (PyObject *)&gempy_iteratorType);
  PyModule_AddObject(m, "Template", (PyObject *)&TemplateType);
  PyModule_AddObject(m, "Alignment", (PyObject *)&AlignmentType);
  PyModule_AddObject(m, "Map", (PyObject *)&MapType);
  PyModule_AddObject(m, "Mismatch", (PyObject *)&MismatchType);
}
