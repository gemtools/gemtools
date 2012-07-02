#include "py_map.h"
#include "py_mismatch.h"

static int Map_init(Map *self, PyObject *args, PyObject *kwds){
    self->map = gt_map_new();
    return 0;
}

static PyObject* Map_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
    Map *self;
    self = (Map *)type->tp_alloc(type, 0);
    Map_init(self, args, kwds);
    return (PyObject *)self;
}

static PyObject* Map_getseq_name(Map *self, void *closure){
    char* tag = gt_map_get_seq_name(self->map);
    if(tag){
        return PyString_FromString(tag);
    }else{
        Py_RETURN_NONE;
    }
}

static int Map_setseq_name(Map *self, PyObject *value, void *closure){
    gt_map_set_seq_name(self->map, PyString_AsString(value));
    return 0;
}

static PyObject* Map_getposition(Map *self, void *closure){
    return PyLong_FromUnsignedLongLong(gt_map_get_position(self->map));
}

static int Map_setposition(Map *self, PyObject *value, void *closure){
    gt_map_set_position(self->map, PyLong_AsUnsignedLongLong(value));
    return 0;
}

static PyObject* Map_getscore(Map *self, void *closure){
    return PyLong_FromUnsignedLongLong(gt_map_get_score(self->map));
}

static int Map_setscore(Map *self, PyObject *value, void *closure){
    gt_map_set_score(self->map, PyLong_AsUnsignedLongLong(value));
    return 0;
}

static PyObject* Map_getbase_lengt(Map *self, void *closure){
    return PyLong_FromUnsignedLongLong(gt_map_get_base_length(self->map));
}

static int Map_setbase_length(Map *self, PyObject *value, void *closure){
    gt_map_set_base_length(self->map, PyLong_AsUnsignedLongLong(value));
    return 0;
}

static PyObject* Map_getdirection(Map *self, void *closure){
    return PyInt_FromLong(gt_map_get_direction(self->map));
}

static int Map_setdirection(Map *self, PyObject *value, void *closure){
    gt_map_set_direction(self->map, PyInt_AsLong(value));
    return 0;
}

static PyObject* Map_getmismatches(Map *self, void *closure)
{
    return py_iterator(0, gt_map_get_num_misms(self->map), gt_map_get_misms, self->map, create_mismatch, 0);
}

static int Map_setmismatches(Map *self, PyObject *value, void *closure){
    PyErr_SetString(PyExc_TypeError, "Setting mismatches is currently not supported");
    return -1;
}

Map* create_map(gt_map* map){
    Map* a = PyObject_New(Map, &MapType);
    a->map = map;
    return a;
}
