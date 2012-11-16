#include "py_map.h"

#include "py_mismatch.h"
#include "py_iterator.h"


//int Map_init(Map *self, PyObject *args, PyObject *kwds){
//    self->map = gt_map_new();
//    return 0;
//}
//
//PyObject* Map_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
//    Map *self;
//    self = (Map *)type->tp_alloc(type, 0);
//    Map_init(self, args, kwds);
//    return (PyObject *)self;
//}
//
//void Map_dealloc(PyObject* self){
//    self->ob_type->tp_free(self);
//}
//
//
//PyObject* Map_getseq_name(Map *self, void *closure){
//    char* tag = gt_map_get_seq_name(self->map);
//    if(tag){
//        PyObject* ret = PyString_FromString(tag);
//        return ret;
//    }else{
//        Py_RETURN_NONE;
//    }
//}
//
//int Map_setseq_name(Map *self, PyObject *value, void *closure){
//    //gt_map_set_seq_name(self->map, PyString_AsString(value));
//    return 0;
//}
//
//PyObject* Map_getposition(Map *self, void *closure){
//    PyObject* ret =  PyLong_FromUnsignedLongLong(gt_map_get_position(self->map));
//    return ret;
//}
//
//int Map_setposition(Map *self, PyObject *value, void *closure){
//    gt_map_set_position(self->map, PyLong_AsUnsignedLongLong(value));
//    return 0;
//}
//
//PyObject* Map_getscore(Map *self, void *closure){
//    PyObject* ret = PyLong_FromUnsignedLongLong(gt_map_get_score(self->map));
//    return ret;
//}
//
//int Map_setscore(Map *self, PyObject *value, void *closure){
//    gt_map_set_score(self->map, PyLong_AsUnsignedLongLong(value));
//    return 0;
//}
//
//// distance
//PyObject* Map_getdistance(Map *self, void *closure){
//    PyObject* ret = PyLong_FromUnsignedLongLong(gt_map_get_distance(self->map));
//    return ret;
//}
//
//PyObject* Map_getbase_length(Map *self, void *closure){
//    PyObject* ret = PyLong_FromUnsignedLongLong(gt_map_get_base_length(self->map));
//    return ret;
//}
//
//int Map_setbase_length(Map *self, PyObject *value, void *closure){
//    gt_map_set_base_length(self->map, PyLong_AsUnsignedLongLong(value));
//    return 0;
//}
//
//PyObject* Map_getdirection(Map *self, void *closure){
//    //PyObject* ret = PyInt_FromLong(gt_map_get_direction(self->map));
//    PyObject* ret = PyInt_FromLong(0);
//    return ret;
//}
//
//int Map_setdirection(Map *self, PyObject *value, void *closure){
//    //gt_map_set_direction(self->map, PyInt_AsLong(value));
//    return 0;
//}
//
//
//// length with indels
//PyObject* Map_getlength(Map *self, void *closure){
//    PyObject* ret =  PyLong_FromUnsignedLongLong(gt_map_get_length(self->map));
//    return ret;
//}
//
//
//// levenshtein
//PyObject* Map_getlevenshtein(Map *self, void *closure){
//    PyObject* ret = PyLong_FromUnsignedLongLong(gt_map_get_levenshtein_distance(self->map));
//    return ret;
//}
//
//
//// global length with indels ?
//PyObject* Map_getglobal_length(Map *self, void *closure){
//    PyObject* ret = PyLong_FromUnsignedLongLong(gt_map_get_global_length(self->map));
//    return ret;
//}
//
//// global distance
//PyObject* Map_getglobal_distance(Map *self, void *closure){
//    PyObject* ret = PyLong_FromUnsignedLongLong(gt_map_get_global_distance(self->map));
//    return ret;
//}
//
//// global score
//PyObject* Map_getglobal_score(Map *self, void *closure){
//    PyObject* ret =  PyLong_FromUnsignedLongLong(gt_map_get_global_score(self->map));
//    return ret;
//}
//
//// global levenshtein
//PyObject* Map_getglobal_levenshtein(Map *self, void *closure){
//    PyObject* ret = PyLong_FromUnsignedLongLong(gt_map_get_global_levenshtein_distance(self->map));
//    return ret;
//}
//
//;
//// number of mismatches
//PyObject* Map_getnum_mismatches(Map *self, void *closure){
//    PyObject* ret = PyLong_FromUnsignedLongLong(gt_map_get_num_misms(self->map));
//    return ret;
//}
//
//PyObject* Map_get_mismatches(PyObject *self, PyObject *closure)
//{
//    Map* mm = (Map*) self;
//    PyObject* ret = create_gempy_iterator(0, gt_map_get_num_misms(mm->map), gt_map_get_misms, mm->map, create_mismatch, 0);
//    return ret;
//}
