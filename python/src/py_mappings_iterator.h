#ifndef GEMPY_PY_MAPPINGS_ITERATOR
#define GEMPY_PY_MAPPINGS_ITERATOR

//#include "Python.h"
//#include "gem_tools.h"
//#include "py_map.h"
//
//typedef struct {
//  PyObject_HEAD
//  gt_template* template;
//  gt_map** map_array;
//  gt_template_maps_iterator maps_iterator;
//  uint64_t position;
//  int64_t jump;
//  uint64_t num_blocks;
//} gempy_alignment_iterator;
//
//typedef struct {
//  PyObject_HEAD
//  gt_alignment* alignment;
//  gt_template* template;
//  uint64_t current;
//  uint64_t total;
//} gempy_alignment_mappings_iterator;
//
//typedef struct {
//  PyObject_HEAD
//  gt_map* map_block;
//  gempy_alignment_iterator* alignment_iterator;
//  int64_t distance;
//  int64_t score;
//} gempy_mappings_iterator;
//
//PyObject* gempy_alignment_iterator_iter(PyObject *self);
//PyObject* gempy_alignment_iterator_iternext(PyObject *self);
//void gempy_alignment_iterator_dealloc(PyObject* self);
//
//PyObject* gempy_alignment_mappings_iterator_iter(PyObject *self);
//PyObject* gempy_alignment_mappings_iterator_iternext(PyObject *self);
//void gempy_alignment_mappings_iterator_dealloc(PyObject* self);
//
//
//PyObject* gempy_mappings_iterator_iter(PyObject *self);
//PyObject* gempy_mappings_iterator_iternext(PyObject *self);
//void gempy_mappings_iterator_dealloc(PyObject* self);
//
//gempy_alignment_mappings_iterator* create_alignment_mappings_iterator(gt_alignment* alignment);
//gempy_alignment_iterator* create_template_mappings_iterator(gt_template* template);
//gempy_mappings_iterator* create_mappings_iterator(gt_map* map_block);
//
#endif
