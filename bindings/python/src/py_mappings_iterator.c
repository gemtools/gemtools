#include "py_mappings_iterator.h"

/*** Iterata over the map blocks of a tempalte ***/
PyObject* gempy_alignment_iterator_iter(PyObject *self){
  Py_INCREF(self);
  return self;
}

PyObject* gempy_alignment_iterator_iternext(PyObject *self){
    gempy_alignment_iterator* p = (gempy_alignment_iterator*) self;
    gempy_mappings_iterator* iterator = NULL;
    p->position += p->jump;
    if(p->jump < 0 || p->position >= p->num_blocks){
        // initialize
        p->jump = gt_template_get_num_blocks(p->template);
        p->position = 0;
        if(gt_template_next_maps(&(p->maps_iterator),&(p->map_array))){
            iterator = create_mappings_iterator(p->map_array[0]);
            iterator->alignment_iterator = p;
            gt_mmap_attributes* attrs = gt_template_get_mmap_attr(p->template, 0);
            iterator->distance = attrs->distance;
            iterator->score = attrs->score;
            return (PyObject*) iterator;
        }
    }else{
        // check next position
        p->position += p->jump;
        iterator = create_mappings_iterator(p->map_array[p->position]);
        iterator->alignment_iterator = p;
        // get the attributes
        gt_mmap_attributes* attrs = gt_template_get_mmap_attr(p->template, p->position);
        iterator->distance = attrs->distance;
        iterator->score = attrs->score;
        // reset jump size
        p->jump = gt_template_get_num_blocks(p->template);
        return (PyObject*) iterator;
    }
    /* Raising of standard StopIteration exception with empty value. */
    PyErr_SetNone(PyExc_StopIteration);
    return (PyObject*) NULL;
}

void gempy_alignment_iterator_dealloc(PyObject* self){
    self->ob_type->tp_free(self);
}



/*** Iterata over the map blocks of an alignment ***/
PyObject* gempy_alignment_mappings_iterator_iter(PyObject *self){
  Py_INCREF(self);
  return self;
}

PyObject* gempy_alignment_mappings_iterator_iternext(PyObject *self){
    gempy_alignment_mappings_iterator* p = (gempy_alignment_mappings_iterator*) self;
    if(p->current < p->total){
        gt_map* next_map = gt_alignment_get_map(p->alignment, p->current++);
        gempy_mappings_iterator* ret =  create_mappings_iterator(next_map);
        gt_mmap_attributes* attrs = gt_template_get_mmap_attr(p->template, p->current-1);
        ret->distance = -1;
        ret->score = -1;
        if (attrs){
            ret->distance = attrs->distance;
            ret->score = attrs->score;
        }
        //Py_DECREF(ret);
        return (PyObject*) ret;
    }
    /* Raising of standard StopIteration exception with empty value. */
    PyErr_SetNone(PyExc_StopIteration);
    return (PyObject*) NULL;
}

void gempy_alignment_mappings_iterator_dealloc(PyObject* self){
    self->ob_type->tp_free(self);
}



/**** Iterate over map_blocks ***/
PyObject* gempy_mappings_iterator_iter(PyObject *self){
    Py_INCREF(self);
    return self;
}

PyObject* gempy_mappings_iterator_iternext(PyObject *self){
    gempy_mappings_iterator* p = (gempy_mappings_iterator*) self;
    if(p->map_block != NULL){
        gt_map* curr = p->map_block;
        PyObject* ret = (PyObject*) create_map(p->map_block);
        int64_t junction_distance = gt_map_get_next_block_distance(p->map_block);
        gt_junction_t junction = gt_map_get_next_block_junction(p->map_block);
        p->map_block = gt_map_get_next_block(p->map_block);
        // jump to the next block ?
        if(p->map_block == NULL && p->alignment_iterator != NULL){
            gempy_alignment_iterator* ali_iter = p->alignment_iterator;
            uint64_t next_pos = ali_iter->position + 1;
            if(next_pos < ali_iter->num_blocks){
                // there is a next block
                curr = ali_iter->map_array[ali_iter->position];
                p->map_block = ali_iter->map_array[next_pos];
                ali_iter->jump--; // reduce jump size
                ali_iter->position = next_pos;
                if(p->map_block != NULL){
                    junction = INSERT;
                    if(gt_map_get_direction(curr) == REVERSE){
                        junction_distance = gt_map_get_position(curr) - gt_map_get_position(p->map_block);
                    }else{
                        junction_distance = gt_map_get_position(p->map_block) - gt_map_get_position(curr);
                    }
                }
            }
        }

        // return tuple of Map, junction, distance
        //return ret;
        PyObject* junc = PyInt_FromLong(junction);
        PyObject* dist = PyInt_FromLong(junction_distance);
        PyObject* tup = PyTuple_Pack(3, ret, junc, dist);
        Py_DECREF(ret);
        Py_DECREF(junc);
        Py_DECREF(dist);
        //Py_DECREF(tup);
        return tup;
    }
    /* Raising of standard StopIteration exception with empty value. */
    PyErr_SetNone(PyExc_StopIteration);
    return (PyObject*) NULL;
}

void gempy_mappings_iterator_dealloc(PyObject* self){
    self->ob_type->tp_free(self);
}

