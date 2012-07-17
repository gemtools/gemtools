#include "py_mappings_iterator.h"

/*** Iterata over the map blocks of a tempalte ***/
PyObject* gempy_alignment_iterator_iter(PyObject *self){
  Py_INCREF(self);
  return self;
}

PyObject* gempy_alignment_iterator_iternext(PyObject *self){
    gempy_alignment_iterator* p = (gempy_alignment_iterator*) self;
    gt_template* template = p->template;
    gempy_mappings_iterator* iterator = NULL;
    if(p->end_position < p->num_blocks){
        iterator = create_mappings_iterator(p->map_array[p->end_position++]);
        iterator->alignment_iterator = p;
        //Py_INCREF(iterator);
        //Py_DECREF(iterator);
        return (PyObject*) iterator;
    }else{
        p->end_position = 1;
        if(gt_template_next_maps(&(p->maps_iterator),&(p->map_array))){
            iterator = create_mappings_iterator(p->map_array[0]);
            iterator->alignment_iterator = p;
            //Py_INCREF(iterator);
            //Py_DECREF(iterator);
            return (PyObject*) iterator;
        }
    }
    /* Raising of standard StopIteration exception with empty value. */
    PyErr_SetNone(PyExc_StopIteration);
    return (PyObject*) NULL;
}

void gempy_alignment_iterator_dealloc(gempy_alignment_iterator* self){
    ////printf("Dealloc alignment_iterator %p\n", self);
    self->ob_type->tp_free((PyObject*)self);
}



/*** Iterata over the map blocks of an alignment ***/
PyObject* gempy_alignment_mappings_iterator_iter(PyObject *self){
  Py_INCREF(self);
  return self;
}

PyObject* gempy_alignment_mappings_iterator_iternext(PyObject *self){
    gempy_alignment_mappings_iterator* p = (gempy_alignment_mappings_iterator*) self;
    //printf("Checking for next mappings iterator... %p\n", self);
    if(p->current < p->total){
        //printf("Creating for next mappings iterator... %p\n", self);
        gt_map* next_map = gt_alignment_get_map(p->alignment, p->current++);
        PyObject* ret =  (PyObject*) create_mappings_iterator(next_map);
        //Py_DECREF(ret);
        return ret;
    }
    //printf("Stopping mappings iterator...\n, %p", self);
    /* Raising of standard StopIteration exception with empty value. */
    PyErr_SetNone(PyExc_StopIteration);
    return (PyObject*) NULL;
}

void gempy_alignment_mappings_iterator_dealloc(gempy_alignment_iterator* self){
    //printf("Dealloc alignment_mappings_iterator %p\n", self);
    self->ob_type->tp_free((PyObject*)self);
}



/**** Iterate over map_blocks ***/
PyObject* gempy_mappings_iterator_iter(PyObject *self){
    Py_INCREF(self);
    return self;
}

PyObject* gempy_mappings_iterator_iternext(PyObject *self){
    gempy_mappings_iterator* p = (gempy_mappings_iterator*) self;
    //printf("Iterator next call...%p %d %d?\n", p, p, p->map_block);
    if(p->map_block != NULL){
        //printf("got block?\n");
        gt_map* curr = p->map_block;
        PyObject* ret = (PyObject*) create_map(p->map_block);
        int64_t junction_distance = gt_map_get_next_block_distance(p->map_block);
        gt_junction_t junction = gt_map_get_next_block_junction(p->map_block);
        p->map_block = gt_map_get_next_block(p->map_block);
        //p->map_block = p->map_block->next_block == NULL ? NULL : p->map_block->next_block->map;
        //printf("Returning block? %d\n", p->map_block);
        // jump to the next block ?
        if(p->map_block == NULL && p->alignment_iterator != NULL){
            //printf("Update block from alignment iterator?\n");
            gempy_alignment_iterator* ali_iter = p->alignment_iterator;
            if(ali_iter->end_position < ali_iter->num_blocks){
                // there is a next block
                curr = ali_iter->map_array[ali_iter->end_position-1];
                p->map_block = ali_iter->map_array[ali_iter->end_position++];
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
    //printf("Stopping iteration\n");
    /* Raising of standard StopIteration exception with empty value. */
    PyErr_SetNone(PyExc_StopIteration);
    return (PyObject*) NULL;
}

void gempy_mappings_iterator_dealloc(gempy_mappings_iterator* self){
    //printf("Dealloc mappings_iterator %p\n", self);
    self->ob_type->tp_free((PyObject*)self);
}

