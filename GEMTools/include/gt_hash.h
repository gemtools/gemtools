/*
 * PROJECT: GEM-Tools library
 * FILE: gt_hash.c
 * DATE: 2/09/2012
 * DESCRIPTION: // TODO
 */

#ifndef GT_HASH_H_
#define GT_HASH_H_

#include "gt_ihash.h"
#include "gt_shash.h"

// TODO: Insert string + duplicate it properly on attr copy

/*
 * Checkers
 */
#define GT_HASH_CHECK(hash) gt_fatal_check(hash==NULL,NULL_HANDLER)

#endif /* GT_HASH_H_ */
