/*
 * PROJECT: GEM-Tools library
 * FILE: gt_commons.c
 * DATE: 01/06/2012
 * DESCRIPTION: // TODO
 */

#include "gt_commons.h"

bool gt_dna[256] =
{
    [0 ... 255] = false,
    ['A'] = true,['C'] = true,['G'] = true,['T'] = true,
    ['a'] = true,['c'] = true,['g'] = true,['t'] = true,
    ['N'] = true
};
