/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2010 Carnegie Mellon University.  All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced
 * Research Projects Agency and the National Science Foundation of the
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */

/**
 * @file vocab_map.h
 * @brief Vocabulary mapping
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#ifndef __VOCAB_MAP_H__
#define __VOCAB_MAP_H__

#include <stdio.h>

#include "dict.h"

/**
 * Vocabulary mapping object.
 */
typedef struct vocab_map_s vocab_map_t;

/**
 * Create a new vocabulary map.
 */
vocab_map_t *vocab_map_init(dict_t *dict);

/**
 * Retain a pointer to a vocabulary map.
 */
vocab_map_t *vocab_map_retain(vocab_map_t *vm);

/**
 * Release a pointer to a vocabulary map.
 */
int vocab_map_free(vocab_map_t *vm);

/**
 * Get the dictionary from the vocabulary map.
 */
dict_t *vocab_map_dict(vocab_map_t *vm);

/**
 * Read a vocabulary map from a file.
 */
int vocab_map_read(vocab_map_t *vm, FILE *fh);

/**
 * Write a vocabulary map to a file.
 */
int vocab_map_write(vocab_map_t *vm, FILE *fh);

/**
 * Map a word to a corresponding pseudo word, if any.
 */
int32 vocab_map_map(vocab_map_t *vm, int32 wid);

/**
 * Map a pseudo word to its constituents.
 */
int32 const *vocab_map_unmap(vocab_map_t *vm, int32 pseudo_wid,
                             int32 *out_n_mapped);

/**
 * Iterator over vocabulary mappings.
 */
typedef struct vocab_map_iter_s vocab_map_iter_t;

/**
 * Iterate over all mappings in a vocabulary map.
 */
vocab_map_iter_t *vocab_map_mappings(vocab_map_t *vm);

/**
 * Move the iterator forward.
 */
vocab_map_iter_t *vocab_map_iter_next(vocab_map_iter_t *itor);

/**
 * Free an iterator early.
 */
void vocab_map_iter_free(vocab_map_iter_t *itor);

/**
 * Get pseudoword and constituents from a vocab map iterator.
 */
int32 const *vocab_map_iter_get(vocab_map_iter_t *itor,
                                int32 *out_pseudo_wid,
                                int32 *out_n_mapped);


#endif /* __VOCAB_MAP_H__ */
