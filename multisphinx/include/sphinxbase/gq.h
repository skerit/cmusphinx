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
 * \file gq.h
 * \brief Generic double-ended queue
 */

#ifndef __SPHINXBASE_GQ_H__
#define __SPHINXBASE_GQ_H__

/* Win32/WinCE DLL gunk */
#include <sphinxbase/sphinxbase_export.h>
#include <sphinxbase/prim_type.h>

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif
#if 0
/* Fool Emacs. */
}
#endif

/**
 * Generic double-ended queue.
 */
typedef struct gq_s gq_t;

/**
 * Create a new, empty queue.
 */
gq_t *gq_init(size_t ent_size);

/**
 * Retain a pointer to a gq_t.
 */
gq_t *gq_retain(gq_t *q);

/**
 * Release a pointer to a gq_t.
 */
int gq_free(gq_t *q);

/**
 * Get the number of elements in a gq_t.
 */
size_t gq_size(gq_t *q);

/**
 * Get the size of an element in a gq_t.
 */
size_t gq_ent_size(gq_t *q);

/**
 * Get the number of elements allocated in a gq_t.
 */
size_t gq_alloc_size(gq_t *q);


/**
 * Append an element to the queue.
 */
void *gq_append(gq_t *q, void const *ent);

/**
 * Prepend an element to the queue.
 */
void *gq_prepend(gq_t *q, void const *ent);

/**
 * Remove one or more elements from the end of the queue.
 */
size_t gq_pop(gq_t *q, size_t n_ent);

/**
 * Remove one or more elements from the start of the queue.
 */
size_t gq_shift(gq_t *q, size_t n_ent);

/**
 * Get a pointer to the start of the queue.
 */
void *gq_head_ptr(gq_t *q);

/**
 * Get a pointer to the end of the queue.
 */
void *gq_tail_ptr(gq_t *q);

/**
 * Get a pointer to the Nth element in the queue
 */
void *gq_void(gq_t *q, size_t idx);

/**
 * Get the first element in the queue.
 */
#define gq_head(q,t) (*(t *)gq_head_ptr(q))

/**
 * Get the last element in the queue.
 */
#define gq_tail(q,t) (*(t *)gq_tail_ptr(q))

/**
 * Get the Nth element in the queue.
 */
#define gq_ent(q,t,i) (*(t *)gq_void(q,i))

#ifdef __cplusplus
}
#endif
#endif /* __SPHINXBASE_GQ_H__ */
