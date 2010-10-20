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
 * \file garray.h
 * \brief Generic expandable arrays.
 */

#ifndef __SPHINXBASE_GARRAY_H__
#define __SPHINXBASE_GARRAY_H__

#include <stdlib.h>
/* Win32/WinCE DLL gunk */
#include <sphinxbase/sphinxbase_export.h>
#include <sphinxbase/prim_type.h>

#ifdef __cplusplus
extern "C" {
#endif
#if 0
/* Fool Emacs. */
}
#endif

typedef struct garray_s garray_t;

/**
 * Create a new garray_t of the given size.
 */
garray_t *garray_init(size_t n_ent, size_t ent_size);

/**
 * Retain a pointer to a garray_t
 */
garray_t *garray_retain(garray_t *gar);

/**
 * Release a pointer to a garray_t
 */
int garray_free(garray_t *gar);

/**
 * Get the number of elements in a garray_t
 */
size_t garray_size(garray_t *gar);

/**
 * Get the number of elements allocated in a garray_t
 */
size_t garray_alloc_size(garray_t *gar);

/**
 * Extend the array to contain the desired number of elements
 */
size_t garray_expand(garray_t *gar, size_t n_ent);

/**
 * Append an element to the array, resizing if necessary.
 */
void *garray_append(garray_t *gar, void *ent);

/**
 * Remove elements from the end of the array.
 */
size_t garray_pop(garray_t *gar, size_t n_ent);

/**
 * Remove elements from the start of the array.
 */
size_t garray_shift(garray_t *gar, size_t n_ent);

/**
 * Remove all elements from the array.
 */
void garray_reset(garray_t *gar);

/**
 * Set memory for elements to zero.
 */
void garray_clear(garray_t *gar, size_t start, size_t n_ent);

/**
 * Copy a subsection of the array to another array.
 */
garray_t *garray_slice(garray_t *gar, size_t start, size_t n_ent);

/**
 * Get a pointer to an element in the array.
 */
void *garray_void(garray_t *gar, size_t idx);

/**
 * Get the underlying data store.
 */
#define garray_ptr(gar, type, idx) ((type *)garray_void(gar, idx))

/**
 * Get an individual element in the array.
 */
#define garray_ent(gar, type, idx) garray_ptr(gar, type, 0)[idx]

#ifdef __cplusplus
}
#endif
#endif /* __SPHINXBASE_GARRAY_H__ */
