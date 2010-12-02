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
 *
 * FIXME: Using these as double-ended queues is suboptimal, but we do
 * it anyway. See <gq.h> for an alternative.
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

#define GARRAY_INVALID_INDEX ((size_t)-1)

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
 * Get the size of an element in a garray_t
 */
size_t garray_ent_size(garray_t *gar);

/**
 * Get the next absolute index in a garray_t
 */
size_t garray_next_idx(garray_t *gar);

/**
 * Get the number of elements allocated in a garray_t
 */
size_t garray_alloc_size(garray_t *gar);

/**
 * Reserve space for elements in an garray_t
 */
size_t garray_reserve(garray_t *gar, size_t n_ent);

/**
 * Extend the array to contain the desired number of elements
 */
size_t garray_expand(garray_t *gar, size_t n_ent);

/**
 * Extend the array up to (but not including) a given absolute index.
 */
size_t garray_expand_to(garray_t *gar, size_t next_idx);

/**
 * Append an element to the array, resizing if necessary.
 */
void *garray_append(garray_t *gar, void const *ent);

/**
 * Insert an element in the array.
 *
 * If the specified index is before the base index of the array, this
 * operation will fail.
 *
 * To ease prototyping things in Python, if the specified index is
 * past the end of the array, this will simply append @a ent rather
 * than failing.
 */
void *garray_insert(garray_t *gar, size_t idx, void const *ent);

/**
 * Delete one or more elements from an array.
 *
 * If @a end is past the end of the array, this will fail - to
 * truncate an array use garray_pop() or garray_pop_from().
 */
int garray_delete(garray_t *gar, size_t start, size_t end);

/**
 * Set an element in the array.
 *
 * If the specified index is past the end of the array or before the
 * base index of the array, this operation will fail.
 */
void *garray_put(garray_t *gar, size_t idx, void const *ent);

/**
 * Comparison function for bisect operations.
 */
typedef int (*garray_cmp_t)(garray_t *gar,
                            void const *a, void const *b,
                            void *udata);

/**
 * Standard comparison function for 32-bit integers.
 */
int garray_cmp_int32(garray_t *gar, void const *a, void const *b, void *udata);

/**
 * Standard comparison function for strings.
 */
int garray_cmp_str(garray_t *gar, void const *a, void const *b, void *udata);

/**
 * Standard comparison function for integer pairs.
 */
int garray_cmp_i32p(garray_t *gar, void const *a, void const *b, void *udata);

/**
 * Standard comparison function for the first element of integer pairs
 */
int garray_cmp_i32p_first(garray_t *gar, void const *a,
                          void const *b, void *udata);

/**
 * Set the comparison function to be used for bisect operations.
 */
void garray_set_cmp(garray_t *gar, garray_cmp_t cmp, void *udata);

/**
 * Find the leftmost position for inserting an element in a sorted array.
 */
size_t garray_bisect_left(garray_t *gar, void *ent);

/**
 * Find the rightmost position for inserting an element in a sorted array.
 *
 * It is recommended to use this for element insertion since this may
 * result in less memory being moved.
 */
size_t garray_bisect_right(garray_t *gar, void *ent);

/**
 * Find the first matching element in a sorted array.
 *
 * If no matching element is found, this function will return the next
 * available index in the array, as returned by garray_next_idx()
 */
size_t garray_find_first(garray_t *gar, void *ent);

/**
 * Heapify an array in-place.
 *
 * Note that contrary to the Python version, this creates a max-heap,
 * where gar[0] is the greatest value in the heap.
 */
void garray_heapify(garray_t *gar);

/**
 * Sort an array in-place using heapsort.
 */
void garray_sort(garray_t *gar);

/**
 * Sort an array in-place using mergesort.
 */
void garray_mergesort(garray_t *gar);


/**
 * Remove elements from the end of the array.
 */
size_t garray_pop(garray_t *gar, size_t n_ent);

/**
 * Remove elements from an absolute index to the end of the array.
 */
size_t garray_pop_from(garray_t *gar, size_t first_idx);

/**
 * Remove elements from the start of the array.
 */
size_t garray_shift(garray_t *gar, size_t n_ent);

/**
 * Remove elements up to an absolute index in the array.
 */
size_t garray_shift_from(garray_t *gar, size_t first_idx);

/**
 * Set the base index of the array.
 *
 * This sets a base index which will be subtracted from all indexing
 * operations into the array.
 *
 * You could use this to emulate FORTRAN but that's not actually what
 * it's for.  Instead, it is used for tables where the array index
 * corresponds to a non-decreasing series of identifers, such as
 * backpointer indices.
 */
size_t garray_set_base(garray_t *gar, size_t base_idx);

/**
 * Get the base index of the array.
 */
size_t garray_base(garray_t *gar);

/**
 * Remove all elements from the array.
 */
void garray_reset(garray_t *gar);

/**
 * Remove all elements from the array and update its base index.
 */
void garray_reset_to(garray_t *gar, size_t base_idx);

/**
 * Set memory for elements to zero.
 */
void garray_clear(garray_t *gar, size_t start, size_t n_ent);

/**
 * Copy a subsection of the array to another array.
 */
garray_t *garray_slice(garray_t *gar, size_t start, size_t n_ent);

/**
 * Move a subsection of the array to another position.
 */
size_t garray_move(garray_t *gar, size_t dest, size_t src, size_t n_ent);

/**
 * Get a pointer to an element in the array.
 */
void *garray_void(garray_t *gar, size_t idx);

/**
 * Get the element index of a pointer in the array.
 */
size_t garray_idx(garray_t *gar, void *ent);

/**
 * Get the underlying data store.
 */
#define garray_ptr(gar, type, idx) ((type *)garray_void(gar, idx))

/**
 * Get an individual element in the array.
 */
#define garray_ent(gar, type, idx) (*garray_ptr(gar, type, idx))

#ifdef __cplusplus
}
#endif
#endif /* __SPHINXBASE_GARRAY_H__ */
