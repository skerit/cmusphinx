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
 * @file ring_queue.h
 * @brief Circular arrays that support synchronized queue operations.
 */

#ifndef __RING_QUEUE_H__
#define __RING_QUEUE_H__

typedef struct ring_queue_s ring_queue_t;

/**
 * Flags passed to ring_queue_init()
 */
enum ring_queue_flags_e {
	RING_QUEUE_SYNCHRONIZED = (1<<0), /**< Synchronize accesses. */
	RING_QUEUE_GROWABLE = (1<<1)      /**< Grow the queue when full. */
};

/**
 * Create a ring queue.
 *
 * @param n_items Number of items in the queue.
 * @param item_size Size of each item.
 * @param flags Configuration flags.
 * @return Newly allocated queue, or NULL on failure.
 */
SPHINXBASE_EXPORT
ring_queue_t *ring_queue_init(size_t n_items, size_t item_size, int flags);

/**
 * Retain a pointer to a ring queue.
 */
SPHINXBASE_EXPORT
ring_queue_t *ring_queue_retain(ring_queue_t *rq);

/**
 * Release a pointer to a ring queue.
 */
SPHINXBASE_EXPORT
int ring_queue_free(ring_queue_t *rq);




#endif /* __RING_QUEUE_H__ */
