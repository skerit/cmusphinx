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
 * \file gq.c
 * \brief Generic double-ended queue
 */

#include <assert.h>
#include <string.h>

#include <sphinxbase/gq.h>
#include <sphinxbase/err.h>
#include <sphinxbase/ckd_alloc.h>

#include "gq.h"
#include "garray_internal.h"

struct gq_s {
    garray_t base;
    size_t head;
    size_t count;
    int order;
};

gq_t *
gq_init(size_t ent_size)
{
    gq_t *q = ckd_calloc(1, sizeof(*q));
    q->order = 3;
    garray_setup(&q->base, 1<<q->order, ent_size);
    q->head = gq_alloc_size(q) / 2;
    return q;
}

gq_t *
gq_retain(gq_t *q)
{
    return (gq_t *)garray_retain(&q->base);
}

int
gq_free(gq_t *q)
{
    return garray_free(&q->base);
}

size_t
gq_index(gq_t *q, size_t idx)
{
    return idx & ((1 << q->order) - 1);
}

size_t
gq_size(gq_t *q)
{
    return q->count;
}

size_t
gq_ent_size(gq_t *q)
{
    return garray_ent_size(&q->base);
}

size_t
gq_alloc_size(gq_t *q)
{
    return 1 << q->order;
}

static size_t
gq_expand(gq_t *q)
{
    size_t prev = gq_alloc_size(q);
    size_t next;

    ++q->order;
    next = gq_alloc_size(q);
    garray_expand(&q->base, next);

    /* Push back head of queue if necessary. */
    if (q->head + q->count > prev) {
	size_t n_move = prev - q->head;
	garray_move(&q->base, next - n_move, q->head, n_move);
	q->head = next - n_move;
    }

    return next;
}

void *
gq_append(gq_t *q, void const *ent)
{
    void *dest;

    if (q->count == gq_alloc_size(q))
	gq_expand(q);
    ++q->count;
    E_DEBUG(2,("tail %d count %d size %d\n",
	       gq_index(q, q->head + q->count - 1), q->count, gq_alloc_size(q)));
    dest = gq_tail_ptr(q);
    memcpy(dest, ent, garray_ent_size(&q->base));
    return dest;
}

void *
gq_prepend(gq_t *q, void const *ent)
{
    void *dest;

    if (q->count == gq_alloc_size(q))
	gq_expand(q);
    ++q->count;
    if (q->head == 0)
	q->head = gq_alloc_size(q) - 1;
    else
	--q->head;
    E_DEBUG(2,("head %d count %d size %d\n",
	       q->head, q->count, gq_alloc_size(q)));
    dest = gq_head_ptr(q);
    memcpy(dest, ent, garray_ent_size(&q->base));
    return dest;
}

size_t
gq_pop(gq_t *q, size_t n_ent)
{
    if (n_ent > q->count)
	n_ent = q->count;
    q->count -= n_ent;
    return q->count;
}

size_t
gq_shift(gq_t *q, size_t n_ent)
{
    if (n_ent > q->count)
	n_ent = q->count;
    q->count -= n_ent;
    q->head += n_ent;
    if (q->head >= gq_alloc_size(q))
	q->head -= gq_alloc_size(q);
    return q->count;
}

void *
gq_head_ptr(gq_t *q)
{
    if (q->count == 0)
	return NULL;
    return garray_void(&q->base, q->head);
}

void *
gq_tail_ptr(gq_t *q)
{
    if (q->count == 0)
	return NULL;
    return garray_void(&q->base, gq_index(q, q->head + q->count - 1));
}

void *
gq_void(gq_t *q, size_t idx)
{
    if (idx >= q->count)
	return NULL;
    return garray_void(&q->base, gq_index(q, idx));
}
