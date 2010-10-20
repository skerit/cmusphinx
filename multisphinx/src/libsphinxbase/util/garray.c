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

#include "garray.h"
#include "err.h"
#include "ckd_alloc.h"

#include <string.h>

struct garray_s {
    int refcount;
    void *ent;
    size_t ent_size;
    size_t n_ent, n_ent_alloc;
};

garray_t *
garray_init(size_t n_ent, size_t ent_size)
{
    garray_t *gar = ckd_calloc(1, sizeof(*gar));
    gar->refcount = 1;
    if (n_ent == 0)
        gar->n_ent_alloc = 8; /* Arbitrary number. */
    else
        gar->n_ent_alloc = n_ent;
    gar->ent_size = ent_size;
    gar->n_ent = n_ent;
    gar->ent = ckd_malloc(gar->n_ent_alloc * gar->ent_size);
    return gar;
}

garray_t *
garray_retain(garray_t *gar)
{
    ++gar->refcount;
    return gar;
}

int
garray_free(garray_t *gar)
{
    if (gar == NULL)
        return 0;
    if (--gar->refcount > 0)
        return gar->refcount;

    ckd_free(gar->ent);
    ckd_free(gar);
    return 0;
}

size_t
garray_size(garray_t *gar)
{
    return gar->n_ent;
}

size_t
garray_alloc_size(garray_t *gar)
{
    return gar->n_ent_alloc;
}

size_t
garray_expand(garray_t *gar, size_t n_ent)
{
    if (n_ent > gar->n_ent_alloc) {
        while (n_ent > gar->n_ent_alloc)
            gar->n_ent_alloc *= 2;
        gar->ent = ckd_realloc(gar->ent, gar->n_ent_alloc * gar->ent_size);
    }
    gar->n_ent = n_ent;
    return gar->n_ent;
}

void *
garray_void(garray_t *gar, size_t idx)
{
    return (char *)gar->ent + idx * gar->ent_size;
}

void *
garray_append(garray_t *gar, void *ent)
{
    garray_expand(gar, gar->n_ent + 1);
    memcpy(garray_void(gar, gar->n_ent - 1),
           ent, gar->ent_size);
    return garray_void(gar, gar->n_ent - 1);
}

size_t
garray_pop(garray_t *gar, size_t n_ent)
{
    if (n_ent > gar->n_ent)
        gar->n_ent = 0;
    else
        gar->n_ent -= n_ent;
    return gar->n_ent;
}

void
garray_reset(garray_t *gar)
{
    gar->n_ent = 0;
}

size_t
garray_shift(garray_t *gar, size_t n_ent)
{
    if (n_ent > gar->n_ent)
        gar->n_ent = 0;
    else
        gar->n_ent -= n_ent;
    if (gar->n_ent == 0)
        return 0;
    memmove(gar->ent, garray_void(gar, n_ent),
            gar->n_ent * gar->ent_size);
    return gar->n_ent;
}

void
garray_clear(garray_t *gar, size_t start, size_t n_ent)
{
    memset(garray_void(gar, start), 0, n_ent * gar->ent_size);
}

garray_t *
garray_slice(garray_t *gar, size_t start, size_t n_ent)
{
    garray_t *gar2;

    if (start + n_ent > gar->n_ent)
        return NULL;
    gar2 = garray_init(n_ent, gar->ent_size);
    memcpy(gar2->ent, garray_void(gar, start), n_ent * gar->ent_size);
    return gar2;
}
