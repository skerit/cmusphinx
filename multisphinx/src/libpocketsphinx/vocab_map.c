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
 * @file vocab_map.c
 * @brief Vocabulary mapping
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#include <sphinxbase/pio.h>
#include <sphinxbase/strfuncs.h>
#include <sphinxbase/ckd_alloc.h>

#include "vocab_map.h"

/**
 * Implementation of the vocabulary mapping object.
 */
struct vocab_map_s {
    int refcount;
    dict_t *dict;
};

vocab_map_t *
vocab_map_init(dict_t *dict)
{
    vocab_map_t *vm;
    
    vm = ckd_calloc(1, sizeof(*vm));
    vm->dict = dict_retain(dict);
    vm->refcount = 1;

    return vm;
}

vocab_map_t *
vocab_map_retain(vocab_map_t *vm)
{
    ++vm->refcount;
    return vm;
}

int
vocab_map_free(vocab_map_t *vm)
{
    if (vm == NULL)
        return 0;
    if (--vm->refcount > 0)
        return vm->refcount;
    dict_free(vm->dict);
    ckd_free(vm);
    return 0;
}

int
vocab_map_read(vocab_map_t *vm, FILE *fh)
{
    lineiter_t *li;

    for (li = lineiter_start(fh); li; li = lineiter_next(li)) {
        size_t nwords;
        char **wptr;

        string_trim(li->buf, STRING_BOTH);
        nwords = str2words(li->buf, NULL, 0);
        wptr = ckd_calloc(nwords, sizeof(*wptr));
        str2words(li->buf, wptr, nwords);
        ckd_free(wptr);
    }
    return 0;
}

int
vocab_map_write(vocab_map_t *vm, FILE *fh)
{
    return 0;
}

int32
vocab_map_map(vocab_map_t *vm, int32 wid)
{
    return -1;
}

int32 const *
vocab_map_unmap(vocab_map_t *vm, int32 pseudo_wid,
                             int32 *out_n_mapped)
{
    return NULL;
}
