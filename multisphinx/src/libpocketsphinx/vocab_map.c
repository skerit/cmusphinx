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
#include <sphinxbase/garray.h>
#include <sphinxbase/ckd_alloc.h>

#include "vocab_map.h"

/**
 * Implementation of the vocabulary mapping object.
 */
struct vocab_map_s {
    int refcount;
    dict_t *dict;
    int gendict;
    garray_t *pseudos; /**< Map pseudo words to words. */
    garray_t *words;   /**< Map words to pseudo words. */
    garray_t *wids;    /**< Word IDs for word sets. */
};

vocab_map_t *
vocab_map_init(dict_t *dict)
{
    vocab_map_t *vm;
    
    vm = ckd_calloc(1, sizeof(*vm));
    vm->refcount = 1;

    if (dict) {
        vm->dict = dict_retain(dict);
        vm->gendict = FALSE;
    }
    else {
        vm->dict = dict_init(NULL, NULL);
        vm->gendict = TRUE;
    }

    vm->pseudos = garray_init(0, sizeof(i32p_t));
    vm->words = garray_init(0, sizeof(i32p_t));
    vm->wids = garray_init(0, sizeof(int32));

    garray_set_cmp(vm->pseudos, garray_cmp_i32p_first, NULL);
    garray_set_cmp(vm->words, garray_cmp_i32p_first, NULL);

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
    garray_free(vm->pseudos);
    garray_free(vm->words);
    garray_free(vm->wids);
    ckd_free(vm);
    return 0;
}

dict_t *
vocab_map_dict(vocab_map_t *vm)
{
    return vm->dict;
}

int
vocab_map_read(vocab_map_t *vm, FILE *fh)
{
    lineiter_t *li;

    for (li = lineiter_start(fh); li; li = lineiter_next(li)) {
        size_t nwords, i;
        char **wptr;
        i32p_t pseudo;
        int32 nmapped;

        string_trim(li->buf, STRING_BOTH);
        nwords = str2words(li->buf, NULL, 0);
        wptr = ckd_calloc(nwords, sizeof(*wptr));
        str2words(li->buf, wptr, nwords);

        if ((pseudo.a = dict_wordid(vm->dict, wptr[0])) == BAD_S3WID)
            if (vm->gendict)
                pseudo.a = dict_add_word(vm->dict, wptr[0], NULL, 0);
        if (pseudo.a == BAD_S3WID) {
            E_ERROR("Skipping unknown pseudo-word %s\n", wptr[0]);
            goto next_line;
        }
        pseudo.b = garray_next_idx(vm->wids);
        garray_append(vm->pseudos, &pseudo);

        nmapped = 0;
        garray_append(vm->wids, &nmapped);
        for (i = 1; i < nwords; ++i) {
            i32p_t word;

            if ((word.a = dict_wordid(vm->dict, wptr[i])) == BAD_S3WID)
                if (vm->gendict)
                    word.a = dict_add_word(vm->dict, wptr[i], NULL, 0);
            if (word.a == BAD_S3WID) {
                E_ERROR("Skipping unknown word %s\n", wptr[i]);
                continue;
            }
            word.b = pseudo.a;
            garray_append(vm->words, &word);
            garray_append(vm->wids, &word.a);
            ++nmapped;
        }
        garray_put(vm->wids, pseudo.b, &nmapped);
    next_line:
        ckd_free(wptr);
    }
    garray_sort(vm->pseudos);
    garray_sort(vm->words);

    return 0;
}

int
vocab_map_write(vocab_map_t *vm, FILE *fh)
{
    size_t i;

    for (i = 0; i < garray_next_idx(vm->pseudos); ++i) {
        i32p_t *pseudo = garray_ptr(vm->pseudos, i32p_t, i);
        int32 nmapped, i;

        fprintf(fh, "%s", dict_wordstr(vm->dict, pseudo->a));
        nmapped = garray_ent(vm->wids, int32, pseudo->b);
        for (i = 0; i < nmapped; ++i) {
            int32 wid = garray_ent(vm->wids, int32, pseudo->b + i + 1);
            fprintf(fh, " %s", dict_wordstr(vm->dict, wid));
        }
        fprintf(fh, "\n");
    }

    return 0;
}

int32
vocab_map_map(vocab_map_t *vm, int32 wid)
{
    size_t pos;
    i32p_t ent;

    ent.a = wid;
    pos = garray_find_first(vm->words, &ent);
    if (pos == garray_next_idx(vm->words))
        return -1;
    return garray_ptr(vm->words, i32p_t, pos)->b;
}

int32 const *
vocab_map_unmap(vocab_map_t *vm, int32 pseudo_wid,
                int32 *out_n_mapped)
{
    size_t pos;
    i32p_t ent;

    ent.a = pseudo_wid;
    pos = garray_find_first(vm->pseudos, &ent);
    if (pos == garray_next_idx(vm->pseudos))
        return NULL;
    ent = garray_ent(vm->pseudos, i32p_t, pos);
    if (out_n_mapped)
        *out_n_mapped = garray_ent(vm->wids, int32, ent.b);
    return garray_ptr(vm->wids, int32, ent.b + 1);
}

struct vocab_map_iter_s {
    vocab_map_t *vm;
    size_t pos;
};

vocab_map_iter_t *
vocab_map_mappings(vocab_map_t *vm)
{
    vocab_map_iter_t *itor;

    if (garray_next_idx(vm->pseudos) == 0)
        return NULL;

    itor = ckd_calloc(1, sizeof(*itor));
    itor->vm = vm;
    itor->pos = 0;
    return itor;
}

vocab_map_iter_t *
vocab_map_iter_next(vocab_map_iter_t *itor)
{
    if (++itor->pos == garray_next_idx(itor->vm->pseudos)) {
        vocab_map_iter_free(itor);
        return NULL;
    }
    return itor;
}

void
vocab_map_iter_free(vocab_map_iter_t *itor)
{
    ckd_free(itor);
}

int32 const *
vocab_map_iter_get(vocab_map_iter_t *itor,
                   int32 *out_pseudo_wid,
                   int32 *out_n_mapped)
{
    i32p_t ent;

    ent = garray_ent(itor->vm->pseudos, i32p_t, itor->pos);
    if (out_pseudo_wid) *out_pseudo_wid = ent.a;
    if (out_n_mapped)
        *out_n_mapped = garray_ent(itor->vm->wids, int32, ent.b);
    return garray_ptr(itor->vm->wids, int32, ent.b + 1);
}
