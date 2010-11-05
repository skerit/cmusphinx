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
 * @file ngram_trie.c
 * @brief Mutable trie implementation of N-Gram language models
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#include "ngram_trie.h"

struct ngram_trie_s {
    int refcount;
    dict_t *dict;
    logmath_t *lmath;
};

struct ngram_trie_node_s {
};

struct ngram_trie_iter_s {
};

ngram_trie_t *
ngram_trie_init(dict_t *d, logmath_t *lmath)
{
    ngram_trie_t *t;

    t = ckd_calloc(1, sizeof(*t));
    t->refcount = 1;
    t->dict = dict_retain(d);
    t->lmath = logmath_retain(lmath);

    return t;
}

ngram_trie_t *
ngram_trie_retain(ngram_trie_t *t)
{
    ++t->refcount;
    return t;
}

int
ngram_trie_free(ngram_trie_t *t)
{
    if (t == NULL)
        return 0;
    if (--t->refcount > 0)
        return t->refcount;
    dict_free(t->dict);
    logmath_free(t->lmath);
    ckd_free(t);
    return 0;
}

int
ngram_trie_read_arpa(ngram_trie_t *t, FILE *arpafile)
{
    return 0;
}

int
ngram_trie_write_arpa(ngram_trie_t *t, FILE *arpafile)
{
    return 0;
}

ngram_trie_node_t *
ngram_trie_ngram(ngram_trie_t *t, int32 w, ...)
{
    return NULL;
}

ngram_trie_node_t *
ngram_trie_ngram_v(ngram_trie_t *t, int32 w,
                   int32 *hist, int32 n_hist)
{
    return NULL;
}

int32
ngram_trie_prob(ngram_trie_t *t, int *n_used, int32 w, ...)
{
    return 0;
}

int32
ngram_trie_prob_v(ngram_trie_t *t, int *n_used, int32 w,
                  int32 *hist, int32 n_hist)
{
    return 0;
}

ngram_trie_iter_t *
ngram_trie_ngrams(ngram_trie_t *t, int n)
{
    return NULL;
}

ngram_trie_iter_t *
ngram_trie_successors(ngram_trie_node_t *h)
{
    return NULL;
}

ngram_trie_node_t *
ngram_trie_successor(ngram_trie_node_t *h, int32 w)
{
    return NULL;
}

int
ngram_trie_delete_successor(ngram_trie_node_t *h, int32 w)
{
    return 0;
}

ngram_trie_node_t *
ngram_trie_add_successor(ngram_trie_node_t *h, int32 w)
{
    return NULL;
}

int
ngram_trie_add_successor_ngram(ngram_trie_node_t *h,
                               ngram_trie_node_t *w)
{
    return 0;
}

ngram_trie_node_t *
ngram_trie_backoff(ngram_trie_t *t,
                   ngram_trie_node_t *ng)
{
    return NULL;
}

int32
ngram_trie_successor_prob(ngram_trie_t *t,
                          ngram_trie_node_t *h, int32 w)
{
    return 0;
}

int32
ngram_trie_calc_bowt(ngram_trie_t *t, ngram_trie_node_t *h)
{
    return 0;
}
