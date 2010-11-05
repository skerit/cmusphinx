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

#include <sphinxbase/garray.h>
#include <sphinxbase/listelem_alloc.h>

#define MIN_LOGPROB 1e-20

/**
 * N-Gram trie.
 */
struct ngram_trie_s {
    int refcount;
    dict_t *dict;      /**< Dictionary which maps words to IDs. */
    logmath_t *lmath;  /**< Log-math used for input/output. */
    int shift;         /**< Shift applied internally to log values. */
    int zero;          /**< Minimum allowable log value. */
    int n;             /**< Maximum N-Gram order. */

    ngram_trie_node_t *root;
    listelem_alloc_t *node_alloc;

#if 0 /* Future memory-optimized version... */
    garray_t *nodes;       /**< Flat array of nodes. */
    garray_t *successors;  /**< Flat array of successor pointers. */
#endif
};

/**
 * N-Gram structure
 */
struct ngram_trie_node_s {
    int32 word;
    int16 log_prob;
    int16 log_bowt;
#if 0 /* Future memory-optimized version... */
    int32 history;    /**< Index of parent node. */
    int32 successors; /**< Index of child nodes. */
#endif
    ngram_trie_node_t *history;
    garray_t *successors;
};

/**
 * Iterator over N-Grams in a trie.
 */
struct ngram_trie_iter_s {
    ngram_trie_t *t;        /**< Parent trie. */
    ngram_trie_node_t *cur; /**< Current node (in this successor array). */
    int32 pos;              /**< Position in cur->successors. */
    int nostop;             /**< Continue to next node at same level. */
};

static ngram_trie_node_t *
ngram_trie_alloc_node(ngram_trie_t *t)
{
    ngram_trie_node_t *node;

    node = listelem_malloc(t->node_alloc);
    node->word = -1;
    node->log_prob = 0;
    node->log_bowt = 0;
    node->history = NULL;
    node->successors = garray_init(0, sizeof(ngram_trie_node_t *));

    return node;
}

ngram_trie_t *
ngram_trie_init(dict_t *d, logmath_t *lmath)
{
    ngram_trie_t *t;

    t = ckd_calloc(1, sizeof(*t));
    t->refcount = 1;
    t->dict = dict_retain(d);
    t->lmath = logmath_retain(lmath);

    /* Determine proper shift to fit min_logprob in 16 bits. */
    t->zero = logmath_log(lmath, MIN_LOGPROB);
    while (t->zero < -32768) {
        t->zero >>= 1;
        ++t->shift;
    }

#if 0
    /* Create arrays and add the root node. */
    t->nodes = garray_init(1, sizeof(ngram_trie_node_t));
    t->successors = garray_init(0, sizeof(int32));
    t->root = garray_ptr(t->nodes, ngram_trie_node_t, 0);
    t->root->word = -1;
    t->root->log_prob = 0;
    t->root->log_bowt = 0;
    t->root->history = -1;
    t->root->successors = 0; /* No successors yet. */
#endif

    t->node_alloc = listelem_alloc_init(sizeof(ngram_trie_node_t));
    t->root = ngram_trie_alloc_node(t);

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
#if 0
    garray_free(t->nodes);
    garray_free(t->successors);
#endif
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
ngram_trie_root(ngram_trie_t *t)
{
    return t->root;
}

ngram_trie_node_t *
ngram_trie_ngram(ngram_trie_t *t, char const *w, ...)
{
    ngram_trie_node_t *node;
    char const *h;
    va_list args;
    int n_hist;
    int32 *hist;
    int32 wid;

    wid = dict_wordid(t->dict, w);
    va_start(args, w);
    n_hist = 0;
    while ((h = va_arg(args, char const *)) != NULL)
        ++n_hist;
    va_end(args);
    hist = ckd_calloc(n_hist, sizeof(*hist));
    va_start(args, w);
    n_hist = 0;
    while ((h = va_arg(args, char const *)) != NULL) {
        hist[n_hist] = dict_wordid(t->dict, h);
        ++n_hist;
    }
    va_end(args);

    node = ngram_trie_ngram_v(t, wid, hist, n_hist);
    ckd_free(hist);
    return node;
}

ngram_trie_node_t *
ngram_trie_ngram_v(ngram_trie_t *t, int32 w,
                   int32 *hist, int32 n_hist)
{
    ngram_trie_node_t *node;

    node = ngram_trie_root(t);
    if (n_hist > t->n - 1)
        n_hist = t->n - 1;
    while (n_hist > 0) {
        int32 nextwid = hist[n_hist - 1];
        if ((node = ngram_trie_successor(t, node, nextwid)) == NULL)
            return NULL;
    }

    return ngram_trie_successor(t, node, w);
}

int32
ngram_trie_prob(ngram_trie_t *t, int *n_used, char const *w, ...)
{
    char const *h;
    va_list args;
    int n_hist;
    int32 *hist;
    int32 wid;
    int32 prob;

    wid = dict_wordid(t->dict, w);
    va_start(args, w);
    n_hist = 0;
    while ((h = va_arg(args, char const *)) != NULL)
        ++n_hist;
    va_end(args);
    hist = ckd_calloc(n_hist, sizeof(*hist));
    va_start(args, w);
    n_hist = 0;
    while ((h = va_arg(args, char const *)) != NULL) {
        hist[n_hist] = dict_wordid(t->dict, h);
        ++n_hist;
    }
    va_end(args);

    prob = ngram_trie_prob_v(t, n_used, wid, hist, n_hist);
    ckd_free(hist);
    return prob;
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
    ngram_trie_iter_t *itor;
    ngram_trie_node_t *h;
    int i;

    /* Find first N-1-Gram */
    h = t->root;
    for (i = 1; i < n; ++i) {
        h = garray_ent(h->successors, ngram_trie_node_t *, 0);
        if (h == NULL)
            return NULL;
    }

    /* Create an iterator with nostop=TRUE */
    itor = ckd_calloc(1, sizeof(*itor));
    itor->cur = h;
    itor->pos = 0;
    itor->nostop = FALSE;

    return itor;
}

ngram_trie_iter_t *
ngram_trie_successors(ngram_trie_t *t, ngram_trie_node_t *h)
{
    ngram_trie_iter_t *itor;

    /* Create an iterator with nostop=FALSE */
    itor = ckd_calloc(1, sizeof(*itor));
    itor->cur = h;
    itor->pos = 0;
    itor->nostop = FALSE;

    return itor;
}

void
ngram_trie_iter_free(ngram_trie_iter_t *itor)
{
    ckd_free(itor);
}

static ngram_trie_node_t *
ngram_trie_next_node(ngram_trie_t *t, ngram_trie_node_t *ng)
{
    ngram_trie_node_t *h = ng->history;
    size_t pos;

    if (h == NULL)
        return NULL;
    /* Locate ng in h->successors. */
    pos = garray_bisect_left(h->successors, &ng);
    assert(pos < garray_next_idx(h->successors));
    assert(ng == garray_ent(h->successors, ngram_trie_node_t *, pos));
    ++pos;
    if (pos == garray_next_idx(h->successors)) {
        h = ngram_trie_next_node(t, h);
        if (h == NULL)
            return NULL;
        return garray_ent(h->successors, ngram_trie_node_t *, 0);
    }
    else
        return garray_ent(h->successors, ngram_trie_node_t *, pos);
}

ngram_trie_iter_t *
ngram_trie_iter_next(ngram_trie_iter_t *itor)
{
    ++itor->pos;
    if (itor->pos >= garray_next_idx(itor->cur->successors)) {
        if (itor->nostop) {
            itor->cur = ngram_trie_next_node(itor->t, itor->cur);
            if (itor->cur == NULL) {
                ngram_trie_iter_free(itor);
                return NULL;
            }
            itor->pos = 0;
        }
        else  {
            ngram_trie_iter_free(itor);
            return NULL;
        }
    }
    return itor;
}

ngram_trie_iter_t *
ngram_trie_iter_up(ngram_trie_iter_t *itor)
{
    itor->cur = itor->cur->history;
    if (itor->cur == NULL) {
        ngram_trie_iter_free(itor);
        return NULL;
    }
    return itor;
}

ngram_trie_iter_t *
ngram_trie_iter_down(ngram_trie_iter_t *itor)
{
    itor->cur = garray_ent(itor->cur->successors,
                           ngram_trie_node_t *, itor->pos);
    assert(itor->cur != NULL);
    if (garray_next_idx(itor->cur->successors) == 0) {
        ngram_trie_iter_free(itor);
        return NULL;
    }
    itor->pos = 0;
    return itor;
}

ngram_trie_node_t *
ngram_trie_iter_get(ngram_trie_iter_t *itor)
{
    if (itor->pos >= garray_next_idx(itor->cur->successors))
        return NULL;
    else
        return garray_ent(itor->cur->successors,
                          ngram_trie_node_t *, itor->pos);
}

ngram_trie_node_t *
ngram_trie_iter_get_parent(ngram_trie_iter_t *itor)
{
    return itor->cur;
}

static size_t
ngram_trie_successor_pos(ngram_trie_t *t, ngram_trie_node_t *h, int32 w)
{
    ngram_trie_node_t *png, ng;
    ng.word = w;
    png = &ng;
    return garray_bisect_left(h->successors, &png);
}

ngram_trie_node_t *
ngram_trie_successor(ngram_trie_t *t, ngram_trie_node_t *h, int32 w)
{
    size_t pos;

    pos = ngram_trie_successor_pos(t, h, w);
    if (pos >= garray_next_idx(h->successors))
        return NULL;
    return garray_ent(h->successors, ngram_trie_node_t *, pos);
}

int
ngram_trie_delete_successor(ngram_trie_t *t, ngram_trie_node_t *h, int32 w)
{
    ngram_trie_node_t *ng;
    size_t pos;

    /* Bisect the successor array. */
    pos = ngram_trie_successor_pos(t, h, w);
    if (pos >= garray_next_idx(h->successors))
        return -1;
    ng = garray_ent(h->successors, ngram_trie_node_t *, pos);
    /* Delete it. */
    listelem_free(t->node_alloc, ng);
    garray_delete(h->successors, pos, pos+1);
    return 0;
}

ngram_trie_node_t *
ngram_trie_add_successor(ngram_trie_t *t, ngram_trie_node_t *h, int32 w)
{
    ngram_trie_node_t *ng;

    /* Create a new node. */
    ng = listelem_malloc(t->node_alloc);

    /* Bisect the successor array. */

    /* Insert it.*/
    return NULL;
}

int
ngram_trie_add_successor_ngram(ngram_trie_t *t,
                               ngram_trie_node_t *h,
                               ngram_trie_node_t *w)
{
    /* Bisect the successor array. */

    /* Update w->history. */

    /* Insert it.*/
    return 0;
}

ngram_trie_node_t *
ngram_trie_backoff(ngram_trie_t *t,
                   ngram_trie_node_t *ng)
{
    /* Extract word IDs from ng's history. */

    /* Look up the backoff N-Gram. */
    return NULL;
}

int32
ngram_trie_successor_prob(ngram_trie_t *t,
                          ngram_trie_node_t *h, int32 w)
{
    /* Extract word IDs from ng's history. */
    /* Call ngram_trie_prob_v */
    return 0;
}

int32
ngram_trie_calc_bowt(ngram_trie_t *t, ngram_trie_node_t *h)
{
    return 0;
}
