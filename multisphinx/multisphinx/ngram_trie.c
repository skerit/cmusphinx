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

#include <string.h>
#include <ctype.h>
#include <math.h>

#include <sphinxbase/pio.h>
#include <sphinxbase/garray.h>
#include <sphinxbase/gq.h>
#include <sphinxbase/strfuncs.h>
#include <sphinxbase/listelem_alloc.h>

#include "ngram_trie.h"

#define MIN_LOGPROB 1e-20

/**
 * N-Gram trie.
 */
struct ngram_trie_s {
    int refcount;
    dict_t *dict;      /**< Dictionary which maps words to IDs. */
    int gendict;       /**< Is the dictionary generated from the unigram? */
    logmath_t *lmath;  /**< Log-math used for input/output. */
    int shift;         /**< Shift applied internally to log values. */
    int zero;          /**< Minimum allowable log value. */
    int n;             /**< Maximum N-Gram order. */
    garray_t *counts;  /**< N-Gram counts. */

    int32 start_wid;   /**< Word ID for special word <s> (a non-event) */
    int32 finish_wid;  /**< Word ID for special word </s> (has no valid
                          successors) */

    ngram_trie_node_t *root;
    listelem_alloc_t *node_alloc;
};

/**
 * N-Gram structure
 */
struct ngram_trie_node_s {
    int32 word;
    int16 log_prob;
    int16 log_bowt;
    ngram_trie_node_t *history;
    ngram_trie_node_t *backoff;
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

static int
ngram_trie_nodeptr_cmp(garray_t *gar, void const *a, void const *b, void *udata)
{
    ngram_trie_node_t *na, *nb;
    ngram_trie_t *t;

    na = *(ngram_trie_node_t **)a;
    nb = *(ngram_trie_node_t **)b;
    t = (ngram_trie_t *)udata;

    return strcmp(dict_wordstr(t->dict, na->word),
                  dict_wordstr(t->dict, nb->word));
}

ngram_trie_node_t *
ngram_trie_node_alloc(ngram_trie_t *t)
{
    ngram_trie_node_t *node;

    node = listelem_malloc(t->node_alloc);
    node->word = -1;
    node->log_prob = 0;
    node->log_bowt = 0;
    node->history = NULL;
    node->successors = NULL; /* Allocate lazily */
    node->backoff = (ngram_trie_node_t *)-1; /* Calculate lazily */
    return node;
}

void
ngram_trie_node_free(ngram_trie_t *t, ngram_trie_node_t *node)
{
    garray_free(node->successors);
    listelem_free(t->node_alloc, node);
}

ngram_trie_t *
ngram_trie_init(dict_t *d, logmath_t *lmath)
{
    ngram_trie_t *t;

    t = ckd_calloc(1, sizeof(*t));
    t->refcount = 1;
    if (d) {
        t->dict = dict_retain(d);
        t->gendict = FALSE;
    }
    else {
        t->dict = dict_init(NULL, NULL);
        t->gendict = TRUE;
    }
    t->start_wid = dict_wordid(t->dict, S3_START_WORD);
    t->finish_wid = dict_wordid(t->dict, S3_FINISH_WORD);
    assert(t->finish_wid != BAD_S3WID);
    t->lmath = logmath_retain(lmath);

    /* Determine proper shift to fit min_logprob in 16 bits. */
    t->zero = logmath_log(lmath, MIN_LOGPROB);
    while (t->zero < -32768) {
        t->zero >>= 1;
        ++t->shift;
    }

    t->counts = garray_init(0, sizeof(int));
    t->node_alloc = listelem_alloc_init(sizeof(ngram_trie_node_t));
    t->root = ngram_trie_node_alloc(t);

    return t;
}

ngram_trie_t *
ngram_trie_retain(ngram_trie_t *t)
{
    ++t->refcount;
    return t;
}

static void
free_successor_arrays(ngram_trie_node_t *node)
{
    size_t pos;

    if (node->successors == NULL)
        return;
    for (pos = 0; pos < garray_size(node->successors); ++pos)
        free_successor_arrays(garray_ent(node->successors, ngram_trie_node_t *, pos));
    garray_free(node->successors);
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
    free_successor_arrays(t->root);
    listelem_alloc_free(t->node_alloc);
    garray_free(t->counts);
    ckd_free(t);
    return 0;
}

dict_t *
ngram_trie_dict(ngram_trie_t *t)
{
    return t->dict;
}

logmath_t *
ngram_trie_logmath(ngram_trie_t *t)
{
    return t->lmath;
}

int32
ngram_trie_zero(ngram_trie_t *t)
{
    return t->zero;
}


ngram_trie_node_t *
ngram_trie_root(ngram_trie_t *t)
{
    return t->root;
}

int
ngram_trie_n(ngram_trie_t *t)
{
    return t->n;
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
                   int32 const *hist, int32 n_hist)
{
    ngram_trie_node_t *node;

    node = ngram_trie_root(t);
    if (n_hist > t->n - 1)
        n_hist = t->n - 1;

#if 0
    E_INFO("Looking for N-Gram %s |",
           dict_wordstr(t->dict, w));
    int i;
    for (i = 0; i < n_hist; ++i) {
        E_INFOCONT(" %s", dict_wordstr(t->dict, hist[n_hist - 1 - i]));
    }
    E_INFOCONT("\n");
#endif

    while (n_hist > 0) {
        int32 nextwid = hist[n_hist - 1];
        if ((node = ngram_trie_successor(t, node, nextwid)) == NULL)
            return NULL;
        --n_hist;
    }

    return ngram_trie_successor(t, node, w);
}

ngram_trie_node_t *
ngram_trie_ngram_init(ngram_trie_t *t, char const *w, ...)
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

    node = ngram_trie_ngram_init_v(t, wid, hist, n_hist);
    ckd_free(hist);
    return node;
}

ngram_trie_node_t *
ngram_trie_ngram_init_v(ngram_trie_t *t, int32 w,
                        int32 const *hist, int32 n_hist)
{
    ngram_trie_node_t *node, *nextnode;

    node = ngram_trie_root(t);
    if (n_hist > t->n - 1)
        t->n = n_hist + 1;

#if 0
    E_INFO("Adding N-Gram %s |",
           dict_wordstr(t->dict, w));
    int i;
    for (i = 0; i < n_hist; ++i) {
        E_INFOCONT(" %s", dict_wordstr(t->dict, hist[n_hist - 1 - i]));
    }
    E_INFOCONT("\n");
#endif

    while (n_hist > 0) {
        int32 nextwid = hist[n_hist - 1];
        if ((nextnode = ngram_trie_successor(t, node, nextwid)) == NULL)
            nextnode = ngram_trie_add_successor(t, node, nextwid);
        node = nextnode;
        --n_hist;
    }

    if ((nextnode = ngram_trie_successor(t, node, w)) == NULL)
        nextnode = ngram_trie_add_successor(t, node, w);
    return nextnode;
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

ngram_trie_iter_t *
ngram_trie_ngrams(ngram_trie_t *t, int n)
{
    ngram_trie_iter_t *itor;
    ngram_trie_node_t *cur = NULL;
    gq_t *q;

    /* Depth-first search (preorder traversal) in trie to find first
     * N-1-Gram that has N-Gram successors.  FIXME: Should just have
     * generic trie traversal functions as they will come in handy. */
    q = gq_init(sizeof(ngram_trie_node_t *));
    /* Note that we are actually using q as a stack.  It does not care. */
    gq_append(q, &t->root);
    while (gq_size(q)) {
        ngram_trie_node_t *h = gq_tail(q, ngram_trie_node_t *);
        size_t i;
        gq_pop(q, 1);
        if (h->successors == NULL)
            continue;
        if (ngram_trie_node_n(t, h) == n - 1) {
            cur = h;
            break;
        }
        for (i = garray_size(h->successors); i > 0; --i)
            gq_append(q, garray_void(h->successors, i-1));
    }
    gq_free(q);

    if (cur == NULL)
        return NULL;

    /* Create an iterator with nostop=TRUE */
    itor = ckd_calloc(1, sizeof(*itor));
    itor->t = t;
    itor->cur = cur;
    itor->pos = 0;
    itor->nostop = TRUE;

    assert(itor->cur->successors);

    return itor;
}

ngram_trie_node_t **
ngram_trie_successors_unchecked(ngram_trie_t *t, ngram_trie_node_t *h, size_t *out_nsucc)
{
    size_t nsucc;
    if (h->successors == NULL
        || (nsucc = garray_size(h->successors)) == 0) {
        if (out_nsucc)
            *out_nsucc = 0;
        return NULL;
    }
    if (out_nsucc)
        *out_nsucc = nsucc;
    return garray_void(h->successors, 0);
}

ngram_trie_iter_t *
ngram_trie_successors(ngram_trie_t *t, ngram_trie_node_t *h)
{
    ngram_trie_iter_t *itor;

    if (h->successors == NULL || garray_size(h->successors) == 0)
        return NULL;

    /* Create an iterator with nostop=FALSE */
    itor = ckd_calloc(1, sizeof(*itor));
    itor->t = t;
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
    assert(h->successors != NULL);
    pos = garray_bisect_left(h->successors, &ng);
    assert(pos < garray_next_idx(h->successors));
    /* Duplicates are not allowed, make it more obvious than just a
     * fairly meaningless assert. */
    if (ng != garray_ent(h->successors, ngram_trie_node_t *, pos)) {
        E_ERROR("Duplicate nodes for word %d = %s\n",
                ng->word, dict_wordstr(t->dict, ng->word));
        assert(FALSE);
    }
    ++pos;
    if (pos == garray_next_idx(h->successors)) {
        h = ngram_trie_next_node(t, h);
        if (h == NULL)
            return NULL;
        while (h->successors == NULL) {
            h = ngram_trie_next_node(t, h);
            if (h == NULL)
                return NULL;
        }
        h = garray_ent(h->successors, ngram_trie_node_t *, 0);
        return h;
    }
    else {
        h = garray_ent(h->successors, ngram_trie_node_t *, pos);
        return h;
    }
}

ngram_trie_iter_t *
ngram_trie_iter_next(ngram_trie_iter_t *itor)
{
    ++itor->pos;
    assert(itor->cur != NULL);
    if (itor->pos >= garray_next_idx(itor->cur->successors)) {
        if (itor->nostop) {
            itor->cur = ngram_trie_next_node(itor->t, itor->cur);
            if (itor->cur == NULL) {
                ngram_trie_iter_free(itor);
                return NULL;
            }
            /* Have to advance to the next one with successors. */
            while (itor->cur->successors == NULL) {
                itor->cur = ngram_trie_next_node(itor->t, itor->cur);
                if (itor->cur == NULL) {
                    ngram_trie_iter_free(itor);
                    return NULL;
                }
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
    assert(itor->cur != NULL);
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

int32
ngram_trie_node_word(ngram_trie_t *t, ngram_trie_node_t *node)
{
    return node->word;
}

void
ngram_trie_node_set_word(ngram_trie_t *t, ngram_trie_node_t *node, int32 wid)
{
    node->word = wid;
}

void
ngram_trie_node_params(ngram_trie_t *t,
                       ngram_trie_node_t *node,
                       int32 *out_log_prob,
                       int32 *out_log_bowt)
{
    assert(node != NULL);
    if (out_log_prob) *out_log_prob = node->log_prob << t->shift;
    if (out_log_bowt) *out_log_bowt = node->log_bowt << t->shift;
}

void
ngram_trie_node_set_params(ngram_trie_t *t,
                           ngram_trie_node_t *node,
                           int32 log_prob,
                           int32 log_bowt)
{
    assert(node != NULL);
    node->log_prob = log_prob >> t->shift;
    node->log_bowt = log_bowt >> t->shift;
}

void
ngram_trie_node_params_raw(ngram_trie_t *t,
                           ngram_trie_node_t *node,
                           int16 *out_log_prob,
                           int16 *out_log_bowt)
{
    assert(node != NULL);
    if (out_log_prob) *out_log_prob = node->log_prob;
    if (out_log_bowt) *out_log_bowt = node->log_bowt;
}

void
ngram_trie_node_set_params_raw(ngram_trie_t *t,
                               ngram_trie_node_t *node,
                               int16 log_prob,
                               int16 log_bowt)
{
    assert(node != NULL);
    node->log_prob = log_prob;
    node->log_bowt = log_bowt;
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
    ngram_trie_node_t *node;
    size_t pos;

#if 0
    E_INFO("Looking for successor %s in node with head word %s\n",
           dict_wordstr(t->dict, w), h->word == -1 
           ? "<root>" : dict_wordstr(t->dict, h->word));
#endif
    if (h->successors == NULL)
        return NULL;
    pos = ngram_trie_successor_pos(t, h, w);
    if (pos >= garray_next_idx(h->successors))
        return NULL;
    node = garray_ent(h->successors, ngram_trie_node_t *, pos);
    if (node->word != w)
        return NULL;
    return node;
}

int
ngram_trie_delete_successor(ngram_trie_t *t, ngram_trie_node_t *h, int32 w)
{
    ngram_trie_node_t *ng;
    size_t pos;

    if (h->successors == NULL)
        return -1;
    /* Bisect the successor array. */
    pos = ngram_trie_successor_pos(t, h, w);
    if (pos >= garray_next_idx(h->successors))
        return -1;
    ng = garray_ent(h->successors, ngram_trie_node_t *, pos);
    /* Delete it. */
    ngram_trie_node_free(t, ng);
    if (garray_delete(h->successors, pos, pos+1) < 0)
        return -1;

    return 0;
}

ngram_trie_node_t *
ngram_trie_add_successor(ngram_trie_t *t, ngram_trie_node_t *h, int32 w)
{
    ngram_trie_node_t *ng;
    int n;

    ng = ngram_trie_node_alloc(t);
    ng->word = w;
    ng->history = h;
    assert(ng->word >= 0);
    if (h->successors == NULL) {
        h->successors = garray_init(1, sizeof(ngram_trie_node_t *));
        garray_set_cmp(h->successors, &ngram_trie_nodeptr_cmp, t);
        garray_put(h->successors, 0, &ng);
    }
    else {
        size_t pos;
        pos = garray_bisect_right(h->successors, &ng);
        garray_insert(h->successors, pos, &ng);
    }
    n = ngram_trie_node_n(t, ng);
    if (n > t->n) {
        E_INFO("Updated N to %d\n", n);
        t->n = n;
    }
    return ng;
}


int
ngram_trie_add_successor_ngram(ngram_trie_t *t,
                               ngram_trie_node_t *h,
                               ngram_trie_node_t *w)
{
    int n;

    assert(w->word >= 0);
    assert(w->log_prob <= 0);
    if (h->successors == NULL) {
        h->successors = garray_init(1, sizeof(ngram_trie_node_t *));
        garray_set_cmp(h->successors, &ngram_trie_nodeptr_cmp, t);
        garray_put(h->successors, 0, &w);
    }
    else {
        size_t pos;
        pos = garray_bisect_right(h->successors, &w);
        garray_insert(h->successors, pos, &w);
    }
    w->history = h;
    w->backoff = (ngram_trie_node_t *)-1;
    n = ngram_trie_node_n(t, w);
    if (n > t->n) {
        E_INFO("Updated N to %d\n", n);
        t->n = n;
    }

    return 0;
}

int
ngram_trie_rename_successor(ngram_trie_t *t,
                            ngram_trie_node_t *h,
                            ngram_trie_node_t *w,
                            int32 new_wid)
{

    size_t pos;

    assert(w->word >= 0);
    assert(new_wid >= 0);
    assert(w->log_prob <= 0);

    pos = garray_bisect_left(h->successors, &w);
    if (pos == garray_next_idx(h->successors)) {
        E_ERROR("Could not find successor to rename\n");
        return -1;
    }
    assert(w == garray_ent(h->successors, ngram_trie_node_t *, pos));
    garray_delete(h->successors, pos, pos + 1);

    w->word = new_wid;
    pos = garray_bisect_right(h->successors, &w);
    garray_insert(h->successors, pos, &w);

    return 0;
}

int
ngram_trie_node_n(ngram_trie_t *t, ngram_trie_node_t *ng)
{
    ngram_trie_node_t *h;
    int n = 0;
    for (h = ng->history; h; h = h->history)
        ++n;
    return n;
}

int32
ngram_trie_node_get_word_hist(ngram_trie_t *t,
                              ngram_trie_node_t *ng,
                              int32 *out_hist)
{
    ngram_trie_node_t *h;
    int32 n_hist;

    if (ng->history == NULL)
        return -1;
    n_hist = 0;
    for (h = ng->history; h->word != -1; h = h->history) {
        if (out_hist)
            out_hist[n_hist] = h->word;
        ++n_hist;
    }
    return n_hist;
}

ngram_trie_node_t *
ngram_trie_backoff(ngram_trie_t *t,
                   ngram_trie_node_t *ng)
{
    ngram_trie_node_t *bong;
    int32 *hist;
    int32 n_hist;

    if (ng->backoff != (ngram_trie_node_t *)-1)
        return ng->backoff;

    /* Extract word IDs from ng's history. */
    n_hist = ngram_trie_node_get_word_hist(t, ng, NULL);
    hist = ckd_calloc(n_hist, sizeof(*hist));
    ngram_trie_node_get_word_hist(t, ng, hist);

    /* Look up the backoff N-Gram. */
    bong = ngram_trie_ngram_v(t, ng->word, hist, n_hist - 1);
    ckd_free(hist);

    ng->backoff = bong;
    return bong;
}

int32
ngram_trie_bowt_v(ngram_trie_t *t, int32 w,
                  int32 const *hist, int32 n_hist)
{
    ngram_trie_node_t *ng;

#if 0
    E_INFO("Getting backoff weight for N-Gram %s |",
           dict_wordstr(t->dict, w));
    int i;
    for (i = 0; i < n_hist; ++i) {
        E_INFOCONT(" %s", dict_wordstr(t->dict, hist[n_hist - 1 - i]));
    }
    E_INFOCONT("\n");
#endif

    if ((ng = ngram_trie_ngram_v(t, w, hist, n_hist)) != NULL)
        return ng->log_bowt << t->shift;
    else
#if 0 /* While this seems like it would be correct, it isn't what the
       * other LM code does. */
        if (n_hist > 0)
        return ngram_trie_bowt_v(t, hist[0], hist + 1, n_hist - 1);
        else
#endif
        return 0;
}

int32
ngram_trie_prob_v(ngram_trie_t *t, int *n_used, int32 w,
                  int32 const *hist, int32 n_hist)
{
    ngram_trie_node_t *ng;

#if 0
    E_INFO("Scoring N-Gram %s |",
           dict_wordstr(t->dict, w));
    int i;
    for (i = 0; i < n_hist; ++i) {
        E_INFOCONT(" %s", dict_wordstr(t->dict, hist[n_hist - 1 - i]));
    }
    E_INFOCONT("\n");
#endif

    if (n_used) *n_used = n_hist + 1;
    if ((ng = ngram_trie_ngram_v(t, w, hist, n_hist)) != NULL) {
        assert(ng->word == w);
        return ng->log_prob << t->shift;
    }
    else if (n_hist > 0) {
        int32 bong, pong;

        if (n_used) --*n_used;
        bong = ngram_trie_prob_v(t, n_used, w, hist, n_hist - 1);
        pong = ngram_trie_bowt_v(t, hist[0], hist + 1, n_hist - 1);
        return bong + pong;
    }
    else {
        if (n_used) *n_used = 0;
        return t->zero << t->shift;
    }
}

int32
ngram_trie_successor_prob(ngram_trie_t *t,
                          ngram_trie_node_t *h, int32 w)
{
    int32 prob;
    int32 *hist;
    int32 n_hist;

    /* Extract word IDs from ng's history. */
    n_hist = ngram_trie_node_get_word_hist(t, h, NULL) + 1;
    hist = ckd_calloc(n_hist, sizeof(*hist));
    ngram_trie_node_get_word_hist(t, h, hist + 1);
    hist[0] = h->word;
    
    prob = ngram_trie_prob_v(t, NULL, w, hist, n_hist);
    ckd_free(hist);
    return prob;
}

int32
ngram_trie_calc_bowt(ngram_trie_t *t, ngram_trie_node_t *h)
{
    ngram_trie_iter_t *itor;
    double nom, dnom;

    assert(t->finish_wid != BAD_S3WID);
    nom = dnom = 1.0;
    for (itor = ngram_trie_successors(t, h); itor;
         itor = ngram_trie_iter_next(itor)) {
        ngram_trie_node_t *ng = ngram_trie_iter_get(itor);
        ngram_trie_node_t *bong;
        int32 log_prob;
        ngram_trie_node_params(t, ng, &log_prob, NULL);
        nom -= logmath_exp(t->lmath, log_prob);
        if ((bong = ngram_trie_backoff(t, ng)) != NULL) {
            ngram_trie_node_params(t, bong, &log_prob, NULL);
            dnom -= logmath_exp(t->lmath, log_prob);
        }
    }

    if (nom == 0) {
        /* Backoff is futile. */
        return t->zero;
    }
    else if (h->word == t->finish_wid) {
        /* Backoff is impossible. */
        return 0;
    }
    else if (nom < 0 || dnom <= 0) {
        /* Oh, something is actually wrong. */
        E_ERROR("Bad backoff weight for ");
        ngram_trie_node_print(t, h, err_get_logfp());
        E_INFOCONT(": %f / %f\n", nom, dnom);
        return -1;
    }

    return logmath_log(t->lmath, nom / dnom);
}

#define MY_EPSILON 0.01
int32
ngram_trie_node_validate(ngram_trie_t *t, ngram_trie_node_t *h)
{
    ngram_trie_iter_t *ui;
    double tprob;

    tprob = 0.0;
    for (ui = ngram_trie_ngrams(t, 1); ui;
         ui = ngram_trie_iter_next(ui)) {
        int32 wid = ngram_trie_node_word(t, ngram_trie_iter_get(ui));
        int32 log_prob = ngram_trie_successor_prob(t, h, wid);
        double prob = logmath_exp(t->lmath, log_prob);
        tprob += prob;
    }
    if (fabs(tprob - 1.0) > MY_EPSILON) {
        E_ERROR("Validation failed, P(.|H) = %f\n", tprob);
        return 1;
    }
    return logmath_log(t->lmath, tprob);
}

static int
skip_arpa_header(lineiter_t *li)
{
    while (li) {
        string_trim(li->buf, STRING_BOTH);
        if (0 == strcmp(li->buf, "\\data\\")) {
            break;
        }
        li = lineiter_next(li);
    }
    if (li == NULL) {
        E_ERROR("Unexpected end of file when reading ARPA format");
        return -1;
    }
    return 0;
}

static int
read_ngram_counts(lineiter_t *li, garray_t *counts)
{
    int one = 1;

    /* Reset and add the number of zerograms (there is one of them) */
    garray_reset(counts);
    garray_append(counts, &one);

    for (;li;li = lineiter_next(li)) {
        string_trim(li->buf, STRING_BOTH);
        if (strlen(li->buf) == 0)
            break;
        if (0 == strncmp(li->buf, "ngram ", 6)) {
            char *n, *c;
            int ni, ci;
            n = li->buf + 6;
            if (n == NULL) {
                E_ERROR("Invalid N-Gram count line when reading ARPA format");
                return -1;
            }
            c = strchr(n, '=');
            if (c == NULL || c[1] == '\0') {
                E_ERROR("Invalid N-Gram count line when reading ARPA format");
                return -1;
            }
            E_INFO("%s\n", li->buf);
            *c++ = '\0';
            ni = atoi(n);
            ci = atoi(c);
            garray_expand_to(counts, ni + 1);
            garray_put(counts, ni, &ci);
        }
    }
    return garray_size(counts) - 1;
}

static ngram_trie_node_t *
add_ngram_line(ngram_trie_t *t, char *buf, int n,
               char **wptr, int32 *wids,
               ngram_trie_node_t **last_history)
{
    int nwords = str2words(buf, NULL, 0);
    double prob, bowt;
    int i;
    char *libuf = ckd_salloc(buf);
    ngram_trie_node_t *node;

    assert(nwords <= n + 2);
    if (nwords < n + 1) {
        E_ERROR("Expected at least %d fields for %d-Gram\n", n+1, n);
        return NULL;
    }
    str2words(buf, wptr, nwords);

    prob = atof_c(wptr[0]);
    if (nwords == n + 2)
        bowt = atof_c(wptr[n + 1]);
    else
        bowt = 0.0;

    wids[0] = dict_wordid(t->dict, wptr[n]);
    /* Add a unigram word to the dictionary if we are generating one. */
    if (wids[0] == BAD_S3WID) {
        if (t->gendict)
            wids[0] = dict_add_word(t->dict, wptr[n], NULL, 0);
        else {
            E_WARN("Unknown unigram %s in ARPA file, skipping\n");
            return NULL;
        }
        assert(wids[0] != BAD_S3WID);
    }
    for (i = 1; i < n; ++i) {
        wids[i] = dict_wordid(t->dict, wptr[n-i]);
        if (wids[i] == BAD_S3WID) {
            E_WARN("Unknown unigram %s in ARPA file, skipping\n");
            return NULL;
        }
    }
#if 0
    E_INFO("Line is %s N-Gram is %s |",
           libuf,
           dict_wordstr(t->dict, wids[0]));
    for (i = 1; i < n; ++i) {
        E_INFOCONT(" %s", dict_wordstr(t->dict, wids[n-i]));
    }
    E_INFOCONT("\n");
#endif
    ckd_free(libuf);

    /* Determine if this N-Gram has the same history as the previous one. */
    if (n == 1) {
        /* Always the same. */
        assert(*last_history == t->root);
    }
    else {
        /* Unfortunately, we can't (well, we shouldn't) assume that
         * the N-Grams are sorted, so we have to look up the entire
         * history. */
        ngram_trie_node_t *h = *last_history;
        for (i = 1; h != NULL && i < n; ++i) {
            if (h->word != wids[i])
                break;
            h = h->history;
        }
        /* Not an exact match, have to get a new one. */
        if (i < n)
            *last_history = ngram_trie_ngram_v(t, wids[1], wids + 2, n - 2);
        if (*last_history == NULL) {
            E_WARN("Unknown history for N-Gram: %s |",
                   dict_wordstr(t->dict, wids[0]));
            for (i = 1; i < n; ++i) {
                E_INFOCONT(" %s", dict_wordstr(t->dict, wids[n-i]));
            }
            E_INFOCONT("\n");
            return NULL;
        }
    }

    /* Fall through and add a successor to last_history. */
    node = ngram_trie_add_successor(t, *last_history, wids[0]);
    node->log_prob = logmath_log10_to_log(t->lmath, prob) >> t->shift;
    node->log_bowt = logmath_log10_to_log(t->lmath, bowt) >> t->shift;
    return node;
}

static int
read_ngrams(ngram_trie_t *t, lineiter_t *li, int n)
{
    /* Pre-allocate these. */
    char **wptr = ckd_calloc(n + 2, sizeof(*wptr));
    int32 *wids = ckd_calloc(n, sizeof(*wids));
    ngram_trie_node_t *last_history;
    int ngcount;

    if (n == 1)
        last_history = t->root; /* always the same for 1-grams */
    else
        last_history = NULL;
    ngcount = 0;

    for (;li;li = lineiter_next(li)) {
        string_trim(li->buf, STRING_BOTH);
        /* Skip blank lines. */
        if (strlen(li->buf) == 0)
            continue;
        /* No more N-Grams to work with. */
        if (0 == strcmp(li->buf, "\\end\\")) {
            E_INFOCONT(" read %d N-Grams\n", ngcount);
            ckd_free(wids);
            ckd_free(wptr);
            return 0;
        }
        /* Look for an N-Gram start marker. */
        if (li->buf[0] == '\\') {
            char *c;
            int nn;

            if (!isdigit(li->buf[1])) {
                E_ERROR("Expected an N-Gram start marker, got %s", li->buf);
                goto error_out;
            }
            for (c = li->buf + 1; *c && isdigit(*c); ++c)
                ;
            if (0 != strcmp(c, "-grams:")) {
                E_ERROR("Expected an N-Gram start marker, got %s", li->buf);
                goto error_out;
            }
            nn = atoi(li->buf + 1);
            if (nn == n+1) {
                ckd_free(wptr);
                ckd_free(wids);
                E_INFOCONT(" read %d N-Grams\n", ngcount);
                if (ngcount != garray_ent(t->counts, int, n)) {
                    E_WARN("Header claims %d %d-Grams, it's wrong\n",
                           garray_ent(t->counts, int, n), n);
                    garray_ent(t->counts, int, n) = ngcount;
                }
                return nn;
            }
            else if (nn == n) {
                E_INFO("%s", li->buf);
                continue;
            }
            else {
                E_ERROR("Expected %d or %d-grams, got %d (%s)\n",
                        n, n+1, nn, li->buf);
                goto error_out;
            }
        }
        /* Now interpret the line as an N-Gram. */
        if (add_ngram_line(t, li->buf, n, wptr,
                           wids, &last_history) == NULL) {
            goto error_out;
        }
        ++ngcount;
    }
    E_ERROR("Expected \\end\\ or an N-Gram marker\n");
error_out:
    ckd_free(wptr);
    ckd_free(wids);
    return -1;
}

int
ngram_trie_read_arpa(ngram_trie_t *t, FILE *arpafile)
{
    lineiter_t *li;
    int n;

    li = lineiter_start(arpafile);

    /* Skip header text. */
    if (skip_arpa_header(li) < 0)
        return -1;

    /* Read N-Gram counts. */
    if ((t->n = read_ngram_counts(li, t->counts)) < 0)
        return -1;

    /* Now read each set of N-Grams for 1..n */
    n = 1;
    while ((n = read_ngrams(t, li, n)) > 1)
        ;
    lineiter_free(li);
    if (n < 0)
        return -1;

    return 0;
}

int ngram_trie_node_print(ngram_trie_t *t,
                          ngram_trie_node_t *ng,
                          FILE *outfh)
{
    int32 *wids;
    int n_wids;

    n_wids = ngram_trie_node_get_word_hist(t, ng, NULL) + 1;
    wids = ckd_calloc(n_wids, sizeof(*wids));
    wids[0] = ng->word;
    n_wids = ngram_trie_node_get_word_hist(t, ng, wids + 1) + 1;
    fprintf(outfh, "%s", dict_wordstr(t->dict, wids[--n_wids]));
    for (--n_wids; n_wids >= 0; --n_wids)
        fprintf(outfh, " %s", dict_wordstr(t->dict, wids[n_wids]));
    ckd_free(wids);
    return 0;
}

int
ngram_trie_update_counts(ngram_trie_t *t)
{
    int n;

    for (n = 1; n <= t->n; ++n) {
        ngram_trie_iter_t *itor;
        size_t n_ngrams = 0;

        for (itor = ngram_trie_ngrams(t, n); itor;
             itor = ngram_trie_iter_next(itor)) {
            ngram_trie_node_t *node = ngram_trie_iter_get(itor);
            /* Detect an increased N-Gram order. */
            if (n == t->n && node->successors != NULL) {
                ++t->n;
                garray_expand_to(t->counts, t->n + 1);
                garray_ent(t->counts, int, t->n) = 0;
            }
            ++n_ngrams;
        }
        /* Detect complete deletion of all N-Grams. */
        if (n_ngrams == 0)
            t->n = n - 1;
        garray_ent(t->counts, int, n) = n_ngrams;
    }
    return 0;
}

int
ngram_trie_write_arpa(ngram_trie_t *t, FILE *arpafile)
{
    int32 *wids;
    int n;

    fprintf(arpafile, "# Written by ngram_trie.c\n");
    fprintf(arpafile, "\\data\\\n");

    ngram_trie_update_counts(t);
    for (n = 1; n <= t->n; ++n)
        fprintf(arpafile, "ngram %d=%d\n", n, garray_ent(t->counts, int, n));

    wids = ckd_calloc(t->n, sizeof(*wids));
    for (n = 1; n <= t->n; ++n) {
        ngram_trie_iter_t *itor;

        fprintf(arpafile, "\n\\%d-grams:\n", n);
        for (itor = ngram_trie_ngrams(t, n); itor;
             itor = ngram_trie_iter_next(itor)) {
            ngram_trie_node_t *ng = ngram_trie_iter_get(itor);
            int n_wids;

            wids[0] = ng->word;
            n_wids = ngram_trie_node_get_word_hist(t, ng, wids + 1) + 1;
            assert(n_wids == n);
            fprintf(arpafile, "%.4f",
                    logmath_log_to_log10(t->lmath, ng->log_prob << t->shift));
            for (--n_wids; n_wids >= 0; --n_wids)
                fprintf(arpafile, " %s", dict_wordstr(t->dict, wids[n_wids]));
            if (ng->log_bowt != 0)
                fprintf(arpafile, " %.4f",
                        logmath_log_to_log10(t->lmath, ng->log_bowt << t->shift));
            fprintf(arpafile, "\n");
        }
    }
    fprintf(arpafile, "\n\\end\\\n");
    ckd_free(wids);
    return 0;
}
