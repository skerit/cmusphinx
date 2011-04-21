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
 * @file ms_lattice.c Word lattices for MultiSphinx.
 */

#include <string.h>

#include <sphinxbase/garray.h>
#include <sphinxbase/bitvec.h>
#include <sphinxbase/gq.h>
#include <sphinxbase/pio.h>
#include <sphinxbase/err.h>
#include <sphinxbase/strfuncs.h>
#include <sphinxbase/ckd_alloc.h>

#include <multisphinx/ngram_trie.h>

#include "ms_lattice.h"

struct ms_lattice_s {
    int refcount;
    /**
     * Log base calculation.
     */
    logmath_t *lmath;
    /**
     * Dictionary for mapping word IDs.
     */
    dict_t *dict;
    /**
     * Do we create new dictionary entries as needed?
     */
    int autodict;
    /**
     * Trie mapping words to language model states.
     */
    ngram_trie_t *lmsids;
    /**
     * List mapping language model states to trie entries.
     */
    garray_t *lms;
    /**
     * List of lattice nodes
     *
     * Lattice nodes are identified by the combination of language model
     * state and start frame - because the sort order may change in the
     * lattice structure we maintain an auxiliary table to preserve this
     * mapping.
     */
    garray_t *node_list;
    /**
     * List of lattice links
     *
     * Since the same link is shared between source and destination
     * node, we maintain a shared pool of links which is referenced
     * indirectly from the node structures.
     */
    garray_t *link_list;
    /**
     * Mapping of lattice node IDs to node list indices.
     */
    nodeid_map_t *node_map;

    int32 start_idx;   /**< Index of start node. */
    int32 end_idx;     /**< Index of end node. */
    int32 next_frame;  /**< One past last frame in lattice (FIXME) */

    /**
     * Total probability of this lattice (normalizer for posteriors)
     */
    int32 norm;
    /**
     * History array used for language model expansion.
     */
    int32 *lmhist;
    /**
     * History array used for language model expansion.
     */
    int32 *lathist;
    /**
     * Size of history arrays.
     */
    int max_n_hist;
};

struct ms_latnode_iter_s {
    ms_lattice_t *l;
    int32 cur;
    /* FIXME: need a type field (shut up, C++ fans) */
    int32 start;     /**< Only for reverse iteration. */
    int32 end;       /**< Only for forward iteration. */
    int32 frame_idx; /**< Only for frame iteration. */
    gq_t *q;
};

ms_lattice_t *
ms_lattice_init(logmath_t *lmath, dict_t *dict)
{
    ms_lattice_t *l = ckd_calloc(1, sizeof(*l));
    l->refcount = 1;
    l->lmath = logmath_retain(lmath);
    if (dict) {
        l->dict = dict_retain(dict);
        l->autodict = FALSE;
    }
    else {
        l->dict = dict_init(NULL, NULL);
        l->autodict = TRUE;
    }
    l->node_list = garray_init(0, sizeof(ms_latnode_t));
    l->link_list = garray_init(0, sizeof(ms_latlink_t));
    l->node_map = nodeid_map_init();
    l->start_idx = l->end_idx = -1;
    l->norm = logmath_get_zero(l->lmath);
    l->lmsids = ngram_trie_init(l->dict, l->lmath);
    /* Create an initial dummy index to catch uninitialized ngram trie
     * nodes. */
    l->lms = garray_init(1, sizeof(ngram_trie_node_t *));
    return l;
}

ms_lattice_t *
ms_lattice_retain(ms_lattice_t *l)
{
    ++l->refcount;
    return l;
}

int
ms_lattice_free(ms_lattice_t *l)
{
    int i;
    if (l == NULL)
        return 0;
    if (--l->refcount > 0)
        return l->refcount;

    ckd_free(l->lmhist);
    ckd_free(l->lathist);
    for (i = 0; i < garray_size(l->node_list); ++i) {
        ms_latnode_t *node = ms_lattice_get_node_idx(l, i);
        garray_free(node->entries);
        garray_free(node->exits);
    }
    garray_free(l->node_list);
    garray_free(l->link_list);
    nodeid_map_free(l->node_map);
    garray_free(l->lms);
    ngram_trie_free(l->lmsids);
    logmath_free(l->lmath);
    dict_free(l->dict);
    ckd_free(l);
    return 0;
}

dict_t *
ms_lattice_dict(ms_lattice_t *l)
{
    return l->dict;
}

logmath_t *
ms_lattice_lmath(ms_lattice_t *l)
{
    return l->lmath;
}

int32
ms_lattice_lmstate_init(ms_lattice_t *l, int32 w,
                        int32 const *hist, int32 n_hist)
{
    ngram_trie_node_t *ng;
    int32 lmstate = garray_size(l->lms);
    if (lmstate == 0xdeadbeef) { /* unlikely to happen */
        garray_expand(l->lms, lmstate + 1);
        lmstate = garray_size(l->lms); /* 0xdeadbef0 */
    }
    ng = ngram_trie_ngram_init_v(l->lmsids, w, hist, n_hist);
    ngram_trie_node_set_params_raw(l->lmsids, ng,
                                   lmstate >> 16,
                                   lmstate & 0xffff);
    garray_append(l->lms, &ng);
    return lmstate;
}

int32
ms_lattice_get_lmstate_idx(ms_lattice_t *l, int32 w,
                           int32 const *hist, int32 n_hist)
{
    ngram_trie_node_t *ng;
    int16 ls_hi, ls_lo;

    ng = ngram_trie_ngram_v(l->lmsids, w, hist, n_hist);
    if (ng == NULL)
        return -1;
    ngram_trie_node_params_raw(l->lmsids, ng,
                               &ls_hi, &ls_lo);
    /* This might be an intermediate node created by
     * ngram_trie_ngram_init_v(), in which case its index will be
     * zero - hence the dummy index in ms_lattice_init(). */
    if (ls_hi == 0 && ls_lo == 0) {
        int32 lmstate = garray_size(l->lms);
        ngram_trie_node_set_params_raw(l->lmsids, ng,
                                       lmstate >> 16,
                                       lmstate & 0xffff);
        garray_append(l->lms, &ng);
        return lmstate;
    }
    else
        return ((int32)ls_hi << 16) | ((int32)ls_lo & 0xffff);
}

int
ms_lattice_get_lmstate_wids(ms_lattice_t *l, int32 idx,
                            int32 *out_w, int32 *out_hist)
{
    ngram_trie_node_t *ng;

    if (idx >= garray_size(l->lms) || idx < 0) {
        /* An index of -1 should result in a word ID of -1.  This is
         * used for epsilon states (FIXME: arguably, abused). */
        if (idx == -1) *out_w = -1;
        return 0;
    }
    ng = garray_ent(l->lms, ngram_trie_node_t *, idx);
    if (out_w) *out_w = ngram_trie_node_word(l->lmsids, ng);
    return ngram_trie_node_get_word_hist(l->lmsids, ng, out_hist);
}

ms_latnode_t *
ms_lattice_node_init(ms_lattice_t *l, int sf, int32 lmstate)
{
    ms_latnode_t *node;
    int32 nodeidx;

    nodeidx = garray_next_idx(l->node_list);
    garray_expand(l->node_list, nodeidx + 1);
    node = garray_ptr(l->node_list, ms_latnode_t, nodeidx);
    nodeid_map_add(l->node_map, sf, lmstate, nodeidx);

    node->exits = node->entries = NULL;
    node->fan = 0;
    node->id.sf = sf;
    node->id.lmstate = lmstate;

    if (sf > l->next_frame)
        l->next_frame = sf;

    return node;
}

int
ms_lattice_get_n_frames(ms_lattice_t *l)
{
    return l->next_frame;
}

int32
ms_lattice_get_idx_node(ms_lattice_t *l, ms_latnode_t *node)
{
    return garray_idx(l->node_list, node);
}

int32
ms_lattice_get_idx_link(ms_lattice_t *l, ms_latlink_t *link)
{
    return garray_idx(l->link_list, link);
}

ms_latlink_t *
ms_lattice_get_link_idx(ms_lattice_t *l, int32 idx)
{
    return garray_ptr(l->link_list, ms_latlink_t, idx);
}

ms_latnode_t *
ms_lattice_get_node_idx(ms_lattice_t *l, int32 idx)
{
    return garray_ptr(l->node_list, ms_latnode_t, idx);
}

ms_latnode_t *
ms_lattice_get_node_id(ms_lattice_t *l, int sf, int32 lmstate)
{
    int32 idx = nodeid_map_map(l->node_map, sf, lmstate);
    if (idx == -1)
        return NULL;
    return ms_lattice_get_node_idx(l, idx);
}

ms_latnode_t *
ms_lattice_set_node_id(ms_lattice_t *l, ms_latnode_t *node,
                       int sf, int32 lmstate)
{
    int32 idx = ms_lattice_get_idx_node(l, node);
    if (nodeid_map_delete(l->node_map, node->id.sf, node->id.lmstate) < 0)
        return NULL;
    if (nodeid_map_add(l->node_map, sf, lmstate, idx) < 0)
        return NULL;
    node->id.sf = sf;
    node->id.lmstate = lmstate;
    return node;
}

int32
ms_lattice_set_start(ms_lattice_t *l, ms_latnode_t *node)
{
    int32 idx = ms_lattice_get_idx_node(l, node);
    l->start_idx = idx;
    return idx;
}

int32
ms_lattice_set_end(ms_lattice_t *l, ms_latnode_t *node)
{
    int32 idx = ms_lattice_get_idx_node(l, node);
    l->end_idx = idx;
    return idx;
}

ms_latnode_t *
ms_lattice_get_start(ms_lattice_t *l)
{
    return ms_lattice_get_node_idx(l, l->start_idx);
}

ms_latnode_t *
ms_lattice_get_end(ms_lattice_t *l)
{
    return ms_lattice_get_node_idx(l, l->end_idx);
}

ms_latlink_t *
ms_lattice_link(ms_lattice_t *l,
                ms_latnode_t *src, ms_latnode_t *dest,
                int32 wid, int32 ascr)
{
    ms_latlink_t *arc;
    int32 lid;
    int32 zero = logmath_get_zero(l->lmath);

    if (src->exits == NULL)
        src->exits = garray_init(0, sizeof(int32));
    if (dest->entries == NULL)
        dest->entries = garray_init(0, sizeof(int32));

    lid = garray_next_idx(l->link_list);
    garray_expand(l->link_list, lid + 1);
    arc = garray_ptr(l->link_list, ms_latlink_t, lid);
    garray_append(src->exits, &lid);
    garray_append(dest->entries, &lid);

    arc->wid = wid;
    arc->ascr = ascr;
    arc->lscr = zero;
    arc->alpha = zero;
    arc->beta = zero;

    arc->src = garray_idx(l->node_list, src);
    assert(arc->src != GARRAY_INVALID_INDEX);
    assert(src == ms_lattice_get_node_idx(l, arc->src));
    arc->dest = garray_idx(l->node_list, dest);
    assert(arc->dest != GARRAY_INVALID_INDEX);
    assert(dest == ms_lattice_get_node_idx(l, arc->dest));

    return arc;
}

static int
get_or_create_wid(ms_lattice_t *l, char *word, int alt)
{
    int32 wid;

    /* HACK, but what can we do... */
    if (0 == strcmp(word, "!SENT_END"))
        strcpy(word, "</s>");
    if (0 == strcmp(word, "!SENT_START"))
        strcpy(word, "<s>");
    /* NOTE: In reality we want to collapse alternate pronunciations,
     * but at this point it's best to be faithful to the input.
     * Therefore we include the alternate pronunciation in the
     * language model state, even though that makes no sense in the
     * grand scheme of things. */
    if (alt != 1) {
        /* First make sure the base word exists. */
        wid = dict_wordid(l->dict, word);
        if (wid == -1 && l->autodict)
            wid = dict_add_word(l->dict, word, NULL, 0);
        if (wid == -1) {
            E_ERROR("No base word for %s(%d) in dictionary\n", word, alt);
            return -1;
        }
        sprintf(word + strlen(word), "(%d)", alt);
    }
    wid = dict_wordid(l->dict, word);
    if (wid == -1 && l->autodict)
        wid = dict_add_word(l->dict, word, NULL, 0);
    return wid;
}

static int
process_htk_node_line(ms_lattice_t *l, lineiter_t *li,
                      garray_t *wptr, int n_wptr, int frate)
{
    ms_latnode_t *node = NULL;
    char *word = NULL;
    int32 nodeidx, wid;
    int i, sf, alt = 1;

    for (i = 0; i < n_wptr; ++i) {
        char *f = garray_ent(wptr, char *, i);
        char *e = strchr(f, '=');
        if (e == NULL) {
            E_ERROR("Invalid field %s in line %d\n",
                    f, (int)li->lineno);
            ckd_free(word);
            return -1;
        }
        *e++ = '\0';
        if (0 == strcmp(f, "I")) {
            nodeidx = atoi(e);
            garray_expand(l->node_list, nodeidx + 1);
            node = garray_ptr(l->node_list, ms_latnode_t, nodeidx);
        }
        else if (0 == strcmp(f, "t")) {
            sf = (int)(atof_c(e) * frate);
        }
        else if (0 == strcmp(f, "W")) {
            word = ckd_calloc(1, strlen(e) + 16);
            strcpy(word, e);
        }
        else if (0 == strcmp(f, "v")) {
            alt = atoi(e);
            if (alt < 1 || alt > 255) {
                E_ERROR("Invalid pronunciation alternate %s at line %d\n",
                        e, (int)li->lineno);
                ckd_free(word);
                return -1;
            }
        }
        else {
            E_WARN("Unknown field type %s in line %d\n",
                   f, (int)li->lineno);
        }
    }
    if (node == NULL) {
        E_ERROR("Found no node ID in line %d\n", (int)li->lineno);
        ckd_free(word);
        return -1;
    }
    wid = get_or_create_wid(l, word, alt);
    if (wid == -1) {
        E_ERROR("Failed to add word %s to dictionary\n", word);
        ckd_free(word);
        return -1;
    }
    node->id.sf = sf;
    /* Enter this word ID as a language model state. */
    node->id.lmstate = ms_lattice_lmstate_init(l, wid, NULL, 0);
    node->fan = 0;
    node->exits = node->entries = NULL;
    nodeid_map_add(l->node_map, sf, wid, nodeidx);
    if (sf > l->next_frame)
        l->next_frame = sf;

    ckd_free(word);
    return 0;
}

static int
process_htk_arc_line(ms_lattice_t *l, lineiter_t *li, garray_t *wptr, int n_wptr)
{
    char *word = NULL;
    ms_latnode_t *src = NULL, *dest = NULL;
    ms_latlink_t *arc = NULL;
    int32 wid, ascr, prob;
    int i, alt;

    for (i = 0; i < n_wptr; ++i) {
        char *f = garray_ent(wptr, char *, i);
        char *e = strchr(f, '=');
        if (e == NULL) {
            E_ERROR("Invalid field %s in line %d\n",
                    f, (int)li->lineno);
            ckd_free(word);
            return -1;
        }
        *e++ = '\0';
        if (0 == strcmp(f, "J")) {
            /* Link ID is irrelevant (for now at least). */
        }
        else if (0 == strcmp(f, "S")) {
            src = ms_lattice_get_node_idx(l, atoi(e));
        }
        else if (0 == strcmp(f, "E")) {
            dest = ms_lattice_get_node_idx(l, atoi(e));
        }
        else if (0 == strcmp(f, "W")) {
            word = ckd_calloc(1, strlen(e) + 16);
            strcpy(word, e);
        }
        else if (0 == strcmp(f, "v")) {
            alt = atoi(e);
            if (alt < 1 || alt > 255) {
                E_ERROR("Invalid pronunciation alternate %s at line %d\n",
                        e, (int)li->lineno);
                ckd_free(word);
                return -1;
            }
        }
        else if (0 == strcmp(f, "a")) {
            ascr = logmath_ln_to_log(l->lmath, atof_c(e));
        }
        else if (0 == strcmp(f, "p")) {
            prob = logmath_log(l->lmath, atof_c(e));
        }
        else {
            E_WARN("Unknown field type %s in line %d\n",
                   f, (int)li->lineno);
        }
    }
    if (src == NULL || dest == NULL) {
        E_ERROR("Found no valid src and dest IDs in line %d\n",
                (int)li->lineno);
        ckd_free(word);
        return -1;
    }
    if (word)
        wid = get_or_create_wid(l, word, alt);
    else {
        /* Get the word ID corresponding to the "lmstate" of the
         * source node (not its actual LM state but we will fix that
         * in expansion) */
        ms_lattice_get_lmstate_wids(l, src->id.lmstate,
                                    &wid, NULL);
    }
    arc = ms_lattice_link(l, src, dest, wid, ascr);

    ckd_free(word);
    return 0;
}

int
ms_lattice_read_htk(ms_lattice_t *l, FILE *fh, int frate)
{
    lineiter_t *li;
    garray_t *wptr;
    int n_wptr;
    int n_nodes = 0, n_arcs = 0;
    int start_idx = -1, end_idx = -1;
    enum {
        HEADER, ENTRIES
    } state = HEADER;

    wptr = garray_init(0, sizeof(char *));
    for (li = lineiter_start(fh); li; li = lineiter_next(li)) {
        string_trim(li->buf, STRING_BOTH);
        if (li->buf[0] == '#')
            continue;
        switch (state) {
        case HEADER:
        {
            char *e;
            if (0 == strncmp(li->buf, "N=", 2)) {
                n_nodes = atoi(li->buf + 2);
                if ((e = strchr(li->buf + 2, '=')) == NULL) {
                    E_ERROR("Invalid node/link count line %d: %s\n",
                            (int)li->lineno, li->buf);
                    goto error_out;
                }
                n_arcs = atoi(e + 1);
                state = ENTRIES;
            }
            else {
                if ((e = strchr(li->buf, '=')) == NULL) {
                    E_ERROR("Invalid header count line %d: %s\n",
                            (int)li->lineno, li->buf);
                    goto error_out;
                }
                *e++ = '\0';
                if (0 == strcmp(li->buf, "start"))
                    start_idx = atoi(e);
                else if (0 == strcmp(li->buf, "end"))
                    end_idx = atoi(e);
            }
            break;
        }
        case ENTRIES: {
            int i;

            n_wptr = str2words(li->buf, NULL, 0);
            garray_expand(wptr, n_wptr);
            str2words(li->buf, garray_void(wptr, 0), n_wptr);
            for (i = 0; i < n_wptr; ++i) {
                char *f = garray_ent(wptr, char *, i);
                if (0 == strncmp(f, "I=", 2)) {
                    if (process_htk_node_line(l, li, wptr, n_wptr, frate) < 0)
                        goto error_out;
                    break;
                }
                if (0 == strncmp(f, "J=", 2)) {
                    if (process_htk_arc_line(l, li, wptr, n_wptr) < 0)
                        goto error_out;
                    break;
                }
            }
            if (i == n_wptr)
                E_WARN("Not a node or arc on line %d\n", (int)li->lineno);
            break;
        }
        }
    }
    if (start_idx == -1) {
        E_WARN("No explicit start node, using first node\n");
        start_idx = 0;
    }
    if (end_idx == -1) {
        end_idx = garray_size(l->node_list) - 1;
        E_WARN("No explicit end node, using last node %d\n", end_idx);
    }
    l->start_idx = start_idx;
    l->end_idx = end_idx;
    garray_free(wptr);
    return 0;

error_out:
    lineiter_free(li);
    garray_free(wptr);
    return -1;
}

int
ms_lattice_write_htk(ms_lattice_t *l, FILE *fh, int frate)
{
    int32 zero = logmath_get_zero(l->lmath);
    int i;

    /* Write a generic header. */
    fprintf(fh, "# Lattice generated by MultiSphinx\n");
    fprintf(fh, "#\n");
    fprintf(fh, "# Header\n");
    fprintf(fh, "#\n");
    fprintf(fh, "VERSION=1.0\n");
    fprintf(fh, "start=%d\n", l->start_idx);
    fprintf(fh, "end=%d\n", l->end_idx);
    fprintf(fh, "#\n");
    fprintf(fh, "N=%u\tL=%u\n",
            (unsigned)garray_size(l->node_list),
            (unsigned)garray_size(l->link_list));
    fprintf(fh, "#\n");
    fprintf(fh, "# Node definitions\n");
    fprintf(fh, "#\n");
    for (i = 0; i < garray_size(l->node_list); ++i) {
        ms_latnode_t *node = garray_ptr(l->node_list, ms_latnode_t, i);
        if (node->id.lmstate != -1) {
            char const *basestr;
            int32 wid;
            ms_lattice_get_lmstate_wids(l, node->id.lmstate, &wid, NULL);
            basestr = dict_basestr(l->dict, wid);
            if (node->id.lmstate == dict_startwid(l->dict))
                basestr = "!SENT_START";
            if (node->id.lmstate == dict_finishwid(l->dict))
                basestr = "!SENT_END";
            fprintf(fh, "I=%d\tt=%.2f\tW=%s\tv=%d\n",
                    i, (double)node->id.sf / frate,
                    basestr,
                    dict_altid(l->dict, node->id.lmstate));
        }
        else
            fprintf(fh, "I=%d\tt=%.2f\n",
                    i, (double)node->id.sf / frate);
    }
    fprintf(fh, "#\n");
    fprintf(fh, "# Link definitions\n");
    fprintf(fh, "#\n");
    for (i = 0; i < garray_size(l->link_list); ++i) {
        ms_latlink_t *link = garray_ptr(l->link_list, ms_latlink_t, i);
        fprintf(fh, "J=%d\tS=%d\tE=%d\ta=%f",
                i, link->src, link->dest,
                logmath_log_to_ln(l->lmath, link->ascr));
        if (link->wid != -1)
            fprintf(fh, "\tW=%s",
                    dict_basestr(l->dict, link->wid));
        if (link->lscr != zero)
            fprintf(fh, "\tl=%f",
                    logmath_log_to_ln(l->lmath, link->lscr));
        if (link->alpha != zero
            && link->beta != zero
            && l->norm != zero)
            fprintf(fh, "\tp=%g",
                    logmath_exp(l->lmath, link->alpha + link->beta - l->norm));
        fprintf(fh, "\n");
    }

    return 0;
}

static int
print_dot_nodeid(ms_lattice_t *l, ms_latnode_t *node, FILE *fh)
{
    if (node->id.lmstate == 0xdeadbeef)
        return -1;
    if (node->id.lmstate != -1) {
        int32 wid;
        ms_lattice_get_lmstate_wids(l, node->id.lmstate, &wid, NULL);
        return fprintf(fh, " \"%s/%d\"",
                       dict_wordstr(l->dict, wid),
                       node->id.sf);
    }
    else
        return fprintf(fh, " \"&epsilon;/%d\"", node->id.sf);
}

int
ms_lattice_write_dot(ms_lattice_t *l, FILE *fh)
{
    int32 zero = logmath_get_zero(l->lmath);
    ms_latnode_t *node;
    int i;

    fprintf(fh, "digraph lattice {\n\trankdir=LR;\n\t");
    fprintf(fh, "\tnode [shape=circle];");
    for (i = 0; i < garray_size(l->node_list); ++i) {
        node = garray_ptr(l->node_list, ms_latnode_t, i);
        if (i != l->end_idx)
            print_dot_nodeid(l, node, fh);
    }
    fprintf(fh, "\n");
    node = ms_lattice_get_end(l);
    fprintf(fh, "\tnode [shape=doublecircle];");
    print_dot_nodeid(l, node, fh);
    fprintf(fh, "\n\n");
    for (i = 0; i < garray_size(l->node_list); ++i) {
        int j;
        node = garray_ptr(l->node_list, ms_latnode_t, i);
        if (node->exits == NULL)
            continue;
        for (j = 0; j < ms_latnode_n_exits(node); ++j) {
            ms_latlink_t *link = ms_latnode_get_exit(l, node, j);
            ms_latnode_t *node2, *node3;
            double weight;

            /* FIXME: Ad hoc behaviour for weights, should be configurable. */
            if (link->alpha != zero
                && link->beta != zero
                && l->norm != zero)
                weight = logmath_exp(l->lmath, link->alpha + link->beta - l->norm);
            else if (link->lscr != zero)
                weight = logmath_log_to_ln(l->lmath, link->lscr);
            else
                weight = logmath_log_to_ln(l->lmath, link->ascr);

            node2 = ms_lattice_get_node_idx(l, link->src);
            if (node2->id.lmstate == 0xdeadbeef)
                continue;
            node3 = ms_lattice_get_node_idx(l, link->dest);
            if (node3->id.lmstate == 0xdeadbeef)
                continue;
            fprintf(fh, "\t");
            print_dot_nodeid(l, node2, fh);
            fprintf(fh, " ->");
            print_dot_nodeid(l, node3, fh);

            if (link->wid != -1)
                fprintf(fh, " [label=\"%s/%.2g\"];\n",
                        dict_wordstr(l->dict, link->wid), weight);
            else
                fprintf(fh, " [label=\"%.2g\"];\n", weight);
        }
    }
    fprintf(fh, "}\n");
    return 0;
}

ms_latnode_iter_t *
ms_lattice_traverse_topo(ms_lattice_t *l,
                         ms_latnode_t *end)
{
    ms_latnode_iter_t *itor;
    int i;

    for (i = 0; i < garray_size(l->node_list); ++i) {
        ms_latnode_t *node = garray_ptr(l->node_list, ms_latnode_t, i);
        node->fan = 0;
    }
    for (i = 0; i < garray_size(l->node_list); ++i) {
        ms_latnode_t *node = garray_ptr(l->node_list, ms_latnode_t, i);
        int j;
        if (node->exits == NULL)
            continue;
        for (j = 0; j < ms_latnode_n_exits(node); ++j) {
            ms_latlink_t *link = ms_latnode_get_exit(l, node, j);
            ms_latnode_t *node2 = ms_lattice_get_node_idx(l, link->dest);
            ++node2->fan;
        }
    }

    itor = ckd_calloc(1, sizeof(*itor));
    itor->l = l;
    itor->cur = ms_lattice_get_idx_node(l, ms_lattice_get_start(l));
    itor->start = -1;
    itor->frame_idx = -1;
    if (end == NULL)
        itor->end = ms_lattice_get_idx_node(l, ms_lattice_get_end(l));
    else
        itor->end = ms_lattice_get_idx_node(l, end);
    itor->q = gq_init(sizeof(int32));

    return itor;
}

ms_latnode_iter_t *
ms_lattice_reverse_topo(ms_lattice_t *l,
                        ms_latnode_t *start)
{
    ms_latnode_iter_t *itor;
    int i;

    for (i = 0; i < garray_size(l->node_list); ++i) {
        ms_latnode_t *node = garray_ptr(l->node_list, ms_latnode_t, i);
        node->fan = 0;
    }
    for (i = 0; i < garray_size(l->node_list); ++i) {
        ms_latnode_t *node = garray_ptr(l->node_list, ms_latnode_t, i);
        int j;
        if (node->exits == NULL)
            continue;
        for (j = 0; j < ms_latnode_n_exits(node); ++j) {
            ms_latlink_t *link = ms_latnode_get_exit(l, node, j);
            ms_latnode_t *node2 = ms_lattice_get_node_idx(l, link->src);
            ++node2->fan;
        }
    }

    itor = ckd_calloc(1, sizeof(*itor));
    itor->l = l;
    if (start == NULL)
        itor->start = ms_lattice_get_idx_node(l, ms_lattice_get_start(l));
    else
        itor->start = ms_lattice_get_idx_node(l, start);
    itor->end = -1;
    itor->frame_idx = -1;
    itor->cur = ms_lattice_get_idx_node(l, ms_lattice_get_end(l));
    itor->q = gq_init(sizeof(int32));

    return itor;
}

ms_latnode_iter_t *
ms_lattice_traverse_frame(ms_lattice_t *l, int frame_idx)
{
    ms_latnode_iter_t *itor;
    int start;

    /* Find the first frame starting in a given node (FIXME: might
     * want to keep track of this in an auxiliary array if this gets
     * slow) */
    for (start = 0; start < garray_size(l->node_list); ++start) {
        ms_latnode_t *node = garray_ptr(l->node_list,
                                        ms_latnode_t, start);
        if (node->id.sf == frame_idx)
            break;
    }
    if (start == garray_size(l->node_list))
        return NULL;

    itor = ckd_calloc(1, sizeof(*itor));
    itor->l = l;
    itor->start = -1;
    itor->end = garray_next_idx(l->node_list);
    itor->frame_idx = frame_idx;
    itor->cur = start;
    itor->q = NULL;

    return itor;
}

ms_latnode_iter_t *
ms_latnode_iter_next(ms_latnode_iter_t *itor)
{
    ms_latnode_t *node;
    garray_t *links;
    int i;

    if (itor->cur == -1) {
        ms_latnode_iter_free(itor);
        return NULL;
    }
    if (itor->cur == itor->start || itor->cur == itor->end) {
        ms_latnode_iter_free(itor);
        return NULL;
    }

    /* Push a frame-based iterator forward/. */
    if (itor->frame_idx != -1) {
        while (++itor->cur != itor->end) {
            node = garray_ptr(itor->l->node_list, ms_latnode_t, itor->cur);
            if (node->id.sf == itor->frame_idx)
                return itor;
        }
        ms_latnode_iter_free(itor);
        return NULL;
    }
    
    /* Figure out which direction we're going. */
    node = ms_lattice_get_node_idx(itor->l, itor->cur);
    if (itor->start != -1)
        links = node->entries;
    else
        links = node->exits;
    if (links == NULL) {
        ms_latnode_iter_free(itor);
        return NULL;
    }

    for (i = 0; i < garray_size(links); ++i) {
        int32 linkid = garray_ent(links, int32, i);
        ms_latlink_t *link = garray_ptr(itor->l->link_list,
                                        ms_latlink_t, linkid);
        ms_latnode_t *next;
        int32 nextid;

        if (itor->start != -1) {
            next = garray_ptr(itor->l->node_list,
                              ms_latnode_t,
                              link->src);
            nextid = link->src;
        }
        else {
            next = garray_ptr(itor->l->node_list,
                              ms_latnode_t,
                              link->dest);
            nextid = link->dest;
        }
        if (--next->fan == 0)
            gq_append(itor->q, &nextid);
    }
    if (gq_size(itor->q) == 0) {
        ms_latnode_iter_free(itor);
        return NULL;
    }
    itor->cur = gq_head(itor->q, int32);
    gq_shift(itor->q, 1);
    return itor;
}

ms_latnode_t *
ms_latnode_iter_get(ms_latnode_iter_t *itor)
{
    return ms_lattice_get_node_idx(itor->l, itor->cur);
}

int32
ms_latnode_iter_get_idx(ms_latnode_iter_t *itor)
{
    return itor->cur;
}

void
ms_latnode_iter_free(ms_latnode_iter_t *itor)
{
    gq_free(itor->q);
    ckd_free(itor);
}

int
ms_latnode_n_entries(ms_latnode_t *node)
{
    return node->entries ? garray_size(node->entries) : 0;
}

int
ms_latnode_n_exits(ms_latnode_t *node)
{
    return node->exits ? garray_size(node->exits) : 0;
}

ms_latlink_t *
ms_latnode_get_entry(ms_lattice_t *l, ms_latnode_t *node, int idx)
{
    int32 linkid = garray_ent(node->entries, int32, idx);
    return garray_ptr(l->link_list, ms_latlink_t, linkid);
}

int32
ms_latnode_get_entry_idx(ms_lattice_t *l,
                         ms_latnode_t *node, int idx)
{
    return garray_ent(node->entries, int32, idx);
}

ms_latlink_t *
ms_latnode_get_exit(ms_lattice_t *l, ms_latnode_t *node, int idx)
{
    int32 linkid = garray_ent(node->exits, int32, idx);
    return garray_ptr(l->link_list, ms_latlink_t, linkid);
}

/**
 * Remove a link from an array (could be implemented generically in
 * garray.c using garray_cmp_int32 and so on, but this is just
 * simpler)
 */
static int
link_array_remove(garray_t *gar, int32 linkid)
{
    int i, j;
    if (gar == NULL)
        return 0;
    for (i = j = 0; i < garray_size(gar); ++i) {
        if (garray_ent(gar, int32, i) == linkid) {
            garray_delete(gar, i, i+1);
            ++j;
        }
    }
    return j;
}

void
ms_latnode_unlink(ms_lattice_t *l, ms_latnode_t *node)
{
    int i;
    for (i = 0; i < ms_latnode_n_entries(node); ++i) {
        ms_latlink_t *link = ms_latnode_get_entry(l, node, i);
        ms_latnode_t *src = ms_lattice_get_node_idx(l, link->src);
        int32 linkid = ms_lattice_get_idx_link(l, link);
        link_array_remove(src->exits, linkid);
    }
    for (i = 0; i < ms_latnode_n_exits(node); ++i) {
        ms_latlink_t *link = ms_latnode_get_exit(l, node, i);
        ms_latnode_t *dest = ms_lattice_get_node_idx(l, link->dest);
        int32 linkid = ms_lattice_get_idx_link(l, link);
        link_array_remove(dest->entries, linkid);
    }
    garray_free(node->exits);
    garray_free(node->entries);
    node->exits = node->entries = NULL;
}

void
ms_latlink_unlink(ms_lattice_t *l, ms_latlink_t *link)
{
    ms_latnode_t *src, *dest;
    int32 linkid = ms_lattice_get_idx_link(l, link);

    src = ms_lattice_get_node_idx(l, link->src);
    dest = ms_lattice_get_node_idx(l, link->dest);
    link_array_remove(src->exits, linkid);
    link_array_remove(dest->entries, linkid);
}

/**
 * Remove unreachable nodes.
 *
 * This just unlinks them, it does not garbage collect.
 */
void
ms_lattice_unlink_unreachable(ms_lattice_t *l)
{
    bitvec_t *forward_reachable, *backward_reachable;
    ms_latnode_iter_t *itor;
    int i;

    forward_reachable = bitvec_alloc(garray_size(l->node_list));
    backward_reachable = bitvec_alloc(garray_size(l->node_list));
    bitvec_set(forward_reachable, l->start_idx);
    bitvec_set(backward_reachable, l->start_idx);
    bitvec_set(forward_reachable, l->end_idx);
    bitvec_set(backward_reachable, l->end_idx);
    for (itor = ms_lattice_traverse_topo(l, NULL);
         itor; itor = ms_latnode_iter_next(itor)) {
        ms_latnode_t *node = ms_latnode_iter_get(itor);
        int32 nodeid = ms_lattice_get_idx_node(l, node);
        assert(nodeid < garray_size(l->node_list));
        bitvec_set(forward_reachable, nodeid);
    }
    for (itor = ms_lattice_reverse_topo(l, NULL);
         itor; itor = ms_latnode_iter_next(itor)) {
        ms_latnode_t *node = ms_latnode_iter_get(itor);
        int32 nodeid = ms_lattice_get_idx_node(l, node);
        assert(nodeid < garray_size(l->node_list));
        bitvec_set(backward_reachable, nodeid);
    }
    for (i = 0; i < garray_size(l->node_list); ++i) {
        ms_latnode_t *node = garray_ptr(l->node_list, ms_latnode_t, i);
        if (node->id.lmstate == 0xdeadbeef
            || !(bitvec_is_set(forward_reachable, i)
                 && bitvec_is_set(backward_reachable, i))) {
            ms_latnode_unlink(l, node);
        }
    }
    bitvec_free(forward_reachable);
    bitvec_free(backward_reachable);
}



/**
 * Rotate a head word into a language model state.
 */
int
rotate_lmstate(int32 wid, int32 *hist, int n_hist, int max_n_hist)
{
    if (wid == -1)
        return 0;
    if (n_hist > 0)
        memmove(hist + 1, hist, (max_n_hist - 1) * sizeof(*hist));
    hist[0] = wid;
    if (n_hist == max_n_hist)
        return max_n_hist;
    else
        return n_hist + 1;
}

/**
 * Copy outgoing arcs from @a backoff to @a dest.
 */
static void
merge_exits(ms_lattice_t *l, ms_latnode_t *new,
            ms_latnode_t *old, int32 lscr)
{
    int y;
    int32 newwid;
    ms_lattice_get_lmstate_wids(l, new->id.lmstate, &newwid, NULL);

    /* In theory we should merge all arcs which lead to the same
     * language model state, but this causes problems because it
     * causes nodes with alternate pronunciations (or rather, their
     * exit arcs) to be prematurely pruned. */
    for (y = 0; y < ms_latnode_n_exits(old); ++y) {
        ms_latlink_t *ylink = ms_latnode_get_exit(l, old, y);
        ms_latnode_t *ysrc = ms_lattice_get_node_idx(l, ylink->src);
        ms_latnode_t *ydest = ms_lattice_get_node_idx(l, ylink->dest);
        ms_latlink_t *yylink;
        int z, duplicate;

        duplicate = FALSE;
        for (z = 0; z < ms_latnode_n_exits(new); ++z) {
            ms_latlink_t *zlink = ms_latnode_get_exit(l, new, z);
            ms_latnode_t *zsrc = ms_lattice_get_node_idx(l, zlink->src);
            if (dict_basewid(l->dict, ylink->wid)
                == dict_basewid(l->dict, zlink->wid)
                /* See above...  In theory this should just check the
                 * base word IDs. */
                && ylink->dest == zlink->dest
                && ysrc->id.sf == zsrc->id.sf) {
                /* Tropical semiring assumption... */
                if (ylink->ascr > zlink->ascr)
                    zlink->ascr = ylink->ascr;
                duplicate = TRUE;
            }
        }
        if (duplicate)
            continue;
        yylink = ms_lattice_link(l, new, ydest, ylink->wid, ylink->ascr);
        yylink->lscr = lscr;
        E_DEBUG(2,("Created new link %s/%d -> %s/%d / %s\n",
                   newwid == -1 ? "&epsilon;" : dict_basestr(l->dict, newwid),
                   new->id.sf,
                   ydestwid == -1 ? "&epsilon;" : dict_basestr(l->dict, ydestwid),
                   ydest->id.sf, dict_basestr(l->dict, ylink->wid)));
    }
    E_DEBUG(2,("%p %s/%d now has %d exits\n",
               new, 
               newwid == -1 ? "&epsilon;" : dict_wordstr(l->dict, newwid),
               new->id.sf, ms_latnode_n_exits(new)));
}

/**
 * Map dictionary wids to language model wids.
 */
static int32
map_lmwid(dict_t *dict, ngram_model_t *lm, int32 wid)
{
    if (wid == -1)
        return wid;
    return ngram_wid(lm, dict_basestr(dict, wid));
}

static void
update_backoff_arcs(ms_lattice_t *l, ms_latnode_t *node,
                    ngram_model_t *lm, int32 endid)
{
    if (node == ms_lattice_get_start(l)
        || ms_latnode_n_entries(node) > 0) {
        int j;
        /* Assign unigram LM probabilities to backoff's outgoing
         * arcs (we did all the other backoff probs above when we
         * generated the other backoff nodes). */
        for (j = 0; j < ms_latnode_n_exits(node); ++j) {
            int32 linkid = garray_ent(node->exits, int32, j);
            ms_latlink_t *link = garray_ptr(l->link_list, ms_latlink_t, linkid);
            if (link->wid == dict_startwid(l->dict))
                link->lscr = 0;
            else
                link->lscr = ngram_ng_score(lm,
                                            map_lmwid(l->dict, lm, link->wid),
                                            NULL, 0, NULL);
        }
        /* Change it to an epsilon node. */
        ms_lattice_set_node_id(l, node, node->id.sf, -1);
        /* If this is a final node then link it to the new final node
         * created above. */
        if (node == ms_lattice_get_end(l)) { 
            ms_latlink_t *end_link;
            ms_latnode_t *end;
            end = ms_lattice_get_node_idx(l, endid);
            end_link = ms_lattice_link(l, node, end,
                                       dict_finishwid(l->dict),
                                       0 /* FIXME: final_node_ascr. */);
            end_link->lscr = ngram_ng_score(lm, ngram_wid(lm, "</s>"),
                                            NULL, 0, NULL);
        }
    }
    else  {
        /* If the backoff node has no more entries (and is not the
         * start node), then mark it for death. */
        node->id.lmstate = 0xdeadbeef;
    }
}

static ms_latnode_t *
create_lmstate_node(ms_lattice_t *l, ms_latnode_t *node,
                    int32 lmstate, int32 lscr, int32 endid)
{
    int32 nodeid = ms_lattice_get_idx_node(l, node);
    ms_latnode_t *new;

    /* Copy outgoing arcs from original node to the
     * newly created node. */
    new = ms_lattice_node_init(l, node->id.sf, lmstate);
    /* That could cause memory to move so update this pointer. */
    node = ms_lattice_get_node_idx(l, nodeid);

    /* If this is a final node then also link it to the new
     * final node created above. */
    if (node == ms_lattice_get_end(l)) { 
        ms_latlink_t *end_link;
        ms_latnode_t *end;
        end = ms_lattice_get_node_idx(l, endid);
        end_link = ms_lattice_link(l, new, end,
                                   dict_finishwid(l->dict),
                                   0 /* FIXME: final_node_ascr. */);
        end_link->lscr = lscr;
    }

    return new;
}

static int32
build_lmstate(ms_lattice_t *l, ngram_model_t *lm,
              ms_latnode_t *node, int32 node_lmwid,
              int32 entry_idx, int32 *out_lscr, int32 *out_bowt)
{
    ms_latlink_t *entry;
    ms_latnode_t *src;
    int32 src_latwid, src_lmwid;
    int32 entry_latwid, entry_lmwid;
    int32 lmstate;
    ngram_iter_t *ni;
    int n_hist, i;

    /* Skip arcs with deleted source nodes. */
    entry = garray_ptr(l->link_list, ms_latlink_t, entry_idx);
    src = ms_lattice_get_node_idx(l, entry->src);
    if (src->id.lmstate == 0xdeadbeef)
        return 0xdeadbeef;

    /* Construct the maximum language model state for this
     * incoming arc.  This consists of the language model
     * state from its source node plus the arc's base word
     * ID. */
    n_hist = ms_lattice_get_lmstate_wids(l, src->id.lmstate,
                                         &src_latwid, l->lathist);
    src_lmwid = map_lmwid(l->dict, lm, src_latwid);
    for (i = 0; i < n_hist; ++i) {
        l->lmhist[i] = map_lmwid(l->dict, lm, l->lathist[i]);
        E_INFO("%d %s %d -> %d\n", i,
               dict_wordstr(l->dict, l->lathist[i]),
               l->lathist[i], l->lmhist[i]);
    }

    /* Build language model state. */
    rotate_lmstate(src_latwid, l->lathist, n_hist, l->max_n_hist);
    n_hist = rotate_lmstate(src_lmwid, l->lmhist, n_hist, l->max_n_hist);
    for (i = 0; i < n_hist; ++i)
        assert(l->lmhist[i] >= 0);
    /* Map link word ID to base ID and lm ID. */
    entry_latwid = dict_basewid(l->dict, entry->wid);
    entry_lmwid = map_lmwid(l->dict, lm, entry->wid);
    /* Rotate in incoming arc word. */
    rotate_lmstate(entry_latwid, l->lathist, n_hist, l->max_n_hist);
    n_hist = rotate_lmstate(entry_lmwid, l->lmhist, n_hist, l->max_n_hist);
    for (i = 0; i < n_hist; ++i)
        assert(l->lmhist[i] >= 0);
    /* Set up pointers for language model state iteration. */
    lmstate = -1;
    *out_lscr = *out_bowt = 0;
    while (n_hist > 0) {
        /* First check if it exists in the language model. */
        if ((ni = ngram_ng_iter(lm, node_lmwid, l->lmhist,
                                n_hist)) != NULL) {
            /* Yes, create an lmstate for the HISTORY ONLY. */
            if ((lmstate = ms_lattice_get_lmstate_idx
                 (l, l->lathist[0], l->lathist + 1, n_hist - 1)) == -1)
                lmstate = ms_lattice_lmstate_init(l, l->lathist[0],
                                                  l->lathist + 1, n_hist - 1);
            ngram_iter_get(ni, out_lscr, NULL);
            /* We are the knights who say... */
            ngram_iter_free(ni);
            break;
        }
        else {
            /* If not, find the backoff weight (if any) and
             * back off the history.  Repeat this until we
             * find an N-Gram or we run out of history. */
            assert(n_hist > 0);
            if ((ni = ngram_ng_iter
                 (lm, l->lmhist[0], l->lmhist + 1, n_hist - 1)) != NULL) {
                ngram_iter_get(ni, NULL, out_bowt);
                ngram_iter_free(ni);
            }
            else
                *out_bowt = 0;
            if (--n_hist == 0)
                lmstate = -1;
        }
    }
    return lmstate;
}

static ms_latlink_t *
expand_entry(ms_lattice_t *l, int32 nodeid, int32 linkid,
             int32 endid, int32 lmstate, int32 lscr, int32 bowt)
{
    ms_latlink_t *link;
    ms_latnode_t *node, *new;
    int k, duplicate;

    /* If we found an appropriate non-null language model
     * state, look for a node with the desired language model
     * state in this frame, and create one if it does not yet
     * exist. */
    node = ms_lattice_get_node_idx(l, nodeid);
    /* FIXME: This is problematic because we have abused lmstates for
     * node words.  I can't think of a situation where that would
     * actually cause a problem but it could happen. */
    new = ms_lattice_get_node_id(l, node->id.sf, lmstate);
    if (new == NULL) {
        new = create_lmstate_node(l, node, lmstate, lscr, endid);
        /* Creating a node can move memory so update this pointer. */
        node = ms_lattice_get_node_idx(l, nodeid);
    }
    merge_exits(l, new, node, lscr);
    /* This incoming link might be a duplicate of one already
     * entering the given node, if so, mark it as such, so we
     * can remove it below. */
    duplicate = FALSE;
    link = garray_ptr(l->link_list, ms_latlink_t, linkid);
    for (k = 0; new->entries
             && k < ms_latnode_n_entries(new); ++k) {
        ms_latlink_t *xlink = ms_latnode_get_entry(l, new, k);
        if (link->wid == xlink->wid && link->src == xlink->src) {
            /* Tropical semiring assumption... */
            if (link->ascr > xlink->ascr)
                xlink->ascr = link->ascr;
            duplicate = TRUE;
        }
    }
    if (duplicate)
        return NULL;
    else {
        ms_latlink_t *link;
        /* Otherwise, repoint link at the newly created/found
         * node and add the necessary backoff weight, if any. */
        link = garray_ptr(l->link_list, ms_latlink_t, linkid);
        link->dest = ms_lattice_get_idx_node(l, new);
        link->lscr += bowt;
        if (new->entries == NULL)
            new->entries = garray_init(0, sizeof(int32));
        garray_append(new->entries, &linkid);
        return link;
    }
}

static void
expand_node(ms_lattice_t *l, ngram_model_t *lm,
            ms_latnode_t *node, int32 endid)
{
    garray_t *keep_entries = garray_init(0, sizeof(int32)); 
    garray_t *duplicate_entries = garray_init(0, sizeof(int32));
    int j, n_entries;
    int32 nodeid, node_lmwid;

    /* Word IDs have already been pushed to arcs. */
    n_entries = ms_latnode_n_entries(node);
    nodeid = ms_lattice_get_idx_node(l, node);
    /* Input lattice has word IDs on the nodes - get the word ID. */
    ms_lattice_get_lmstate_wids(l, node->id.lmstate, &node_lmwid, NULL);
    node_lmwid = map_lmwid(l->dict, lm, node_lmwid);
    /* Keep track of which entry links are moved to newly created
     * nodes, and which are deleted as duplicates. */
    /* Expand this node with unique incoming N-gram histories. */
    for (j = 0; j < n_entries; ++j) {
        ms_latlink_t *link;
        int32 lscr, bowt;
        int32 lmstate;
        int32 linkid;

        node = ms_lattice_get_node_idx(l, nodeid);
        linkid = garray_ent(node->entries, int32, j);
        lmstate = build_lmstate(l, lm, node,
                                node_lmwid, linkid,
                                &lscr, &bowt);
        /* This means the source of this arc is a deleted node. */
        if (lmstate == 0xdeadbeef)
            continue;

        if (lmstate != -1)
            link = expand_entry(l, nodeid, linkid, endid, lmstate, lscr, bowt);
        else {
            /* If we run out of history do nothing, the original node
             * becomes a null backoff node, and this incoming arc gets
             * a backoff weight added to it. */
            link = garray_ptr(l->link_list, ms_latlink_t, linkid);
            link->lscr += bowt;
            garray_append(keep_entries, &linkid);
        }
        if (link == NULL)
            garray_append(duplicate_entries, &linkid);
    }

    /* Remove duplicate incoming arcs. */
    for (j = 0; j < garray_size(duplicate_entries); ++j) {
        int32 linkid = garray_ent(duplicate_entries, int32, j);
        ms_latlink_t *link = garray_ptr(l->link_list, ms_latlink_t, linkid);
        ms_latlink_unlink(l, link);
    }
    garray_free(duplicate_entries);
    /* Remove copied incoming arcs from backoff node. */
    node = garray_ptr(l->node_list, ms_latnode_t, nodeid);
    garray_free(node->entries);
    if (garray_size(keep_entries) > 0)
        node->entries = keep_entries;
    else {
        node->entries = NULL;
        garray_free(keep_entries);
    }
    update_backoff_arcs(l, node, lm, endid);
}

static int
ms_lattice_alloc_hist(ms_lattice_t *l, int max_n_hist)
{
    if (max_n_hist > l->max_n_hist) {
        l->lmhist = ckd_realloc(l->lmhist,
                                max_n_hist * sizeof(*l->lmhist));
        l->lathist = ckd_realloc(l->lathist,
                                 max_n_hist * sizeof(*l->lathist));
        memset(l->lmhist, -1, max_n_hist * sizeof(*l->lmhist));
        memset(l->lathist, -1, max_n_hist * sizeof(*l->lathist));
    }
    l->max_n_hist = max_n_hist;
    return l->max_n_hist;
}

int
ms_lattice_expand(ms_lattice_t *l, ngram_model_t *lm)
{
    ms_latnode_t *end;
    ms_latnode_iter_t *itor;
    int32 endid, i;

    /* Storage for generating language model state.  We need to keep
     * parallel LM states for language model and dictionary word IDs.
     * The lattice uses dictionary word IDs, but we need language
     * model word IDs to get language model scores.  Once expansion is
     * complete the LM word IDs aren't needed. */
    ms_lattice_alloc_hist(l, ngram_model_get_size(lm) - 1);

    /* New final node (in theory, there are a bunch of different final
     * language model histories, but since nothing follows this there
     * is no need to distinguish them) */
    end = ms_lattice_node_init
        (l, ms_lattice_get_n_frames(l),
         ms_lattice_get_lmstate_idx(l, dict_finishwid(l->dict), NULL, 0));
    /* Have to use its ID because adding nodes can move memory. */
    endid = ms_lattice_get_idx_node(l, end);

    /* Traverse nodes. */
    for (itor = ms_lattice_traverse_topo(l, NULL);
         itor; itor = ms_latnode_iter_next(itor)) {
        expand_node(l, lm, ms_latnode_iter_get(itor), endid);
    }

    /* Update final node. */
    end = ms_lattice_get_node_idx(l, endid);
    ms_lattice_set_end(l, end);

    /* Remove nodes marked for death. */
    for (i = 0; i < garray_size(l->node_list); ++i) {
        ms_latnode_t *node = ms_lattice_get_node_idx(l, i);
        if (node->id.lmstate == 0xdeadbeef)
            ms_latnode_unlink(l, node);
    }

    /* Clean up unreachable nodes. */
    ms_lattice_unlink_unreachable(l);

    return 0;
}

int32
ms_lattice_forward(ms_lattice_t *l, int32 inv_aw)
{
    ms_latnode_iter_t *itor;
    ms_latnode_t *end;
    int i;

    for (itor = ms_lattice_traverse_topo(l, NULL);
         itor; itor = ms_latnode_iter_next(itor)) {
        ms_latnode_t *n = ms_latnode_iter_get(itor);
        int i, j;
        for (i = 0; i < ms_latnode_n_exits(n); ++i) {
            ms_latlink_t *wx = ms_latnode_get_exit(l, n, i);
            int32 forward = logmath_get_zero(l->lmath);
            for (j = 0; j < ms_latnode_n_entries(n); ++j) {
                ms_latlink_t *vx = ms_latnode_get_entry(l, n, j);
                forward = logmath_add(l->lmath, forward, vx->alpha);
            }
            if (ms_latnode_n_entries(n) == 0)
                forward = 0;
            wx->alpha = forward + wx->lscr + wx->ascr / inv_aw;
            E_DEBUG(2,("<Link: %s %d -> %d> alpha = %f\n",
                       dict_basestr(l->dict, wx->wid),
                       ms_lattice_get_node_idx(l, wx->src)->id.sf,
                       ms_lattice_get_node_idx(l, wx->dest)->id.sf,
                       logmath_log_to_ln(l->lmath, wx->alpha)));
        }
    }
    end = ms_lattice_get_end(l);
    l->norm = logmath_get_zero(l->lmath);
    for (i = 0; i < ms_latnode_n_entries(end); ++i) {
        ms_latlink_t *vx = ms_latnode_get_entry(l, end, i);
        l->norm = logmath_add(l->lmath, l->norm, vx->alpha);
    }
    E_DEBUG(2,("norm = %d = %f\n", l->norm, logmath_log_to_ln(l->lmath, l->norm)));
    return l->norm;
}

int32
ms_lattice_backward(ms_lattice_t *l, int32 inv_aw)
{
    ms_latnode_iter_t *itor;

    for (itor = ms_lattice_reverse_topo(l, NULL);
         itor; itor = ms_latnode_iter_next(itor)) {
        ms_latnode_t *n = ms_latnode_iter_get(itor);
        int i, j;
        for (i = 0; i < ms_latnode_n_entries(n); ++i) {
            ms_latlink_t *vx = ms_latnode_get_entry(l, n, i);
            vx->beta = logmath_get_zero(l->lmath);
            for (j = 0; j < ms_latnode_n_exits(n); ++j) {
                ms_latlink_t *wx = ms_latnode_get_exit(l, n, j);
                vx->beta = logmath_add(l->lmath, vx->beta,
                                       wx->beta + wx->lscr + wx->ascr / inv_aw);
            }
            if (ms_latnode_n_exits(n) == 0)
                vx->beta = 0;
            E_DEBUG(2,("<Link: %s %d -> %d> beta = %f\n",
                       dict_basestr(l->dict, vx->wid),
                       ms_lattice_get_node_idx(l, vx->src)->id.sf,
                       ms_lattice_get_node_idx(l, vx->dest)->id.sf,
                       logmath_log_to_ln(l->lmath, vx->beta)));
        }
    }
    return l->norm;
}

int
ms_latnode_print(FILE *fh, ms_lattice_t *l, ms_latnode_t *n)
{
    if (n->id.lmstate == 0xdeadbeef)
        return fprintf(fh, "0xdeadbeef");
    if (n->id.lmstate != -1) {
        int32 wid, i, n_hist;
        n_hist = ms_lattice_get_lmstate_wids
            (l, n->id.lmstate, &wid, l->lmhist);
        fprintf(fh, "%s", dict_wordstr(l->dict, wid));
        for (i = 0; i < n_hist; ++i)
            fprintf(fh, ",%s", dict_wordstr(l->dict, l->lmhist[i]));
        return fprintf(fh, "/%d", n->id.sf);
    }
    else
        return fprintf(fh, "&epsilon;/%d", n->id.sf);
}

int
ms_latlink_print(FILE *fh, ms_lattice_t *l, ms_latlink_t *vx)
{
    return fprintf
        (fh, "<Link: %s %d -> %d>",
         dict_basestr(l->dict, vx->wid),
         ms_lattice_get_node_idx(l, vx->src)->id.sf,
         ms_lattice_get_node_idx(l, vx->dest)->id.sf);
}
