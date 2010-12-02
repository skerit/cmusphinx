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

#include <sphinxbase/pio.h>
#include <sphinxbase/err.h>
#include <sphinxbase/strfuncs.h>
#include <sphinxbase/ckd_alloc.h>

#include "ms_lattice.h"

struct ms_lattice_s {
    int refcount;
    /**
     * Log base calculation.
     */
    logmath_t *lmath;
    /**
     * Dictionary for mapping language model states.
     */
    dict_t *dict;
    /**
     * Do we create new dictionary entries as needed?
     */
    int autodict;
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
    int32 start_idx;
    int32 end_idx;
    /**
     * Total probability of this lattice (normalizer for posteriors)
     */
    int32 norm;
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
    if (l == NULL)
        return 0;
    if (--l->refcount > 0)
        return l->refcount;

    garray_free(l->node_list);
    nodeid_map_free(l->node_map);
    logmath_free(l->lmath);
    dict_free(l->dict);
    ckd_free(l);
    return 0;
}

int32
ms_lattice_lmstate_init(ms_lattice_t *l, char const *name)
{
    int32 idx;

    if ((idx = dict_wordid(l->dict, name)) != -1)
        return idx;
    return dict_add_word(l->dict, name, NULL, 0);
}

int32
ms_lattice_get_lmstate_idx(ms_lattice_t *l, char const *name)
{
    return dict_wordid(l->dict, name);
}

char const *
ms_lattice_get_lmstate_name(ms_lattice_t *l, int32 idx)
{
    return dict_wordstr(l->dict, idx);
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

    return node;
}

int32
ms_lattice_get_idx_node(ms_lattice_t *l, ms_latnode_t *node)
{
    return garray_idx(l->node_list, node);
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
    node->id.lmstate = wid; /* This may get updated later. */
    node->fan = 0;
    node->exits = node->entries = NULL;
    nodeid_map_add(l->node_map, sf, wid, nodeidx);

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
    else
        /* NOTE: push forward word IDs immediately (might as well). */
        wid = src->id.lmstate;
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
    return 0;

error_out:
    lineiter_free(li);
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
            char const *basestr = dict_basestr(l->dict, node->id.lmstate);
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
    if (node->id.lmstate != -1)
        return fprintf(fh, " \"%s/%d\"",
                       dict_wordstr(l->dict, node->id.lmstate),
                       node->id.sf);
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
        for (j = 0; j < garray_size(node->exits); ++j) {
            int32 linkid = garray_ent(node->exits, int32, j);
            ms_latlink_t *link = garray_ptr(l->link_list, ms_latlink_t, linkid);
            ms_latnode_t *node2;
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

            fprintf(fh, "\t");
            node2 = ms_lattice_get_node_idx(l, link->src);
            print_dot_nodeid(l, node2, fh);
            fprintf(fh, " ->");
            node2 = ms_lattice_get_node_idx(l, link->dest);
            print_dot_nodeid(l, node2, fh);

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

int
ms_lattice_bigram_expand(ms_lattice_t *l, ngram_model_t *lm)
{
    return 0;
}

int32
ms_lattice_forward(ms_lattice_t *l, int32 inv_aw)
{
    return l->norm;
}

int32
ms_lattice_backward(ms_lattice_t *l, int32 inv_aw)
{
    return l->norm;
}

