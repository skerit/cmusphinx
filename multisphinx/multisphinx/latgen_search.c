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
 * @file latgen_search.c Lattice generation (as a search pass).
 */

#include <sphinxbase/garray.h>
#include <sphinxbase/err.h>
#include <sphinxbase/filename.h>
#include <sphinxbase/strfuncs.h>
#include <sphinxbase/pio.h>

#include "search_internal.h"
#include "latgen_search.h"
#include "ms_lattice.h"
#include "nodeid_map.h"
#include "hmm.h"

static int latgen_search_decode(search_t *base);
static int latgen_search_free(search_t *base);
static char const *latgen_search_hyp(search_t *base, int32 *out_score);
static int32 latgen_search_prob(search_t *base);
static seg_iter_t *latgen_search_seg_iter(search_t *base, int32 *out_score);
static search_t *latgen_search_init(search_t *other, cmd_ln_t *config, acmod_t *acmod,
            dict2pid_t *d2p);


static searchfuncs_t latgen_funcs = {
    /* name: */   "latgen",
    /* init: */   latgen_search_init,
    /* free: */   latgen_search_free,
    /* decode: */ latgen_search_decode,
    /* hyp: */      latgen_search_hyp,
    /* prob: */     latgen_search_prob,
    /* seg_iter: */ latgen_search_seg_iter,
    /* bptbl: */  NULL,
    /* lmset: */ NULL
};

typedef struct latgen_search_s {
    search_t base;
    ngram_model_t *lm;
    dict2pid_t *d2p;
    logmath_t *lmath;
    arc_buffer_t *input_arcs;
    ms_lattice_t *output_lattice;

    /** Output directory for lattices. */
    char const *outlatdir;
    /** Output directory for arc buffer contents. */
    char const *outarcdir;
    int ctr;

    /** Storage for language model state components. */
    int32 *lmhist;
    /** Allocation size of @a lmhist. */
    int max_n_hist;
    /** List of active node IDs at current frame. */
    garray_t *active_nodes;
    /** Right context ID for all links in current lattice. */
    garray_t *link_rcid;
    /**
     * Original word ID corresponding to each link.
     *
     * We need to maintain this for building the lattice, because
     * links contain base word IDs (or rather language model word IDs)
     * but we need the correct word ID in order to find the correct
     * right context mapping.
     */
    garray_t *link_altwid;
    /** Raw path score for all links in current lattice. */
    garray_t *link_score;

    int32 silpen, fillpen;
} latgen_search_t;

searchfuncs_t const *
latgen_search_query(void)
{
    return &latgen_funcs;
}

search_t *
latgen_search_init(search_t *other,
            cmd_ln_t *config,
            acmod_t *acmod,
	    dict2pid_t *d2p)
{
    latgen_search_t *latgen;
    int32 wip;

    latgen = ckd_calloc(1, sizeof(*latgen));
    search_base_init(&latgen->base, &latgen_funcs,
                     config, NULL, d2p);
    latgen->d2p = dict2pid_retain(d2p);
    latgen->lmath = logmath_retain(acmod->lmath);

    latgen->outarcdir = cmd_ln_str_r(config, "-arcdumpdir");

    /* NOTE: this is one larger than the actual history size for a
     * language model state. */
    latgen->max_n_hist = 16; // FIXME
    latgen->lmhist = ckd_calloc(latgen->max_n_hist,
                                sizeof(*latgen->lmhist));
    latgen->active_nodes = garray_init(0, sizeof(int32));
    latgen->link_rcid = garray_init(0, sizeof(uint8));
    latgen->link_altwid = garray_init(0, sizeof(int32));
    latgen->link_score = garray_init(0, sizeof(int32));

    wip = logmath_log(latgen->lmath,
                           cmd_ln_float32_r(config, "-wip"))
        >> SENSCR_SHIFT;
    latgen->silpen = wip +
        (logmath_log(latgen->lmath,
                     cmd_ln_float32_r(config, "-silprob"))
         >> SENSCR_SHIFT);
    latgen->fillpen = wip +
        (logmath_log(latgen->lmath,
                     cmd_ln_float32_r(config, "-fillprob"))
         >> SENSCR_SHIFT);
	
    return &latgen->base;
}

/**
 * Construct the list of nodes active at this frame.
 */
static int
get_frame_active_nodes(ms_lattice_t *l, garray_t *out_active_nodes,
                       int32 frame_idx)
{
    ms_latnode_iter_t *node_itor;

    garray_reset(out_active_nodes);
    for (node_itor = ms_lattice_traverse_frame(l, frame_idx); node_itor;
         node_itor = ms_latnode_iter_next(node_itor)) {
        int32 node_idx = ms_latnode_iter_get_idx(node_itor);
        garray_append(out_active_nodes, &node_idx);
        E_INFO("Frame %d active node ", frame_idx);
        ms_latnode_print(err_get_logfp(), l,
                         ms_lattice_get_node_idx(l, node_idx));
        E_INFOCONT("\n");
    }
    return garray_size(out_active_nodes);
}

static int
print_lmstate(FILE *fh, ngram_model_t *lm, int32 headwid, int32 *lmhist,
              int n_hist)
{
    int i;

    if (fprintf(fh, "%s", ngram_word(lm, headwid)) < 0)
        return -1;
    for (i = 0; i < n_hist; ++i)
        if (fprintf(fh, ", %s", ngram_word(lm, lmhist[i])) < 0)
            return -1;
    return 0;
}

/**
 * Create a new link in the output lattice.
 */
static ms_latlink_t *
create_new_link(latgen_search_t *latgen,
                int32 srcidx, int32 destidx,
                int32 incoming_linkid,
                sarc_t *arc, int32 headwid, int32 score, int32 rc)
{
    ms_latnode_t *src, *dest;
    ms_latlink_t *link;
    int32 linkid, ascr;

    /* Calculate the acoustic score for this link. */
    ascr = score - arc->lscr;
    if (incoming_linkid != -1)
        ascr -= garray_ent(latgen->link_score, int32,
                           incoming_linkid);

    /* Create the new link. */
    /* FIXME: A matching link may already exist (with a different
     * alternate word ID) in which case we should just take the best
     * ascr. */
    src = ms_lattice_get_node_idx(latgen->output_lattice, srcidx);
    dest = ms_lattice_get_node_idx(latgen->output_lattice, destidx);
    link = ms_lattice_link(latgen->output_lattice,
                           src, dest, headwid, ascr);

    /* Record useful facts about this link for its successors. */
    linkid = ms_lattice_get_idx_link(latgen->output_lattice, link);
    garray_expand(latgen->link_rcid, linkid + 1);
    garray_ent(latgen->link_rcid, uint8, linkid) = rc;
    garray_expand(latgen->link_altwid, linkid + 1);
    garray_ent(latgen->link_altwid, int32, linkid) = arc->arc.wid;
    garray_expand(latgen->link_score, linkid + 1);
    garray_ent(latgen->link_score, int32, linkid) = score;

    return link;
}

/**
 * Create a backoff link.
 */
static ms_latlink_t *
create_backoff_link(latgen_search_t *latgen,
                    int32 linkid, int32 destidx,
                    int32 bowt)
{
    ms_latlink_t *link = ms_lattice_get_link_idx
        (latgen->output_lattice, linkid);
    ms_latlink_t *link2;
    ms_latnode_t *src, *dest;
    int32 link2id;

    E_INFO("Duplicating incoming link ");
    ms_latlink_print(err_get_logfp(),
                     latgen->output_lattice,
                     link);
    E_INFOCONT("\n");

    /* Create the new link. */
    /* FIXME: A matching link may already exist (with a different
     * alternate word ID) in which case we should just take the best
     * ascr. */
    src = ms_lattice_get_node_idx(latgen->output_lattice, link->src);
    dest = ms_lattice_get_node_idx(latgen->output_lattice, destidx);
    assert(src != NULL);
    assert(dest != NULL);
    link2 = ms_lattice_link(latgen->output_lattice,
                           src, dest, link->wid, link->ascr);
    link2->lscr = link->lscr + bowt;

    /* Record useful facts about this link for its successors. */
    link2id = ms_lattice_get_idx_link(latgen->output_lattice, link2);
    garray_expand(latgen->link_rcid, link2id + 1);
    garray_ent(latgen->link_rcid, uint8, link2id)
        = garray_ent(latgen->link_rcid, uint8, linkid);
    garray_expand(latgen->link_altwid, link2id + 1);
    garray_ent(latgen->link_altwid, int32, link2id)
        = garray_ent(latgen->link_altwid, int32, linkid);
    garray_expand(latgen->link_score, link2id + 1);
    garray_ent(latgen->link_score, int32, link2id)
        = garray_ent(latgen->link_score, int32, link2id);

    return link2;
}

static int32
find_best_incoming(latgen_search_t *latgen, int32 nodeidx,
                   sarc_t *arc, int32 *out_incoming_score)
{
    ms_latnode_t *node;
    int32 incoming_linkid, incoming_score;
    int i, ciphone;

    node = ms_lattice_get_node_idx(latgen->output_lattice, nodeidx);
    ciphone = dict_first_phone(latgen->d2p->dict,
                               arc->arc.wid);
    incoming_linkid = -1;
    incoming_score = WORST_SCORE;
    for (i = 0; i < ms_latnode_n_entries(node); ++i) {
        int32 linkid = ms_latnode_get_entry_idx
            (latgen->output_lattice, node, i);
        int rcid = garray_ent(latgen->link_rcid, uint8, linkid);
        int32 score = garray_ent(latgen->link_score, int32, linkid);

        /* No multiple right contexts: everything matches. */
        if (rcid == NO_RC) {
            if (score BETTER_THAN incoming_score) {
                incoming_linkid = linkid;
                incoming_score = score;
            }
        }
        /* Try to match ciphone against the right context ID. */
        else {
            int32 linkwid = garray_ent(latgen->link_altwid, int32, linkid);
            xwdssid_t *rssid = dict2pid_rssid
                (latgen->d2p,
                 dict_last_phone(latgen->d2p->dict, linkwid),
                 dict_second_last_phone(latgen->d2p->dict, linkwid));
            if (rssid->cimap[ciphone] == rcid) {
                if (score BETTER_THAN incoming_score) {
                    incoming_linkid = linkid;
                    incoming_score = score;
                }
            }
        }
    }
    *out_incoming_score = incoming_score;
    return incoming_linkid;
}

/**
 * Create lattice links for a given node and arc.
 *
 * 1) Find the incoming arc corresponding to the initial phone of
 *    this arc's word.
 * 2) Record the starting path score.
 * 3) Find the language model state for this arc's target.
 * 4) Find or create a node for that lmstate in the target frame.
 * 5) Create links to that node for all right contexts of this arc.
 */
static int
create_outgoing_links_one(latgen_search_t *latgen,
                          int32 nodeidx,
                          sarc_t *arc)
{
    ms_latnode_t *node;
    int32 srcidx, destidx;
    int32 incoming_score, incoming_linkid;
    int32 linkwid, lscr, bowt, headwid;
    int32 src_lmstate, dest_lmstate, bo_lmstate;
    int frame_idx, n_hist, i;
    int n_links = 0;

    node = ms_lattice_get_node_idx(latgen->output_lattice, nodeidx);
    src_lmstate = node->id.lmstate;
    frame_idx = node->id.sf;
    linkwid = dict_basewid(latgen->d2p->dict, arc->arc.wid);
    /* Might be modified if we have backoff. */
    srcidx = ms_lattice_get_idx_node(latgen->output_lattice, node);

    E_INFO("Original source node ");
    ms_latnode_print(err_get_logfp(), latgen->output_lattice, node);
    E_INFOCONT(" arc %s/%d\n", dict_wordstr(latgen->d2p->dict, linkwid),
               arc->arc.dest + 1);

    /* Find best incoming link (to get starting path score) */
    incoming_linkid = find_best_incoming
        (latgen, nodeidx, arc, &incoming_score);
    /* Start node has no incoming link (duh) */
    if (frame_idx > 0 && incoming_linkid == -1) {
        /* No matching incoming link for this arc, so we don't know
         * what its acoustic score should be.  Drop it on the floor.*/
        E_INFO("No incoming link found, skipping this arc\n");
        return 0;
    }

    /* Find language model state and language model score for
     * destination node.  */
    if (dict_filler_word(latgen->d2p->dict, arc->arc.wid)) {
        /* Except if this arc is a filler word in which case we
         * propagate the source node's LM state (which we assume is
         * valid by induction). */
        dest_lmstate = src_lmstate;
        if (arc->arc.wid == dict_silwid(latgen->d2p->dict))
            lscr = latgen->silpen;
        else
            lscr = latgen->fillpen;
        bowt = 0;

        /* NOTE: This is just for debugging purposes, all we actually
         * need is the language model state ID. */
        n_hist =
            ms_lattice_get_lmstate_wids(latgen->output_lattice,
                                        dest_lmstate,
                                        &headwid, latgen->lmhist);
        E_INFO("Filler %s, propagating language model state ",
               dict_wordstr(latgen->d2p->dict, arc->arc.wid));
        print_lmstate(err_get_logfp(), latgen->lm, headwid,
                      latgen->lmhist, n_hist);
        E_INFOCONT("\n");
    }
    else if (linkwid == dict_startwid(latgen->d2p->dict)) {
        /* Or if the arc is the start word ID in which case the
         * destination LM state is just <s>.  (this is done in part
         * because <s> is not actually in the LM). */
        dest_lmstate = ms_lattice_get_lmstate_idx
            (latgen->output_lattice, linkwid, NULL, 0);
        if (dest_lmstate == -1)
            dest_lmstate = ms_lattice_lmstate_init
                (latgen->output_lattice, linkwid, NULL, 0);
        lscr = bowt = 0;
        E_INFO("Start word <s> destination lmstate id %d\n", dest_lmstate);
    }
    else {
        /* This is the most complicated part of incremental lattice
         * construction.  At ths point we know what the previous and
         * next language model states should be if the language model
         * contains the N-Gram linking them.  If that N-Gram does not
         * exist, the previous language model state is invalid for
         * this arc, so we back it off until we find an appropriate
         * one.  The next language model state now becomes the
         * extension of this backed-off state with the arc's word, and
         * the incoming link is copied to a new backoff node and its
         * language model score is weighted with the backoff weight
         * obtained when backing off.
         *
         * As an example using trigrams, consider a source state (A,B)
         * and new arc C.  If (A,B,C) exists in the language model
         * then we rotate the lmstate and create a destination node
         * (B,C) with an arc to it.
         *
         * If (A,B,C) does not exist then we find the backoff weight
         * for (A,B) and back the incoming language model state off to
         * (B,), creating a new node for this and adding the backoff
         * weight to the best incoming arc which we copy to its
         * entries.  Then, if (B,C) exists in the language model, we
         * rotate the lmstate, creating a destination node (B,C) with
         * an arc from the backoff state to it.
         *
         * If (B,C) does not exist in the language model, we find the
         * backoff weight for (B,) and back the incoming language
         * model state off to &epsilon, creating a new node for this
         * and adding the backoff weight to the best incoming arc.
         * Then we rotate the lmstate creating a destination node (C,)
         * with an arc from the backoff state to it.
         */

        /* Get source language model state components. */
        n_hist =
            ms_lattice_get_lmstate_wids(latgen->output_lattice,
                                        src_lmstate,
                                        &headwid, latgen->lmhist);
        E_INFO("Source language model state %d ", src_lmstate);
        print_lmstate(err_get_logfp(), latgen->lm, headwid, latgen->lmhist, n_hist);
        E_INFOCONT("\n");
        /* Construct the target N-Gram (language model state + arc) */
        n_hist = 
            rotate_lmstate(headwid, latgen->lmhist, n_hist,
                           /* NOTE: This is bigger than an actual lm state. */
                           latgen->max_n_hist);
        headwid = linkwid;
        bo_lmstate = src_lmstate;
        while (n_hist >= 0) {
            ngram_iter_t *ni;
            E_INFO("Looking for N-Gram ");
            print_lmstate(err_get_logfp(), latgen->lm, headwid,
                          latgen->lmhist, n_hist);
            E_INFOCONT("\n");
            ni = ngram_ng_iter(latgen->lm, headwid, latgen->lmhist, n_hist);
            if (ni != NULL) {
                ngram_iter_get(ni, &lscr, NULL);
                ngram_iter_free(ni);
                E_INFO("Found: lscr %d\n", lscr);
                /* Destination language model state is this N-Gram,
                 * truncated to max_n_hist - 1 if necessary. */
                if (n_hist == latgen->max_n_hist)
                    --n_hist;
                dest_lmstate = ms_lattice_get_lmstate_idx
                    (latgen->output_lattice, headwid, latgen->lmhist, n_hist);
                if (dest_lmstate == -1)
                    dest_lmstate = ms_lattice_lmstate_init
                        (latgen->output_lattice, headwid,
                         latgen->lmhist, n_hist);
                E_INFO("Destination language model state %d: ", dest_lmstate);
                print_lmstate(err_get_logfp(), latgen->lm, headwid,
                              latgen->lmhist, n_hist);
                E_INFOCONT("\n");
                break;
            }
            else if (n_hist > 0) {
                --n_hist;
                E_INFO("Not found, looking for backoff N-Gram ");
                print_lmstate(err_get_logfp(), latgen->lm,
                              latgen->lmhist[0],
                              latgen->lmhist + 1, n_hist);
                E_INFOCONT("\n");
                /* Find backoff weight */
                ni = ngram_ng_iter(latgen->lm, latgen->lmhist[0],
                                   latgen->lmhist + 1, n_hist);
                if (ni != NULL) {
                    ngram_iter_get(ni, NULL, &bowt);
                    ngram_iter_free(ni);
                }
                else
                    bowt = 0;
                E_INFO("Backoff weight: %d\n", bowt);
                /* Now update source language model state. */
                if (n_hist == 0)
                    bo_lmstate = -1;
                else {
                    bo_lmstate = ms_lattice_get_lmstate_idx
                        (latgen->output_lattice, latgen->lmhist[0],
                         latgen->lmhist + 1, n_hist - 1);
                    if (bo_lmstate == -1)
                        bo_lmstate = ms_lattice_lmstate_init
                            (latgen->output_lattice, latgen->lmhist[0],
                             latgen->lmhist + 1, n_hist - 1);
                }
            }
            else {
                /* This implies that a unigram for headwid was not
                 * found, which should not happen. */
                E_ERROR("Unigram %s not found\n",
                        dict_wordstr(latgen->d2p->dict, headwid));
                assert(n_hist > 0);
            }
        }

        /* Find or create a source node for the backoff language model
         * state. */
        if (bo_lmstate == -1) {
            E_INFO("Backoff language model state &epsilon;\n");
            if ((node = ms_lattice_get_node_id
                 (latgen->output_lattice, frame_idx, -1)) == NULL) {
                node = ms_lattice_node_init
                    (latgen->output_lattice, frame_idx, -1);
            }
            srcidx = ms_lattice_get_idx_node(latgen->output_lattice, node);
        }
        else if (bo_lmstate != src_lmstate) {
            n_hist =
                ms_lattice_get_lmstate_wids(latgen->output_lattice,
                                            bo_lmstate,
                                            &headwid, latgen->lmhist);
            E_INFO("Backoff language model state %d ", bo_lmstate);
            print_lmstate(err_get_logfp(), latgen->lm, headwid,
                          latgen->lmhist, n_hist);
            E_INFOCONT("\n");
            if ((node = ms_lattice_get_node_id
                 (latgen->output_lattice, frame_idx, bo_lmstate)) == NULL) {
                node = ms_lattice_node_init
                    (latgen->output_lattice, frame_idx, bo_lmstate);
            }
            srcidx = ms_lattice_get_idx_node(latgen->output_lattice, node);
        }
        else {
            /* Source node is unchanged. */
            srcidx = nodeidx;
        }

        /* If we switched the source node then copy over the incoming link. */
        if (srcidx != nodeidx && incoming_linkid != -1)
            create_backoff_link(latgen, incoming_linkid,
                                srcidx, bowt >> SENSCR_SHIFT);
    }

    /* Find or create a destination node. */
    if ((node = ms_lattice_get_node_id
         /* NOTE: bptbl indices are inclusive, ours are not. */
         (latgen->output_lattice, arc->arc.dest + 1, dest_lmstate)) == NULL) {
        node = ms_lattice_node_init
            (latgen->output_lattice, arc->arc.dest + 1, dest_lmstate);
    }
    destidx = ms_lattice_get_idx_node(latgen->output_lattice, node);
    assert(destidx >= 0);
    E_INFO("Destination node ");
    ms_latnode_print(err_get_logfp(), latgen->output_lattice, node);
    E_INFOCONT("\n");

    /* For all right contexts create a link to dest. */
    if (dict_pronlen(latgen->d2p->dict, arc->arc.wid) == 1) {
        ms_latlink_t *link;
        link = create_new_link(latgen, srcidx, destidx, incoming_linkid,
                               arc, linkwid, arc->score, NO_RC);
        link->lscr = lscr >> SENSCR_SHIFT;
        E_INFO("Created non-rc link ");
        ms_latlink_print(err_get_logfp(),
                         latgen->output_lattice, link);
        E_INFOCONT("\n");
        ++n_links;
    }
    else {
        for (i = 0; i < arc_buffer_max_n_rc(search_input_arcs(latgen)); ++i) {
            ms_latlink_t *link;
            /* FIXME: Important observation, it seems that these are
             * always sequential, i.e. if we exit one of them we end
             * up exiting all of them.  Not sure this is generally
             * true but maybe it's useful for compression. */
            /* FIXME: ALSO, many of the rc deltas are exactly the
             * same.  So if we could compress them further (maybe
             * here) we could save a lot of spurious arcs. */
            if (!bitvec_is_set(arc->rc_bits, i))
                continue;
            /* FIXME: ALSO some of these scores are positive which
             * means that we are screwing up somehow. */
            link = create_new_link
                (latgen, srcidx, destidx, incoming_linkid,
                 arc, linkwid,
                 arc_buffer_get_rcscore(search_input_arcs(latgen), arc, i), i);
            link->lscr = lscr >> SENSCR_SHIFT;
            E_INFO("Created rc %d link ", i);
            ms_latlink_print(err_get_logfp(),
                             latgen->output_lattice, link);
            E_INFOCONT(" %d %d\n", link->ascr, link->lscr);
            ++n_links;
        }
    }

    return n_links;
}

/**
 * Create lattice links for a given arc.
 */
static int
create_outgoing_links(latgen_search_t *latgen,
                      sarc_t *arc)
{
    int n_links = 0;
    int i;

    for (i = 0; i < garray_size(latgen->active_nodes); ++i) {
        int32 nodeidx = garray_ent(latgen->active_nodes, int32, i);
        n_links += create_outgoing_links_one(latgen, nodeidx, arc);
    }

    return n_links;
}

static int
latgen_search_process_arcs(latgen_search_t *latgen,
                           sarc_t *itor, int32 frame_idx, FILE *arcfh)
{
    int n_arc;

    /* Get source nodes for these arcs. */
#if 0
    if (get_frame_active_nodes(latgen->output_lattice,
                               latgen->active_nodes, frame_idx) == 0)
        return 0;
#endif

    /* Iterate over all arcs exiting in this frame */
    for (n_arc = 0; itor; itor = (sarc_t *)arc_buffer_iter_next
             (search_input_arcs(latgen), &itor->arc)) {
        /* See note in arc_buffer.h... */
        if (itor->arc.src != frame_idx)
            break;
        if (arcfh) {
            rcdelta_t const *deltas, *d;
            int i;
            fprintf(arcfh, "%s %d %d %d %d",
                    dict_wordstr(latgen->d2p->dict, itor->arc.wid),
                    itor->arc.src, itor->arc.dest,
                    itor->score, itor->lscr);
            d = deltas = arc_buffer_get_rcdeltas(search_input_arcs(latgen), itor);
            for (i = 0; i < arc_buffer_max_n_rc(search_input_arcs(latgen)); ++i) {
                if (bitvec_is_set(itor->rc_bits, i))
                    fprintf(arcfh, " %d:%u", i, *d++);
            }
            fprintf(arcfh, "\n");
        }

        /* Create new outgoing links for each source node. */
        /* n_arc += create_outgoing_links(latgen, itor); */
        ++n_arc;
    }

    return n_arc;
}

static int
latgen_search_cleanup_frame(latgen_search_t *latgen, int32 frame_idx)
{
    int i, dead;

    /* Refresh the list of nodes in this frame. */
    if (get_frame_active_nodes(latgen->output_lattice,
                               latgen->active_nodes, frame_idx) == 0)
        return 0;

    /* Reap dangling nodes. */
    dead = 0;
    for (i = 0; i < garray_size(latgen->active_nodes); ++i) {
        int32 nodeidx = garray_ent(latgen->active_nodes, int32, i);
        ms_latnode_t *node = ms_lattice_get_node_idx
            (latgen->output_lattice, nodeidx);
        int32 wid;

        ms_lattice_get_lmstate_wids(latgen->output_lattice,
                                    node->id.lmstate, &wid,
                                    latgen->lmhist);
        if (node->id.sf != 0 && ms_latnode_n_entries(node) == 0) {
            ms_latnode_unlink(latgen->output_lattice, node);
            ++dead;
        }
        else if (wid != dict_finishwid(latgen->d2p->dict)
                 && ms_latnode_n_exits(node) == 0) {
            ms_latnode_unlink(latgen->output_lattice, node);
            ++dead;
        }
        else {
            /* For nodes that remain, reap dangling links. */
            /* A dangling link is an incoming link which does not
             * match any outgoing right contexts. */
            /* So basically we create a bitmap of initial phones for
             * each outgoing link.  Then for each incoming link we check that . */
        }
    }
    E_INFO("Cleaned up %d nodes in frame %d\n", dead, frame_idx);

    return 0;
}

static int
latgen_search_decode(search_t *base)
{
    latgen_search_t *latgen = (latgen_search_t *)base;
    int frame_idx;
    FILE *arcfh;

    frame_idx = 0;
    E_INFO("waiting for arc buffer start\n");
    if (arc_buffer_consumer_start_utt(search_input_arcs(latgen), -1) < 0)
        return -1;
    base->uttid = arc_buffer_uttid(search_input_arcs(latgen));

    /* Create lattice and initial epsilon node. */
    latgen->output_lattice = ms_lattice_init(latgen->lmath,
                                             search_dict(base));
    ms_lattice_node_init(latgen->output_lattice, 0, -1);

    /* Reset some internal arrays. */
    garray_reset(latgen->link_rcid);
    garray_reset(latgen->link_altwid);
    garray_reset(latgen->link_score);

    /* Start logging arcs. */
    if (latgen->outarcdir) {
        char *outfile;
        char *basedir;
        outfile = string_join(latgen->outarcdir, "/",
                              base->uttid, ".arc", NULL);
        basedir = ckd_salloc(outfile);
        path2dirname(outfile, basedir);
        build_directory(basedir);
        if ((arcfh = fopen(outfile, "w")) == NULL)
            E_FATAL_SYSTEM("WTF %s", outfile);
        ckd_free(basedir);
        ckd_free(outfile);
    }

    /* Process frames full of arcs. */
    while (arc_buffer_consumer_wait(search_input_arcs(latgen), -1) >= 0) {
        ptmr_start(&base->t);
        while (1) {
            arc_t *itor;
            int n_arc;

            /* Grab arcs from the input buffer. */
            arc_buffer_lock(search_input_arcs(latgen));
            itor = arc_buffer_iter(search_input_arcs(latgen), frame_idx);
            if (itor == NULL) {
                arc_buffer_unlock(search_input_arcs(latgen));
                break;
            }
            n_arc = latgen_search_process_arcs(latgen, (sarc_t *)itor, frame_idx, arcfh);
            arc_buffer_unlock(search_input_arcs(latgen));

            /* Release arcs, we don't need them anymore. */
            arc_buffer_consumer_release(search_input_arcs(latgen), frame_idx);
            /* Remove any inaccessible nodes in this frame. */
            latgen_search_cleanup_frame(latgen, frame_idx);
            ++frame_idx;
        }
        ptmr_stop(&base->t);
        if (arc_buffer_eou(search_input_arcs(latgen))) {
            E_INFO("latgen: got EOU\n");
            arc_buffer_consumer_end_utt(search_input_arcs(latgen));
            if (arcfh) fclose(arcfh);
            return frame_idx;
        }
    }
    if (arcfh) fclose(arcfh);
    return -1;
}

static int
latgen_search_free(search_t *base)
{
    latgen_search_t *latgen = (latgen_search_t *)base;

    logmath_free(latgen->lmath);
    dict2pid_free(latgen->d2p);
    ngram_model_free(latgen->lm);
    ckd_free(latgen->lmhist);
    garray_free(latgen->active_nodes);
    return 0;
}

/**
 * Bestpath search over the lattice.
 */
static char const *
latgen_search_hyp(search_t *base, int32 *out_score)
{
    return NULL;
}

/**
 * Bestpath search over the lattice.
 */
static seg_iter_t *
latgen_search_seg_iter(search_t *base, int32 *out_score)
{
    return NULL;
}

/**
 * Forward-backward calculation over the lattice.
 */
static int32
latgen_search_prob(search_t *base)
{
    return 0;
}
