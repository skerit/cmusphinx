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

#include "ps_search.h"
#include "arc_buffer.h"
#include "nodeid_map.h"
#include "ms_lattice.h"

static int latgen_search_decode(ps_search_t *base);
static int latgen_search_free(ps_search_t *base);
static char const *latgen_search_hyp(ps_search_t *base, int32 *out_score);
static int32 latgen_search_prob(ps_search_t *base);
static ps_seg_t *latgen_search_seg_iter(ps_search_t *base, int32 *out_score);

static ps_searchfuncs_t latgen_funcs = {
    /* name: */   "latgen",
    /* free: */   latgen_search_free,
    /* decode: */ latgen_search_decode,
    /* hyp: */      latgen_search_hyp,
    /* prob: */     latgen_search_prob,
    /* seg_iter: */ latgen_search_seg_iter,
};

typedef struct latgen_search_s {
    ps_search_t base;
    ngram_model_t *lm;
    dict2pid_t *d2p;
    logmath_t *lmath;
    arc_buffer_t *input_arcs;
    ms_lattice_t *output_lattice;
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
} latgen_search_t;

ps_search_t *
latgen_init(cmd_ln_t *config,
	    dict2pid_t *d2p,
            ngram_model_t *lm,
	    arc_buffer_t *input_arcs)
{
    latgen_search_t *latgen;

    latgen = ckd_calloc(1, sizeof(*latgen));
    ps_search_init(&latgen->base, &latgen_funcs,
                   config, NULL, d2p->dict, d2p);
    latgen->d2p = dict2pid_retain(d2p);
    latgen->input_arcs = arc_buffer_retain(input_arcs);
    latgen->lmath = logmath_retain(ngram_model_get_lmath(lm));
    latgen->lm = ngram_model_retain(lm);

    latgen->max_n_hist = ngram_model_get_size(lm) - 1;
    latgen->lmhist = ckd_calloc(latgen->max_n_hist,
                                sizeof(*latgen->lmhist));
    latgen->active_nodes = garray_init(0, sizeof(int32));
    latgen->link_rcid = garray_init(0, sizeof(uint8));
    latgen->link_altwid = garray_init(0, sizeof(int32));
    latgen->link_score = garray_init(0, sizeof(int32));
	
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
        E_INFO("Frame %d added node ", frame_idx);
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
                          bitvec_t *active_incoming_links,
                          sarc_t *arc)
{
    ms_latnode_t *node;
    int32 srcidx, destidx;
    int32 incoming_score, incoming_linkid;
    int32 lmstate, bo_lmstate, headwid, lscr, bowt;
    int ciphone, n_hist, i;
    xwdssid_t *rssid;
    int n_links = 0;

    /* Find best incoming link (to get starting path score) */
    ciphone = dict_first_phone(latgen->d2p->dict, arc->arc.wid);
    incoming_linkid = -1;
    incoming_score = WORST_SCORE;
    node = ms_lattice_get_node_idx(latgen->output_lattice, nodeidx);
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
            rssid = dict2pid_rssid
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
    /* Start node has no incoming link (duh) */
    assert(node->id.sf == 0 || incoming_linkid != -1);
    /* Now mark this incoming link as active. */
    if (i < ms_latnode_n_entries(node))
        bitvec_set(active_incoming_links, i);

    /* Create new language model state.  This is the language model
     * state of the source node plus the arc's base word ID.  */
    n_hist =
        ms_lattice_get_lmstate_wids(latgen->output_lattice,
                                    node->id.lmstate,
                                    &headwid, latgen->lmhist);
    /* First rotate the previous head word into the history. */
    rotate_lmstate(headwid, latgen->lmhist, n_hist,
                   latgen->max_n_hist);
    /* headwid + latgen->lmhist are now the raw lmstate. */
    headwid = dict_basewid(latgen->d2p->dict, arc->arc.wid);
    E_INFO("Raw language model state ");
    print_lmstate(err_get_logfp(), latgen->lm, headwid, latgen->lmhist, n_hist);
    E_INFOCONT("\n");

    /* Create or find a destination node with the above language model
     * state.  If none exists, back off until one is found.  (FIXME: Do
     * we have to create backoff nodes for each stage of backoff or
     * just the last one?)
     */
    lmstate = bo_lmstate = -1;
    lscr = bowt = 0;
    while (n_hist >= 0) {
        ngram_iter_t *ni = NULL;
        /* <s> is not an N-gram but it is an acceptable LM state. */
        if (headwid == dict_startwid(latgen->d2p->dict)
            || (ni = ngram_ng_iter(latgen->lm, headwid,
                                    latgen->lmhist, n_hist)) != NULL) {
            /* Create or find the relevant lmstate. */
            if ((lmstate = ms_lattice_get_lmstate_idx
                 (latgen->output_lattice, headwid,
                  latgen->lmhist, n_hist)) == -1)
                lmstate = ms_lattice_lmstate_init
                    (latgen->output_lattice, headwid,
                     latgen->lmhist, n_hist);
            if (ni) {
                ngram_iter_get(ni, &lscr, NULL);
                ngram_iter_free(ni);
            }
            break;
        }
        else {
            assert(n_hist > 0);
            /* Get the backoff weight. */
            if ((ni = ngram_ng_iter
                 (latgen->lm, latgen->lmhist[0],
                  latgen->lmhist + 1, n_hist - 1)) != NULL) {
                ngram_iter_get(ni, NULL, &bowt);
                ngram_iter_free(ni);
            }
            else
                bowt = 0;
            --n_hist;
        }
    }
    /* Epsilon is not appropriate for the destination node. */
    assert(lmstate != -1);
    E_INFO("Final language model state ");
    print_lmstate(err_get_logfp(), latgen->lm, headwid, latgen->lmhist, n_hist);
    E_INFOCONT("\n");
    /* Find or create a source node for the backoff language model
     * state.  If we create a source node copy incoming arcs, adding
     * the backoff weight to their lscr.  (FIXME: Verify that this is
     * the correct thing to do) */
    if (n_hist == 0) {
        E_INFO("Backoff language model state &epsilon bowt %d\n", bowt);
        if ((node = ms_lattice_get_node_id
             (latgen->output_lattice, node->id.sf, -1)) == NULL) {
            node = ms_lattice_node_init
                (latgen->output_lattice, node->id.sf, -1);
        }
    }
    else {
        E_INFO("Backoff language model state ");
        print_lmstate(err_get_logfp(), latgen->lm, latgen->lmhist[0],
                      latgen->lmhist + 1, n_hist - 1);
        E_INFOCONT("bowt %d\n", bowt);
        if ((node = ms_lattice_get_node_id
             (latgen->output_lattice, node->id.sf, bo_lmstate)) == NULL) {
            node = ms_lattice_node_init
                (latgen->output_lattice, node->id.sf, bo_lmstate);
        }
    }
    srcidx = ms_lattice_get_idx_node(latgen->output_lattice, node);
    assert(srcidx >= 0);
    if (srcidx != nodeidx) {
        E_INFO("New source node:");
        ms_latnode_print(err_get_logfp(), latgen->output_lattice, node);
        E_INFOCONT("\n");
    }

    /* Find or create a destination node. */
    if ((node = ms_lattice_get_node_id
         /* NOTE: bptbl indices are inclusive, ours are not. */
         (latgen->output_lattice, arc->arc.dest + 1, lmstate)) == NULL) {
        node = ms_lattice_node_init
            (latgen->output_lattice, arc->arc.dest + 1, lmstate);
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
                               arc, headwid, arc->score, NO_RC);
        link->lscr = lscr; /* FIXME: See above re. bowt */
        E_INFO("Created non-rc link ");
        ms_latlink_print(err_get_logfp(),
                         latgen->output_lattice, link);
        E_INFOCONT("\n");
        ++n_links;
    }
    else {
        for (i = 0; i < arc_buffer_max_n_rc(latgen->input_arcs); ++i) {
            if (bitvec_is_set(arc->rc_bits, i)) {
                ms_latlink_t *link;
                link = create_new_link
                    (latgen, srcidx, destidx, incoming_linkid,
                     arc, headwid,
                     arc_buffer_get_rcscore(latgen->input_arcs, arc, i), i);
                link->lscr = lscr; /* FIXME: See above re. bowt */
                E_INFO("Created rc %d link ", i);
                ms_latlink_print(err_get_logfp(),
                                 latgen->output_lattice, link);
                E_INFOCONT("\n");
                ++n_links;
            }
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
        ms_latnode_t *node = ms_lattice_get_node_idx
            (latgen->output_lattice, nodeidx);
        bitvec_t *active_links;
        int node_n_links;

        /* FIXME: Should allocate this in latgen and grow as needed. */
        active_links = bitvec_alloc(ms_latnode_n_entries(node));
        node_n_links = create_outgoing_links_one(latgen, nodeidx,
                                                 active_links, arc);
        if (node_n_links == 0) {
            /* This node is a goner, prune it. */
            /* FIXME: Actually it's not clear this will happen. */
        }
        else {
            /* Prune dangling incoming links (with no active right
             * context) */
            /* FIXME: This is a kind of annoying way to have to do
             * this... overzealous encapsulation perhaps? */
            glist_t deadlinks;
            gnode_t *gn;
            int j;

            deadlinks = NULL;
            /* Re-up this pointer as it may have changed (!@#$#) */
            node = ms_lattice_get_node_idx
                (latgen->output_lattice, nodeidx);
            for (j = 0; j < ms_latnode_n_entries(node); ++j) {
                if (bitvec_is_set(active_links, j))
                    continue;
                deadlinks = glist_add_ptr
                    (deadlinks, ms_latnode_get_entry
                     (latgen->output_lattice, node, j));
            }
            for (gn = deadlinks; gn; gn = gnode_next(gn))
                ms_latlink_unlink(latgen->output_lattice,
                                  (ms_latlink_t *)gnode_ptr(gn));
            glist_free(deadlinks);
        }
        bitvec_free(active_links);
        n_links += node_n_links;
    }

    return n_links;
}

static int
latgen_search_process_arcs(latgen_search_t *latgen,
                           sarc_t *itor, int32 frame_idx)
{
    int n_arc;

    /* Get source nodes for these arcs. */
    if (get_frame_active_nodes(latgen->output_lattice,
                               latgen->active_nodes, frame_idx) == 0)
        return 0;

    /* Iterate over all arcs exiting in this frame */
    for (n_arc = 0; itor; itor = (sarc_t *)arc_buffer_iter_next
             (latgen->input_arcs, &itor->arc)) {
        /* See note in arc_buffer.h... */
        if (itor->arc.src != frame_idx)
            break;

        /* Create new outgoing links for each source node. */
        n_arc += create_outgoing_links(latgen, itor);
    }

    return n_arc;
}

static int
latgen_search_decode(ps_search_t *base)
{
    latgen_search_t *latgen = (latgen_search_t *)base;
    int frame_idx;

    frame_idx = 0;
    E_INFO("waiting for arc buffer start\n");
    if (arc_buffer_consumer_start_utt(latgen->input_arcs, -1) < 0)
        return -1;

    /* Create lattice and initial epsilon node. */
    latgen->output_lattice = ms_lattice_init(latgen->lmath,
                                             ps_search_dict(base));
    ms_lattice_node_init(latgen->output_lattice, 0, -1);

    /* Reset some internal arrays. */
    garray_reset(latgen->link_rcid);
    garray_reset(latgen->link_altwid);
    garray_reset(latgen->link_score);

    /* Process frames full of arcs. */
    while (arc_buffer_consumer_wait(latgen->input_arcs, -1) >= 0) {
        ptmr_start(&base->t);
        while (1) {
            arc_t *itor;
            int n_arc;

            arc_buffer_lock(latgen->input_arcs);
            itor = arc_buffer_iter(latgen->input_arcs, frame_idx);
            if (itor == NULL) {
                arc_buffer_unlock(latgen->input_arcs);
                break;
            }
            n_arc = latgen_search_process_arcs(latgen, (sarc_t *)itor, frame_idx);
            E_INFO("Added %d links leaving frame %d\n", n_arc, frame_idx);
            arc_buffer_unlock(latgen->input_arcs);
            ++frame_idx;
        }
        ptmr_stop(&base->t);
        if (arc_buffer_eou(latgen->input_arcs)) {
            E_INFO("latgen: got EOU\n");
            /* Clean up the output lattice. */
            arc_buffer_consumer_end_utt(latgen->input_arcs);
            return frame_idx;
        }
    }
    return -1;
}

static int
latgen_search_free(ps_search_t *base)
{
    latgen_search_t *latgen = (latgen_search_t *)base;

    arc_buffer_free(latgen->input_arcs);
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
latgen_search_hyp(ps_search_t *base, int32 *out_score)
{
    return NULL;
}

/**
 * Bestpath search over the lattice.
 */
static ps_seg_t *
latgen_search_seg_iter(ps_search_t *base, int32 *out_score)
{
    return NULL;
}

/**
 * Forward-backward calculation over the lattice.
 */
static int32
latgen_search_prob(ps_search_t *base)
{
    return 0;
}
