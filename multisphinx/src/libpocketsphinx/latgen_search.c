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
    logmath_t *lmath;
    arc_buffer_t *input_arcs;
    ms_lattice_t *output_lattice;
    int32 incomplete;
} latgen_search_t;

ps_search_t *
latgen_init(cmd_ln_t *config,
	    dict2pid_t *d2p,
            logmath_t *lmath,
	    arc_buffer_t *input_arcs)
{
    latgen_search_t *latgen;

    latgen = ckd_calloc(1, sizeof(*latgen));
    ps_search_init(&latgen->base, &latgen_funcs,
                   config, NULL, d2p->dict, d2p);
    latgen->input_arcs = arc_buffer_retain(input_arcs);
    latgen->lmath = logmath_retain(lmath);
	
    return &latgen->base;
}

static int
latgen_search_process_arcs(latgen_search_t *latgen,
                           sarc_t *itor, int32 frame_idx)
{
    ms_latnode_t *src_incomplete;
    int n_arc;
    /* Find the incomplete node for this start frame.  If said node
     * does not exist, then do nothing, since no future arcs can exist
     * with an end frame of frame_idx-1, and thus no links will be made
     * to any nodes created in this frame (bptbl gc gets rid of a lot
     * of these spurious arcs but not all of them since it only
     * operates locally). */
    if ((src_incomplete = ms_lattice_get_node_id
         (latgen->output_lattice, frame_idx, latgen->incomplete)) == NULL)
        return 0;

    for (n_arc = 0; itor; ++n_arc, 
             itor = (sarc_t *)arc_buffer_iter_next
             (latgen->input_arcs, &itor->arc)) {
        ms_latnode_t *src, *dest;
        int32 lmstate;

        if (itor->arc.src != frame_idx) /* FIXME: iterators don't work like they should */
            continue;
        /* Create or find a "language model state" (actually not a
         * language model state, just a node identifier). */
        if ((lmstate = ms_lattice_get_lmstate_idx
             (latgen->output_lattice, itor->arc.wid, NULL, 0)) == -1)
            lmstate = ms_lattice_lmstate_init
                (latgen->output_lattice, itor->arc.wid, NULL, 0);

        E_INFO("Input arc %s / %d -> %d / %d\n",
               dict_wordstr(ps_search_dict(latgen), itor->arc.wid),
               itor->arc.src, itor->arc.dest + 1, itor->score);
        /* Look for a node to extend with this arc. */
        if ((src = ms_lattice_get_node_id
             (latgen->output_lattice, frame_idx, lmstate)) != NULL) {
            /* Extend said node with this arc, copying relevant input
             * arcs from the incomplete source node. */
        }
        else {
            /* Copy the incomplete node and all its input arcs,
             * generating appropriate acoustic scores for all of them
             * based on the initial phone of this outgoing arc. */
        }

        /* Get or create the incomplete node for the destination
         * frame.  Note that the arc buffer stores *inclusive* end
         * frame indices, since they are derived from backpointers. */
        if ((dest = ms_lattice_get_node_id
             (latgen->output_lattice, itor->arc.dest + 1, latgen->incomplete)) == NULL) {
            dest = ms_lattice_node_init
                (latgen->output_lattice, itor->arc.dest + 1, latgen->incomplete);

            int32 wid;
            ms_lattice_get_lmstate_wids(latgen->output_lattice,
                                        dest->id.lmstate, &wid, NULL);
            E_INFO("Created destination node %s/%d\n",
                   dict_wordstr(ps_search_dict(latgen), wid),
                   dest->id.sf);
        }
        else {
        }
    }
    return n_arc;
}

static int
latgen_search_decode(ps_search_t *base)
{
    latgen_search_t *latgen = (latgen_search_t *)base;
    ms_latnode_t *start;
    int32 incomplete_wid;
    int frame_idx;

    frame_idx = 0;
    E_INFO("waiting for arc buffer start\n");
    if (arc_buffer_consumer_start_utt(latgen->input_arcs, -1) < 0)
        return -1;
    latgen->output_lattice = ms_lattice_init(latgen->lmath,
                                             ps_search_dict(base));
    /* Create a special language model state for incomplete nodes. */
    incomplete_wid = dict_add_word(ps_search_dict(base),
                                   "<incomplete>", NULL, 0);
    latgen->incomplete = ms_lattice_lmstate_init
        (latgen->output_lattice, incomplete_wid, NULL, 0);
    /* And create a start node (which will always be <s>, but we will
     * use <incomplete> anyway in case something changes) */
    start = ms_lattice_node_init(latgen->output_lattice,
                                 0, latgen->incomplete);
    ms_lattice_set_start(latgen->output_lattice, start);
    while (arc_buffer_consumer_wait(latgen->input_arcs, -1) >= 0) {
        /* Process any incoming arcs.  For the time being the arc
         * buffer is treated essentially like a bptbl that has been
         * forward-sorted for us.  N-Gram expansion and
         * determinization are done by the lattice code. */
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
            E_INFO("Added %d arcs leaving frame %d\n", n_arc, frame_idx);
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
