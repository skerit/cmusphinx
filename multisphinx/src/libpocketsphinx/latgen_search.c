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
    ngram_model_t *lm;
    logmath_t *lmath;
    arc_buffer_t *input_arcs;
    ms_lattice_t *output_lattice;
    int32 *lmhist;
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
    latgen->input_arcs = arc_buffer_retain(input_arcs);
    latgen->lmath = logmath_retain(ngram_model_get_lmath(lm));
    latgen->lm = ngram_model_retain(lm);
    latgen->lmhist = ckd_calloc(ngram_model_get_size(lm)-1,
                                sizeof(*latgen->lmhist));
	
    return &latgen->base;
}

static int
latgen_search_process_arcs(latgen_search_t *latgen,
                           sarc_t *itor, int32 frame_idx)
{
    int n_arc;

    /* Iterate over all arcs exiting in this frame */
    for (n_arc = 0; itor; ++n_arc, 
             itor = (sarc_t *)arc_buffer_iter_next
             (latgen->input_arcs, &itor->arc)) {
        /* See note in arc_buffer.h... */
        if (itor->arc.src != frame_idx)
            break;
        /* Construct/find language model states associated with this
         * arc.  This is done by scanning all existing arcs pointing
         * to this start frame and rotating in their words. */
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
    latgen->output_lattice = ms_lattice_init(latgen->lmath,
                                             ps_search_dict(base));

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
    ngram_model_free(latgen->lm);
    ckd_free(latgen->lmhist);

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
