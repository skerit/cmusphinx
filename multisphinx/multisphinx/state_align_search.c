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
 * @file state_align_search.c State (and phone and word) alignment search.
 */

#include <multisphinx/search_internal.h>

#include "state_align_search.h"
#include "tmat.h"

/**
 * State alignment search structure.
 */
struct state_align_search_s {
    search_t base;       /**< Base search structure. */
    hmm_context_t *hmmctx;  /**< HMM context structure. */
    alignment_t *al;     /**< Alignment structure being operated on. */
    hmm_t *hmms;            /**< Vector of HMMs corresponding to phone level. */
    int n_phones;	    /**< Number of HMMs (phones). */

    int frame;              /**< Current frame being processed. */
    int32 best_score;       /**< Best score in current frame. */

    int n_emit_state;       /**< Number of emitting states (tokens per frame) */
    uint16 *tokens;         /**< Tokens (backpointers) for state alignment. */
    int n_tok_alloc;         /**< Number of  tokens allocated. */
};
typedef struct state_align_search_s state_align_search_t;

static int
state_align_search_start(search_t *search)
{
    state_align_search_t *sas = (state_align_search_t *)search;

    /* Activate the initial state. */
    hmm_enter(sas->hmms, 0, 0, 0);

    return 0;
}

static void
renormalize_hmms(state_align_search_t *sas, int frame_idx, int32 norm)
{
    int i;
    for (i = 0; i < sas->n_phones; ++i)
        hmm_normalize(sas->hmms + i, norm);
}

static int32
evaluate_hmms(state_align_search_t *sas, int16 const *senscr, int frame_idx)
{
    int32 bs = WORST_SCORE;
    int i, bi;

    hmm_context_set_senscore(sas->hmmctx, senscr);

    bi = 0;
    for (i = 0; i < sas->n_phones; ++i) {
        hmm_t *hmm = sas->hmms + i;
        int32 score;

        if (hmm_frame(hmm) < frame_idx)
            continue;
        score = hmm_vit_eval(hmm);
        if (score BETTER_THAN bs) {
            bs = score;
            bi = i;
        }
    }
    return bs;
}

static void
prune_hmms(state_align_search_t *sas, int frame_idx)
{
    int nf = frame_idx + 1;
    int i;

    /* Check all phones to see if they remain active in the next frame. */
    for (i = 0; i < sas->n_phones; ++i) {
        hmm_t *hmm = sas->hmms + i;
        if (hmm_frame(hmm) < frame_idx)
            continue;
        hmm_frame(hmm) = nf;
    }
}

static void
phone_transition(state_align_search_t *sas, int frame_idx)
{
    int nf = frame_idx + 1;
    int i;

    for (i = 0; i < sas->n_phones - 1; ++i) {
        hmm_t *hmm, *nhmm;
        int32 newphone_score;

        hmm = sas->hmms + i;
        if (hmm_frame(hmm) != nf)
            continue;

        newphone_score = hmm_out_score(hmm);
        /* Transition into next phone using the usual Viterbi rule. */
        nhmm = hmm + 1;
        if (hmm_frame(nhmm) < frame_idx
            || newphone_score BETTER_THAN hmm_in_score(nhmm)) {
            hmm_enter(nhmm, newphone_score, hmm_out_history(hmm), nf);
        }
    }
}

/* FIXME: Use garray_t for this. */
#define TOKEN_STEP 20
static void
extend_tokenstack(state_align_search_t *sas, int frame_idx)
{
    if ((frame_idx + 1) * sas->n_emit_state >= sas->n_tok_alloc) {
        sas->n_tok_alloc = (frame_idx + TOKEN_STEP) * sas->n_emit_state;
        sas->tokens = ckd_realloc(sas->tokens,
                                  sas->n_tok_alloc * sizeof(*sas->tokens));
    }
    memset(sas->tokens + frame_idx * sas->n_emit_state, 0xff,
           sas->n_emit_state * sizeof(*sas->tokens));
}

static void
record_transitions(state_align_search_t *sas, int frame_idx)
{
    uint16 *tokens;
    int i;

    /* Push another frame of tokens on the stack. */
    extend_tokenstack(sas, frame_idx);
    tokens = sas->tokens + frame_idx * sas->n_emit_state;

    /* Scan all active HMMs */
    for (i = 0; i < sas->n_phones; ++i) {
        hmm_t *hmm = sas->hmms + i;
        int j;

        if (hmm_frame(hmm) < frame_idx)
            continue;
        for (j = 0; j < sas->hmmctx->n_emit_state; ++j) {
            int state_idx = i * sas->hmmctx->n_emit_state + j;
            /* Record their backpointers on the token stack. */
            tokens[state_idx] = hmm_history(hmm, j);
            /* Update backpointer fields with state index. */
            hmm_history(hmm, j) = state_idx;
        }
    }
}

static int
state_align_search_step(search_t *search)
{
    state_align_search_t *sas = (state_align_search_t *)search;
    acmod_t *acmod = search_acmod(search);
    int16 const *senscr;
    int i, frame_idx;

    if ((frame_idx = acmod_consumer_wait(acmod, -1)) < 0) {
        /* Normal end of utterance... */
        if (acmod_eou(acmod))
            return 0;
        /* Or something bad (i.e. cancellation)? */
        return -1;
    }

    /* Calculate senone scores. */
    for (i = 0; i < sas->n_phones; ++i)
        acmod_activate_hmm(acmod, sas->hmms + i);
    if ((senscr = acmod_score(acmod, frame_idx)) == NULL)
        return 0;

    /* Renormalize here if needed. */
    /* FIXME: Make sure to (unit-)test this!!! */
    if ((sas->best_score - 0x300000) WORSE_THAN WORST_SCORE) {
        E_INFO("Renormalizing Scores at frame %d, best score %d\n",
               frame_idx, sas->best_score);
        renormalize_hmms(sas, frame_idx, sas->best_score);
    }
    
    /* Viterbi step. */
    sas->best_score = evaluate_hmms(sas, senscr, frame_idx);
    prune_hmms(sas, frame_idx);

    /* Transition out of non-emitting states. */
    phone_transition(sas, frame_idx);

    /* Generate new tokens from best path results. */
    record_transitions(sas, frame_idx);

    /* Release the frame just searched. */
    acmod_consumer_release(acmod, frame_idx);

    /* Update frame counter */
    sas->frame = frame_idx;

    return 1;
}

static int
state_align_search_finish(search_t *search)
{
    state_align_search_t *sas = (state_align_search_t *)search;
    hmm_t *final_phone = sas->hmms + sas->n_phones - 1;
    alignment_iter_t *itor;
    alignment_entry_t *ent;
    int next_state, next_start, state, frame;

    /* Best state exiting the last frame. */
    next_state = state = hmm_out_history(final_phone);
    if (state == 0xffff) {
        E_ERROR("Failed to reach final state in alignment\n");
        return -1;
    }
    itor = alignment_states(sas->al);
    next_start = sas->frame + 1;
    for (frame = sas->frame - 1; frame >= 0; --frame) {
        state = sas->tokens[frame * sas->n_emit_state + state];
        /* State boundary, update alignment entry for next state. */
        if (state != next_state) {
            itor = alignment_iter_goto(itor, next_state);
            assert(itor != NULL);
            ent = alignment_iter_get(itor);
            ent->start = frame + 1;
            ent->duration = next_start - ent->start;
            E_DEBUG(1,("state %d start %d end %d\n", next_state,
                       ent->start, next_start));
            next_state = state;
            next_start = frame + 1;
        }
    }
    /* Update alignment entry for initial state. */
    itor = alignment_iter_goto(itor, 0);
    assert(itor != NULL);
    ent = alignment_iter_get(itor);
    ent->start = 0;
    ent->duration = next_start;
    E_DEBUG(1,("state %d start %d end %d\n", 0,
               ent->start, next_start));
    alignment_iter_free(itor);
    alignment_propagate(sas->al);

    /* Finalize the input acmod (signals producer) */
    acmod_consumer_end_utt(search->acmod);
    search->total_frames += search->acmod->output_frame;

    return 0;
}

static int
state_align_search_decode(search_t *base)
{
    int nfr, k;

    if (acmod_consumer_start_utt(base->acmod, -1) < 0)
        return -1;

    base->uttid = base->acmod->uttid;
    nfr = 0;
    state_align_search_start(base);
    while ((k = state_align_search_step(base)) > 0) {
        nfr += k;
    }

    /* Abnormal termination (cancellation, error) */
    if (k < 0) {
        if (base->output_arcs)
            arc_buffer_producer_shutdown(base->output_arcs);
        return k;
    }

    state_align_search_finish(base);
    return nfr;
}

static int
state_align_search_free(search_t *search)
{
    state_align_search_t *sas = (state_align_search_t *)search;
    if (sas == NULL)
        return 0;
    ckd_free(sas->hmms);
    ckd_free(sas->tokens);
    hmm_context_free(sas->hmmctx);
    return 0;
}

static char const *
state_align_search_hyp(search_t *base, int32 *out_score)
{
    state_align_search_t *sas = (state_align_search_t *)base;
    alignment_iter_t *itor;
    size_t len;

    if (out_score)
        *out_score = sas->best_score;

    if (sas->al == NULL)
        return NULL;

    ckd_free(base->hyp_str);
    len = 0;
    for (itor = alignment_words(sas->al); itor; itor = alignment_iter_next(itor)) {
        alignment_entry_t *ent = alignment_iter_get(itor);
        char *w = dict_wordstr(search_dict(base), ent->id.wid);
        len += strlen(w) + 2;
    }
    base->hyp_str = ckd_calloc(1, len);
    for (itor = alignment_words(sas->al); itor; itor = alignment_iter_next(itor)) {
        alignment_entry_t *ent = alignment_iter_get(itor);
        char *w = dict_wordstr(search_dict(base), ent->id.wid);
        strcat(base->hyp_str, w);
        strcat(base->hyp_str, " ");
    }
    base->hyp_str[strlen(base->hyp_str)-1] = '\0';

    return base->hyp_str;
}

static int32
state_align_search_prob(search_t *base)
{
    return 0;
}

typedef struct state_align_seg_s {
    struct seg_iter_s base;  /**< Base structure. */
    alignment_iter_t *itor;
} state_align_seg_t;

static void
state_align_bp2itor(state_align_seg_t *itor)
{
    alignment_entry_t *ent = alignment_iter_get(itor->itor);
    itor->base.wid= ent->id.wid;
    itor->base.sf = ent->start;
    itor->base.ef = ent->start + ent->duration - 1;
}

static void
state_align_seg_free(seg_iter_t *seg)
{
    state_align_seg_t *itor = (state_align_seg_t *)seg;

    alignment_iter_free(itor->itor);
    ckd_free(itor);
}

static seg_iter_t *
state_align_seg_next(seg_iter_t *seg)
{
    state_align_seg_t *itor = (state_align_seg_t *)seg;
    if ((itor->itor = alignment_iter_next(itor->itor)) == NULL) {
        state_align_seg_free(seg);
        return NULL;
    }
    state_align_bp2itor(itor);
    return seg;
}

static segfuncs_t state_align_segfuncs = {
    /* seg_next */ state_align_seg_next,
    /* seg_free */ state_align_seg_free
};


static seg_iter_t *
state_align_search_seg_iter(search_t *base, int32 *out_score)
{
    state_align_search_t *sas = (state_align_search_t *)base;
    state_align_seg_t *itor;

    if (sas->al == NULL)
        return NULL;

    if (out_score)
        *out_score = sas->best_score;

    itor = ckd_calloc(1, sizeof(*itor));
    itor->base.vt = &state_align_segfuncs;
    itor->base.lwf = 1.0;
    itor->base.search = base;
    if ((itor->itor = alignment_words(sas->al)) == NULL) {
        state_align_seg_free(&itor->base);
        return NULL;
    }
    state_align_bp2itor(itor);

    return &itor->base;
}

static search_t *state_align_search_init(search_t *other,
                        cmd_ln_t *config,
                        acmod_t *acmod,
                        dict2pid_t *d2p);

static searchfuncs_t state_align_search_funcs = {
    /* name: */   "state_align",
    /* init: */   state_align_search_init,
    /* free: */   state_align_search_free,
    /* decode: */ state_align_search_decode,
    /* hyp: */    state_align_search_hyp,
    /* prob: */   state_align_search_prob,
    /* seg_iter: */state_align_search_seg_iter
};

searchfuncs_t const *
state_align_search_query(void)
{
    return &state_align_search_funcs;
}

int
state_align_search_set_alignment(search_t *base, alignment_t *al)
{
    state_align_search_t *sas = (state_align_search_t *)base;
    alignment_iter_t *itor;
    hmm_t *hmm;

    alignment_free(sas->al);
    sas->al = alignment_retain(al);
    /* Generate HMM vector from phone level of alignment. */
    sas->n_phones = alignment_n_phones(al);
    sas->n_emit_state = alignment_n_states(al);
    sas->hmms = ckd_calloc(sas->n_phones, sizeof(*sas->hmms));
    for (hmm = sas->hmms, itor = alignment_phones(al); itor;
         ++hmm, itor = alignment_iter_next(itor)) {
        alignment_entry_t *ent = alignment_iter_get(itor);
        hmm_init(sas->hmmctx, hmm, FALSE,
                 ent->id.pid.ssid, ent->id.pid.tmatid);
    }
    return 0;
}

static search_t *
state_align_search_init(search_t *other,
                        cmd_ln_t *config,
                        acmod_t *acmod,
                        dict2pid_t *d2p)
{
    state_align_search_t *sas;

    sas = ckd_calloc(1, sizeof(*sas));
    search_base_init(search_base(sas), &state_align_search_funcs,
                   config, acmod, d2p);
    sas->hmmctx = hmm_context_init(bin_mdef_n_emit_state(acmod->mdef),
                                   acmod->tmat->tp, NULL, acmod->mdef->sseq);
    if (sas->hmmctx == NULL) {
        ckd_free(sas);
        return NULL;
    }
    return search_base(sas);
}
