/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008 Carnegie Mellon University.  All rights
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
 * @file ngram_search_fwdflat.c Flat lexicon search.
 */

/* System headers. */
#include <string.h>
#include <assert.h>

/* SphinxBase headers. */
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/listelem_alloc.h>
#include <sphinxbase/err.h>

/* Local headers. */
#include "fwdflat_search.h"

/* Turn this on to dump channels for debugging */
#define __CHAN_DUMP__		0
#if __CHAN_DUMP__
#define chan_v_eval(chan) hmm_dump_vit_eval(&(chan)->hmm, stderr)
#else
#define chan_v_eval(chan) hmm_vit_eval(&(chan)->hmm)
#endif

static int fwdflat_search_start(ps_search_t *base);
static int fwdflat_search_step(ps_search_t *base, int frame_idx);
static int fwdflat_search_finish(ps_search_t *base);
static int fwdflat_search_reinit(ps_search_t *base, dict_t *dict, dict2pid_t *d2p);
static void fwdflat_search_free(ps_search_t *base);
static char const *fwdflat_search_hyp(ps_search_t *base, int32 *out_score);
static int32 fwdflat_search_prob(ps_search_t *base);
static ps_seg_t *fwdflat_search_seg_iter(ps_search_t *base, int32 *out_score);

static ps_searchfuncs_t fwdflat_funcs = {
    /* name: */   "fwdflat",
    /* start: */  fwdflat_search_start,
    /* step: */   fwdflat_search_step,
    /* finish: */ fwdflat_search_finish,
    /* reinit: */ fwdflat_search_reinit,
    /* free: */   fwdflat_search_free,
    /* lattice: */  NULL,
    /* hyp: */      fwdflat_search_hyp,
    /* prob: */     fwdflat_search_prob,
    /* seg_iter: */ fwdflat_search_seg_iter,
};

static void build_fwdflat_chan(fwdflat_search_t *ffs);

static void
fwdflat_search_calc_beams(fwdflat_search_t *ffs)
{
    cmd_ln_t *config;
    acmod_t *acmod;

    config = ps_search_config(ffs);
    acmod = ps_search_acmod(ffs);

    /* Log beam widths. */
    ffs->fwdflatbeam = logmath_log(acmod->lmath,
                                   cmd_ln_float64_r(config, "-fwdflatbeam")) >> SENSCR_SHIFT;
    ffs->fwdflatwbeam = logmath_log(acmod->lmath,
                                    cmd_ln_float64_r(config, "-fwdflatwbeam")) >> SENSCR_SHIFT;

    /* Other things. */
    ffs->pip = logmath_log(acmod->lmath,
                           cmd_ln_float32_r(config, "-pip")) >> SENSCR_SHIFT;
    ffs->silpen = logmath_log(acmod->lmath,
                              cmd_ln_float32_r(config, "-silprob")) >> SENSCR_SHIFT;
    ffs->fillpen = logmath_log(acmod->lmath,
                               cmd_ln_float32_r(config, "-fillprob")) >> SENSCR_SHIFT;
    ffs->min_ef_width = cmd_ln_int32_r(ps_search_config(ffs), "-fwdflatefwid");
    ffs->max_sf_win = cmd_ln_int32_r(ps_search_config(ffs), "-fwdflatsfwin");
}

static void
fwdflat_search_update_widmap(fwdflat_search_t *ffs)
{
    const char **words;
    int32 i, n_words;

    /* It's okay to include fillers since they won't be in the LM */
    n_words = ps_search_n_words(ffs);
    words = ckd_calloc(n_words, sizeof(*words));
    /* This will include alternates, again, that's okay since they aren't in the LM */
    for (i = 0; i < n_words; ++i)
        words[i] = (const char *)dict_wordstr(ps_search_dict(ffs), i);
    ngram_model_set_map_words(ffs->lmset, words, n_words);
    ckd_free(words);
}

ps_search_t *
fwdflat_search_init(cmd_ln_t *config, acmod_t *acmod,
                    dict_t *dict, dict2pid_t *d2p)
{
    fwdflat_search_t *ffs;
    const char *path;

    ffs = ckd_calloc(1, sizeof(*ffs));
    ps_search_init(&ffs->base, &fwdflat_funcs, config, acmod, dict, d2p);
    ffs->hmmctx = hmm_context_init(bin_mdef_n_emit_state(acmod->mdef),
                                   acmod->tmat->tp, NULL, acmod->mdef->sseq);
    if (ffs->hmmctx == NULL) {
        ps_search_free(ps_search_base(ffs));
        return NULL;
    }
    ffs->chan_alloc = listelem_alloc_init(sizeof(internal_node_t));
    ffs->root_chan_alloc = listelem_alloc_init(sizeof(first_node_t));

    /* Calculate various beam widths and such. */
    fwdflat_search_calc_beams(ffs);

    /* Allocate a billion different tables for stuff. */
    ffs->word_chan = ckd_calloc(dict_size(dict),
                                sizeof(*ffs->word_chan));
    ffs->word_active = bitvec_alloc(dict_size(dict));
    ffs->word_idx = ckd_calloc(dict_size(dict),
                               sizeof(*ffs->word_idx));

    ffs->bptbl = bptbl_init(d2p, cmd_ln_int32_r(config, "-latsize"), 256);

    /* Allocate active word list array */
    ffs->active_word_list = ckd_calloc_2d(2, dict_size(dict),
                                          sizeof(**ffs->active_word_list));

    /* Load language model(s) */
    if ((path = cmd_ln_str_r(config, "-lmctl"))) {
        ffs->lmset = ngram_model_set_read(config, path, acmod->lmath);
        if (ffs->lmset == NULL) {
            E_ERROR("Failed to read language model control file: %s\n",
                    path);
            goto error_out;
        }
        /* Set the default language model if needed. */
        if ((path = cmd_ln_str_r(config, "-lmname"))) {
            ngram_model_set_select(ffs->lmset, path);
        }
    }
    else if ((path = cmd_ln_str_r(config, "-lm"))) {
        static const char *name = "default";
        ngram_model_t *lm;

        lm = ngram_model_read(config, path, NGRAM_AUTO, acmod->lmath);
        if (lm == NULL) {
            E_ERROR("Failed to read language model file: %s\n", path);
            goto error_out;
        }
        ffs->lmset = ngram_model_set_init(config,
                                          &lm, (char **)&name,
                                          NULL, 1);
        if (ffs->lmset == NULL) {
            E_ERROR("Failed to initialize language model set\n");
            goto error_out;
        }
    }
    if (ffs->lmset != NULL
        && ngram_wid(ffs->lmset, S3_FINISH_WORD) == ngram_unknown_wid(ffs->lmset)) {
        E_ERROR("Language model/set does not contain </s>, recognition will fail\n");
        goto error_out;
    }

    /* Create word mappings. */
    fwdflat_search_update_widmap(ffs);

    /* Create HMM network. */
    build_fwdflat_chan(ffs);

    return (ps_search_t *)ffs;

error_out:
    fwdflat_search_free((ps_search_t *)ffs);
    return NULL;
}

static void
fwdflat_search_free_all_rc(fwdflat_search_t *ffs, int32 w)
{
    internal_node_t *hmm, *thmm;

    hmm = ffs->word_chan[w]->next;
    hmm_deinit(&ffs->word_chan[w]->hmm);
    listelem_free(ffs->root_chan_alloc, ffs->word_chan[w]);
    while (hmm) {
        thmm = hmm->next;
        hmm_deinit(&hmm->hmm);
        listelem_free(ffs->chan_alloc, hmm);
        hmm = thmm;
    }
    ffs->word_chan[w] = NULL;
}

static void
destroy_fwdflat_chan(fwdflat_search_t *ffs)
{
    int32 wid;

    for (wid = 0; wid < ps_search_n_words(ffs); ++wid) {
        assert(ffs->word_chan[wid] != NULL);
        fwdflat_search_free_all_rc(ffs, wid);
    }
}

static void
fwdflat_search_free(ps_search_t *base)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *)base;

    ps_search_deinit(base);

    destroy_fwdflat_chan(ffs);

    hmm_context_free(ffs->hmmctx);
    listelem_alloc_free(ffs->chan_alloc);
    listelem_alloc_free(ffs->root_chan_alloc);
    ngram_model_free(ffs->lmset);

    ckd_free(ffs->word_idx);
    ckd_free(ffs->word_chan);
    bitvec_free(ffs->word_active);
    bptbl_free(ffs->bptbl);
    ckd_free_2d(ffs->active_word_list);
    ckd_free(ffs);
}

static int
fwdflat_search_reinit(ps_search_t *base, dict_t *dict, dict2pid_t *d2p)
{
    /* Reallocate things that depend on the number of words. */
    fwdflat_search_t *ffs = (fwdflat_search_t *)base;
    int old_n_words;

    old_n_words = ps_search_n_words(ffs);
    if (old_n_words != dict_size(dict)) {
        base->n_words = dict_size(dict);
        ckd_free(ffs->word_idx);
        ckd_free(ffs->word_active);
        ckd_free_2d(ffs->active_word_list);
        ckd_free(ffs->word_chan);

        ffs->word_idx = ckd_calloc(base->n_words,
                                   sizeof(*ffs->word_idx));
        ffs->word_active = bitvec_alloc(base->n_words);
        ffs->active_word_list
            = ckd_calloc_2d(2, base->n_words,
                            sizeof(**ffs->active_word_list));
        ffs->word_chan = ckd_calloc(base->n_words, sizeof(*ffs->word_chan));
    }

    /* Free old dict2pid, dict */
    ps_search_base_reinit(base, dict, d2p);
    /* Update beam widths. */
    fwdflat_search_calc_beams(ffs);
    /* Update word mappings. */
    fwdflat_search_update_widmap(ffs);
    /* Rebuild HMM network */
    build_fwdflat_chan(ffs);

    return 0;
}

static internal_node_t *
fwdflat_search_alloc_all_rc(fwdflat_search_t *ffs, int32 w)
{
    internal_node_t *fhmm, *hmm;
    xwdssid_t *rssid;
    int32 i, tmatid, ciphone;

    /* DICT2PID */
    /* Get pointer to array of triphones for final diphone. */
    assert(!dict_is_single_phone(ps_search_dict(ffs), w));
    ciphone = dict_last_phone(ps_search_dict(ffs),w);
    rssid = dict2pid_rssid(ps_search_dict2pid(ffs),
                           ciphone,
                           dict_second_last_phone(ps_search_dict(ffs),w));
    tmatid = bin_mdef_pid2tmatid(ps_search_acmod(ffs)->mdef, ciphone);
    fhmm = hmm = listelem_malloc(ffs->chan_alloc);
    hmm->rc_id = 0;
    hmm->next = NULL;
    hmm->ciphone = dict_last_phone(ps_search_dict(ffs),w);
    hmm_init(ffs->hmmctx, &hmm->hmm, FALSE, rssid->ssid[0], tmatid);
    E_DEBUG(3,("allocated rc_id 0 ssid %d ciphone %d lc %d word %s\n",
               rssid->ssid[0], hmm->ciphone,
               dict_second_last_phone(ps_search_dict(ffs),w),
               dict_wordstr(ps_search_dict(ffs),w)));
    for (i = 1; i < rssid->n_ssid; ++i) {
        if ((hmm->next == NULL) || (hmm_nonmpx_ssid(&hmm->next->hmm) != rssid->ssid[i])) {
            internal_node_t *thmm;
            thmm = listelem_malloc(ffs->chan_alloc);
            thmm->next = hmm->next;
            hmm->next = thmm;
            hmm = thmm;

            hmm->rc_id = i;
            hmm->ciphone = ciphone;
            hmm_init(ffs->hmmctx, &hmm->hmm, FALSE, rssid->ssid[i], tmatid);
            E_DEBUG(3,("allocated rc_id %d ssid %d ciphone %d lc %d word %s\n",
                       i, rssid->ssid[i], hmm->ciphone,
                       dict_second_last_phone(ps_search_dict(ffs),w),
                       dict_wordstr(ps_search_dict(ffs),w)));
        }
        else
            hmm = hmm->next;
    }

    return fhmm;
}

/**
 * Build HMM network.
 */
static void
build_fwdflat_chan(fwdflat_search_t *ffs)
{
    int32 wid, p;
    first_node_t *rhmm;
    internal_node_t *hmm, *prevhmm;
    dict_t *dict;
    dict2pid_t *d2p;

    dict = ps_search_dict(ffs);
    d2p = ps_search_dict2pid(ffs);

    /* Build word HMMs for each word in the dictionary. */
    for (wid = 0; wid < ps_search_n_words(ffs); ++wid) {
        assert(ffs->word_chan[wid] == NULL);

        /* Multiplex root HMM for first phone (one root per word, flat
         * lexicon). */
        rhmm = listelem_malloc(ffs->root_chan_alloc);
        if (dict_is_single_phone(dict, wid)) {
            rhmm->ciphone = dict_first_phone(dict, wid);
            rhmm->ci2phone = bin_mdef_silphone(ps_search_acmod(ffs)->mdef);
        }
        else {
            rhmm->ciphone = dict_first_phone(dict, wid);
            rhmm->ci2phone = dict_second_phone(dict, wid);
        }
        rhmm->next = NULL;
        hmm_init(ffs->hmmctx, &rhmm->hmm, TRUE,
                 bin_mdef_pid2ssid(ps_search_acmod(ffs)->mdef, rhmm->ciphone),
                 bin_mdef_pid2tmatid(ps_search_acmod(ffs)->mdef, rhmm->ciphone));

        /* HMMs for word-internal phones */
        prevhmm = NULL;
        for (p = 1; p < dict_pronlen(dict, wid) - 1; p++) {
            hmm = listelem_malloc(ffs->chan_alloc);
            hmm->ciphone = dict_pron(dict, wid, p);
            hmm->rc_id = -1;
            hmm->next = NULL;
            hmm_init(ffs->hmmctx, &hmm->hmm, FALSE,
                     dict2pid_internal(d2p,wid,p),
                     bin_mdef_pid2tmatid(ps_search_acmod(ffs)->mdef,
                                         hmm->ciphone));

            if (prevhmm)
                prevhmm->next = hmm;
            else
                rhmm->next = hmm;

            prevhmm = hmm;
        }

        /* Right-context phones */
        if (!dict_is_single_phone(dict, wid)) {
            hmm = fwdflat_search_alloc_all_rc(ffs, wid);
            /* Link in just allocated right-context phones */
            if (prevhmm)
                prevhmm->next = hmm;
            else
                rhmm->next = hmm;
        }
        ffs->word_chan[wid] = rhmm;
    }
}

static int
fwdflat_search_start(ps_search_t *base)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *)base;
    first_node_t *rhmm;
    int i;

    bptbl_reset(ffs->bptbl);
    ffs->oldest_bp = -1;
    for (i = 0; i < ps_search_n_words(ffs); i++)
        ffs->word_idx[i] = NO_BP;

    /* Start search with <s> */
    rhmm = (first_node_t *) ffs->word_chan[ps_search_start_wid(ffs)];
    hmm_enter(&rhmm->hmm, 0, NO_BP, 0);
    ffs->active_word_list[0][0] = ps_search_start_wid(ffs);
    ffs->n_active_word[0] = 1;

    ffs->best_score = 0;
    ffs->renormalized = FALSE;

    ffs->bptbl->n_frame = 0;
    ffs->st.n_fwdflat_chan = 0;
    ffs->st.n_fwdflat_words = 0;
    ffs->st.n_fwdflat_word_transition = 0;
    ffs->st.n_senone_active_utt = 0;

    return 0;
}

static int32
update_oldest_bp(fwdflat_search_t *ffs, hmm_t *hmm)
{
    int j;

    for (j = 0; j < hmm->n_emit_state; ++j)
        if (hmm_score(hmm, j) BETTER_THAN WORST_SCORE)
            if (hmm_history(hmm, j) < ffs->oldest_bp)
                ffs->oldest_bp = hmm_history(hmm, j);
    if (hmm_out_score(hmm) BETTER_THAN WORST_SCORE)
        if (hmm_out_history(hmm) < ffs->oldest_bp)
            ffs->oldest_bp = hmm_out_history(hmm);

    return ffs->oldest_bp;
}

static void
compute_fwdflat_sen_active(fwdflat_search_t *ffs, int frame_idx)
{
    int32 i, w;
    int32 *awl;
    first_node_t *rhmm;
    internal_node_t *hmm;

    acmod_clear_active(ps_search_acmod(ffs));
    ffs->oldest_bp = ffs->bptbl->n_ent;

    i = ffs->n_active_word[frame_idx & 0x1];
    awl = ffs->active_word_list[frame_idx & 0x1];

    for (w = *(awl++); i > 0; --i, w = *(awl++)) {
        rhmm = (first_node_t *)ffs->word_chan[w];
        if (hmm_frame(&rhmm->hmm) == frame_idx) {
            acmod_activate_hmm(ps_search_acmod(ffs), &rhmm->hmm);
            update_oldest_bp(ffs, &rhmm->hmm);
        }

        for (hmm = rhmm->next; hmm; hmm = hmm->next) {
            if (hmm_frame(&hmm->hmm) == frame_idx) {
                acmod_activate_hmm(ps_search_acmod(ffs), &hmm->hmm);
                update_oldest_bp(ffs, &hmm->hmm);
            }
        }
    }
    assert(ffs->oldest_bp < ffs->bptbl->n_ent);
}

static void
fwdflat_eval_chan(fwdflat_search_t *ffs, int frame_idx)
{
    int32 i, w, bestscore;
    int32 *awl;
    first_node_t *rhmm;
    internal_node_t *hmm;

    i = ffs->n_active_word[frame_idx & 0x1];
    awl = ffs->active_word_list[frame_idx & 0x1];
    bestscore = WORST_SCORE;

    ffs->st.n_fwdflat_words += i;

    /* Scan all active words. */
    for (w = *(awl++); i > 0; --i, w = *(awl++)) {
        rhmm = (first_node_t *) ffs->word_chan[w];
        if (hmm_frame(&rhmm->hmm) == frame_idx) {
            int32 score = chan_v_eval(rhmm);
            if ((score BETTER_THAN bestscore) && (w != ps_search_finish_wid(ffs)))
                bestscore = score;
            ffs->st.n_fwdflat_chan++;
        }

        for (hmm = rhmm->next; hmm; hmm = hmm->next) {
            if (hmm_frame(&hmm->hmm) == frame_idx) {
                int32 score = chan_v_eval(hmm);
                if (score BETTER_THAN bestscore)
                    bestscore = score;
                ffs->st.n_fwdflat_chan++;
            }
        }
    }

    ffs->best_score = bestscore;
}

static void
fwdflat_search_save_bp(fwdflat_search_t *ffs, int frame_idx,
                       int32 w, int32 score, int32 path, int32 rc)
{
    int32 bp;

    /* Look for an existing exit for this word in this frame. */
    bp = ffs->word_idx[w];
    if (bp != NO_BP) {
        bp_t *bpe = bptbl_ent(ffs->bptbl, bp);
        /* Keep only the best scoring one (this is a potential source
         * of search errors...) */
        if (bpe->score WORSE_THAN score) {
            if (bpe->bp != path) {
                bpe->bp = path;
                bptbl_fake_lmstate(ffs->bptbl, bp);
            }
            bpe->score = score;
        }
        /* But do keep track of scores for all right contexts, since
         * we need them to determine the starting path scores for any
         * successors of this word exit. */
        ffs->bptbl->bscore_stack[bpe->s_idx + rc] = score;
    }
    else {
        bp_t *bpe = bptbl_enter(ffs->bptbl, w, path, score, rc);
        ffs->word_idx[w] = bptbl_idx(ffs->bptbl, bpe);
        assert(frame_idx == bpe->frame);
    }
}

static void
fwdflat_prune_chan(fwdflat_search_t *ffs, int frame_idx)
{
    int32 i, cf, nf, w, pip, newscore, thresh, wordthresh;
    int32 *awl;
    first_node_t *rhmm;
    internal_node_t *hmm, *nexthmm;

    cf = frame_idx;
    nf = cf + 1;
    i = ffs->n_active_word[cf & 0x1];
    awl = ffs->active_word_list[cf & 0x1];
    bitvec_clear_all(ffs->word_active, ps_search_n_words(ffs));

    thresh = ffs->best_score + ffs->fwdflatbeam;
    wordthresh = ffs->best_score + ffs->fwdflatwbeam;
    pip = ffs->pip;

    /* Scan all active words. */
    for (w = *(awl++); i > 0; --i, w = *(awl++)) {
        rhmm = (first_node_t *) ffs->word_chan[w];
        /* Propagate active root channels */
        if (hmm_frame(&rhmm->hmm) == cf
            && hmm_bestscore(&rhmm->hmm) BETTER_THAN thresh) {
            hmm_frame(&rhmm->hmm) = nf;
            bitvec_set(ffs->word_active, w);

            /* Transitions out of root channel */
            newscore = hmm_out_score(&rhmm->hmm);
            if (rhmm->next) {
                assert(!dict_is_single_phone(ps_search_dict(ffs), w));
                newscore += pip;
                if (newscore BETTER_THAN thresh) {
                    hmm = rhmm->next;
                    /* Enter all right context phones */
                    if (hmm->rc_id >= 0) {
                        for (; hmm; hmm = hmm->next) {
                            if ((hmm_frame(&hmm->hmm) < cf)
                                || (newscore BETTER_THAN hmm_in_score(&hmm->hmm))) {
                                hmm_enter(&hmm->hmm, newscore,
                                          hmm_out_history(&rhmm->hmm), nf);
                            }
                        }
                    }
                    /* Just a normal word internal phone */
                    else {
                        if ((hmm_frame(&hmm->hmm) < cf)
                            || (newscore BETTER_THAN hmm_in_score(&hmm->hmm))) {
                                hmm_enter(&hmm->hmm, newscore,
                                          hmm_out_history(&rhmm->hmm), nf);
                        }
                    }
                }
            }
            else {
                assert(dict_is_single_phone(ps_search_dict(ffs), w));
                /* Word exit for single-phone words */
                if (newscore BETTER_THAN wordthresh) {
                    fwdflat_search_save_bp(ffs, cf, w, newscore,
                                           hmm_out_history(&rhmm->hmm), 0);
                }
            }
        }

        /* Transitions out of non-root channels. */
        for (hmm = rhmm->next; hmm; hmm = hmm->next) {
            if (hmm_frame(&hmm->hmm) >= cf) {
                /* Propagate forward HMMs inside the beam. */
                if (hmm_bestscore(&hmm->hmm) BETTER_THAN thresh) {
                    hmm_frame(&hmm->hmm) = nf;
                    bitvec_set(ffs->word_active, w);

                    newscore = hmm_out_score(&hmm->hmm);
                    /* Word-internal phones */
                    if (hmm->rc_id < 0) {
                        newscore += pip;
                        if (newscore BETTER_THAN thresh) {
                            nexthmm = hmm->next;
                            /* Enter all right-context phones. */
                            if (nexthmm->rc_id >= 0) {
                                 for (; nexthmm; nexthmm = nexthmm->next) {
                                    if ((hmm_frame(&nexthmm->hmm) < cf)
                                        || (newscore BETTER_THAN
                                            hmm_in_score(&nexthmm->hmm))) {
                                        hmm_enter(&nexthmm->hmm,
                                                  newscore,
                                                  hmm_out_history(&hmm->hmm),
                                                  nf);
                                    }
                                }
                            }
                            /* Enter single word-internal phone. */
                            else {
                                if ((hmm_frame(&nexthmm->hmm) < cf)
                                    || (newscore BETTER_THAN
                                        hmm_in_score(&nexthmm->hmm))) {
                                    hmm_enter(&nexthmm->hmm, newscore,
                                              hmm_out_history(&hmm->hmm), nf);
                                }
                            }
                        }
                    }
                    /* Right-context phones - apply word beam and exit. */
                    else {
                        if (newscore BETTER_THAN wordthresh) {
                            fwdflat_search_save_bp(ffs, cf, w, newscore,
                                                   hmm_out_history(&hmm->hmm),
                                                   hmm->rc_id);
                        }
                    }
                }
                /* Zero out inactive HMMs. */
                else if (hmm_frame(&hmm->hmm) != nf) {
                    hmm_clear_scores(&hmm->hmm);
                }
            }
        }
    }
}

static void
fwdflat_word_transition(fwdflat_search_t *ffs, int frame_idx)
{
    int32 cf, nf, b, thresh, pip, i, w, newscore;
    int32 best_silrc_score = 0, best_silrc_bp = 0;      /* FIXME: good defaults? */
    int32 *rcss;
    first_node_t *rhmm;
    int32 *awl;
    dict_t *dict = ps_search_dict(ffs);
    dict2pid_t *d2p = ps_search_dict2pid(ffs);

    cf = frame_idx;
    nf = cf + 1;
    thresh = ffs->best_score + ffs->fwdflatbeam;
    pip = ffs->pip;
    best_silrc_score = WORST_SCORE;

    /* Scan words exited in current frame */
    for (b = bptbl_ef_idx(ffs->bptbl, cf);
         b < bptbl_ef_idx(ffs->bptbl, cf + 1); b++) {
        xwdssid_t *rssid;
        int32 silscore;
        bp_t *ent;

        ent = bptbl_ent(ffs->bptbl, b);
        ffs->word_idx[ent->wid] = NO_BP;

        if (ent->wid == ps_search_finish_wid(ffs))
            continue;

        /* Get the mapping from right context phone ID to index in the
         * right context table and the bptbl->bscore_stack. */
        rcss = ffs->bptbl->bscore_stack + ent->s_idx;
        if (ent->last2_phone == -1)
            rssid = NULL;
        else
            rssid = dict2pid_rssid(d2p, ent->last_phone, ent->last2_phone);

        /* Transition to all successor words. */
        for (w = 0; w < ps_search_n_words(ffs); w++) {
            int32 n_used;

            if (!dict_real_word(ps_search_dict(ffs), w))
                continue;
            if (!ngram_model_set_known_wid(ffs->lmset,
                                           dict_basewid(ps_search_dict(ffs),w)))
                continue;

            /* Get the exit score we recorded in save_bwd_ptr(), or
             * something approximating it. */
            if (rssid)
                newscore = rcss[rssid->cimap[dict_first_phone(dict, w)]];
            else
                newscore = rcss[0];
            if (newscore == WORST_SCORE)
                continue;
            newscore += ngram_tg_score(ffs->lmset,
                                       dict_basewid(dict, w),
                                       ent->real_wid,
                                       ent->prev_real_wid, &n_used)>>SENSCR_SHIFT;
            newscore += pip >> SENSCR_SHIFT;

            /* Enter the next word */
            if (newscore BETTER_THAN thresh) {
                rhmm = (first_node_t *) ffs->word_chan[w];
                if ((hmm_frame(&rhmm->hmm) < cf)
                    || (newscore BETTER_THAN hmm_in_score(&rhmm->hmm))) {
                    hmm_enter(&rhmm->hmm, newscore, b, nf);
                    /* DICT2PID: This is where mpx ssids get introduced. */
                    /* Look up the ssid to use when entering this mpx triphone. */
                    hmm_mpx_ssid(&rhmm->hmm, 0) =
                        dict2pid_ldiph_lc(d2p, rhmm->ciphone, rhmm->ci2phone,
                                          dict_last_phone(dict, ent->wid));
                    assert(IS_S3SSID(hmm_mpx_ssid(&rhmm->hmm, 0)));
                    bitvec_set(ffs->word_active, w);
                }
            }
        }

        /* Get the best exit into silence. */
        if (rssid)
            silscore = rcss[rssid->cimap[ps_search_acmod(ffs)->mdef->sil]];
        else
            silscore = rcss[0];
        if (silscore BETTER_THAN best_silrc_score) {
            best_silrc_score = silscore;
            best_silrc_bp = b;
        }
    }

    /* Transition to <sil> */
    newscore = best_silrc_score + ffs->silpen + pip;
    if ((newscore BETTER_THAN thresh) && (newscore BETTER_THAN WORST_SCORE)) {
        w = ps_search_silence_wid(ffs);
        rhmm = (first_node_t *) ffs->word_chan[w];
        if ((hmm_frame(&rhmm->hmm) < cf)
            || (newscore BETTER_THAN hmm_in_score(&rhmm->hmm))) {
            hmm_enter(&rhmm->hmm, newscore,
                      best_silrc_bp, nf);
            bitvec_set(ffs->word_active, w);
        }
    }
    /* Transition to noise words (FIXME: Depends on dictionary organization) */
    newscore = best_silrc_score + ffs->fillpen + pip;
    if ((newscore BETTER_THAN thresh) && (newscore BETTER_THAN WORST_SCORE)) {
        for (w = ps_search_silence_wid(ffs) + 1; w < ps_search_n_words(ffs); w++) {
            rhmm = (first_node_t *) ffs->word_chan[w];
            /* Noise words that aren't a single phone will have NULL here. */
            if (rhmm == NULL)
                continue;
            if ((hmm_frame(&rhmm->hmm) < cf)
                || (newscore BETTER_THAN hmm_in_score(&rhmm->hmm))) {
                hmm_enter(&rhmm->hmm, newscore,
                          best_silrc_bp, nf);
                bitvec_set(ffs->word_active, w);
            }
        }
    }

    /* Reset initial channels of words that have become inactive even after word trans. */
    i = ffs->n_active_word[cf & 0x1];
    awl = ffs->active_word_list[cf & 0x1];
    for (w = *(awl++); i > 0; --i, w = *(awl++)) {
        rhmm = (first_node_t *) ffs->word_chan[w];
        if (hmm_frame(&rhmm->hmm) == cf) {
            hmm_clear_scores(&rhmm->hmm);
        }
    }
}

static void
fwdflat_renormalize_scores(fwdflat_search_t *ffs, int frame_idx, int32 norm)
{
    first_node_t *rhmm;
    internal_node_t *hmm;
    int32 i, cf, w, *awl;

    cf = frame_idx;

    /* Renormalize individual word channels */
    i = ffs->n_active_word[cf & 0x1];
    awl = ffs->active_word_list[cf & 0x1];
    for (w = *(awl++); i > 0; --i, w = *(awl++)) {
        rhmm = (first_node_t *) ffs->word_chan[w];
        if (hmm_frame(&rhmm->hmm) == cf) {
            hmm_normalize(&rhmm->hmm, norm);
        }
        for (hmm = rhmm->next; hmm; hmm = hmm->next) {
            if (hmm_frame(&hmm->hmm) == cf) {
                hmm_normalize(&hmm->hmm, norm);
            }
        }
    }

    ffs->renormalized = TRUE;
}

static int
fwdflat_search_step(ps_search_t *base, int frame_idx)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *)base;
    int16 const *senscr;
    int32 nf, i, j;
    int32 *nawl;
    int fi;

    /* Activate our HMMs for the current frame if need be. */
    if (!ps_search_acmod(ffs)->compallsen)
        compute_fwdflat_sen_active(ffs, frame_idx);

    /* Compute GMM scores for the current frame. */
    senscr = acmod_score(ps_search_acmod(ffs), &frame_idx);
    ffs->st.n_senone_active_utt += ps_search_acmod(ffs)->n_senone_active;

    /* Mark backpointer table for current frame. */
    fi = bptbl_push_frame(ffs->bptbl, ffs->oldest_bp);
    assert(fi == frame_idx);

    /* If the best score is equal to or worse than WORST_SCORE,
     * recognition has failed, don't bother to keep trying. */
    if (ffs->best_score == WORST_SCORE || ffs->best_score WORSE_THAN WORST_SCORE)
        return 0;
    /* Renormalize if necessary */
    if (ffs->best_score + (2 * ffs->fwdflatbeam) WORSE_THAN WORST_SCORE) {
        E_INFO("Renormalizing Scores at frame %d, best score %d\n",
               frame_idx, ffs->best_score);
        fwdflat_renormalize_scores(ffs, frame_idx, ffs->best_score);
    }

    ffs->best_score = WORST_SCORE;
    hmm_context_set_senscore(ffs->hmmctx, senscr);

    /* Evaluate HMMs */
    fwdflat_eval_chan(ffs, frame_idx);
    /* Prune HMMs and do phone transitions. */
    fwdflat_prune_chan(ffs, frame_idx);
    /* Do word transitions. */
    fwdflat_word_transition(ffs, frame_idx);

    /* Create next active word list */
    nf = frame_idx + 1;
    nawl = ffs->active_word_list[nf & 0x1];
    for (i = 0, j = 0; i < ps_search_n_words(ffs); i++) {
        if (bitvec_is_set(ffs->word_active, i)) {
            *(nawl++) = i;
            j++;
        }
    }
    ffs->n_active_word[nf & 0x1] = j;

    /* Return the number of frames processed. */
    return 1;
}

static int
fwdflat_search_finish(ps_search_t *base)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *)base;
    int32 cf;
    int fi;

    bitvec_clear_all(ffs->word_active, ps_search_n_words(ffs));

    /* This is the number of frames processed. */
    cf = ps_search_acmod(ffs)->output_frame;
    /* Add a mark in the backpointer table for one past the final frame. */
    fi = bptbl_push_frame(ffs->bptbl, ffs->oldest_bp);
    assert(fi == cf);

    /* Print out some statistics. */
    if (cf > 0) {
        E_INFO("%8d words recognized (%d/fr)\n",
               ffs->bptbl->n_ent, (ffs->bptbl->n_ent + (cf >> 1)) / (cf + 1));
        E_INFO("%8d senones evaluated (%d/fr)\n", ffs->st.n_senone_active_utt,
               (ffs->st.n_senone_active_utt + (cf >> 1)) / (cf + 1));
        E_INFO("%8d channels searched (%d/fr)\n",
               ffs->st.n_fwdflat_chan, ffs->st.n_fwdflat_chan / (cf + 1));
        E_INFO("%8d words searched (%d/fr)\n",
               ffs->st.n_fwdflat_words, ffs->st.n_fwdflat_words / (cf + 1));
        E_INFO("%8d word transitions (%d/fr)\n",
               ffs->st.n_fwdflat_word_transition,
               ffs->st.n_fwdflat_word_transition / (cf + 1));
    }

    ffs->done = TRUE;
    return 0;
}

static int32
fwdflat_search_prob(ps_search_t *base)
{
    /* FIXME: Going to estimate this from partial results in the future. */
    return 0;
}

static int
fwdflat_search_find_exit(fwdflat_search_t *ffs, int frame_idx, int32 *out_best_score)
{
    /* End of backpointers for this frame. */
    int end_bpidx;
    int best_exit, bp;
    int32 best_score;

    /* No hypothesis means no exit node! */
    if (ffs->bptbl->n_frame == 0)
        return NO_BP;

    if (frame_idx == -1 || frame_idx >= ffs->bptbl->n_frame)
        frame_idx = ffs->bptbl->n_frame - 1;
    end_bpidx = bptbl_ef_idx(ffs->bptbl, frame_idx);

    best_score = WORST_SCORE;
    best_exit = NO_BP;

    /* Scan back to find a frame with some backpointers in it. */
    while (frame_idx >= 0 && bptbl_ef_idx(ffs->bptbl, frame_idx) == end_bpidx)
        --frame_idx;
    /* This is NOT an error, it just means there is no hypothesis yet. */
    if (frame_idx < 0)
        return NO_BP;

    /* Now find the entry for </s> OR the best scoring entry. */
    for (bp = bptbl_ef_idx(ffs->bptbl, frame_idx); bp < end_bpidx; ++bp) {
        bp_t *bpe = bptbl_ent(ffs->bptbl, bp);
        if (bpe->wid == ps_search_finish_wid(ffs)
            || bpe->score BETTER_THAN best_score) {
            best_score = bpe->score;
            best_exit = bp;
        }
        if (bpe->wid == ps_search_finish_wid(ffs))
            break;
    }

    if (out_best_score) *out_best_score = best_score;
    return best_exit;
}

static char const *
fwdflat_search_hyp(ps_search_t *base, int32 *out_score)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *)base;
    char *c;
    size_t len;
    int bp, bpidx;

    bpidx = fwdflat_search_find_exit(ffs, -1, out_score);
    if (bpidx == NO_BP)
        return NULL;

    bp = bpidx;
    len = 0;
    while (bp != NO_BP) {
        bp_t *be = bptbl_ent(ffs->bptbl, bp);
        E_INFO("bp %d -> %d\n", bp, be->bp);
        assert(be->valid);
        bp = be->bp;
        if (dict_real_word(ps_search_dict(ffs), be->wid))
            len += strlen(dict_basestr(ps_search_dict(ffs), be->wid)) + 1;
    }

    ckd_free(base->hyp_str);
    if (len == 0) {
	base->hyp_str = NULL;
	return base->hyp_str;
    }
    base->hyp_str = ckd_calloc(1, len);

    bp = bpidx;
    c = base->hyp_str + len - 1;
    while (bp != NO_BP) {
        bp_t *be = bptbl_ent(ffs->bptbl, bp);
        size_t len;

        bp = be->bp;
        if (dict_real_word(ps_search_dict(ffs), be->wid)) {
            len = strlen(dict_basestr(ps_search_dict(ffs), be->wid));
            c -= len;
            memcpy(c, dict_basestr(ps_search_dict(ffs), be->wid), len);
            if (c > base->hyp_str) {
                --c;
                *c = ' ';
            }
        }
    }

    return base->hyp_str;
}

static int32
fwdflat_search_exit_score(fwdflat_search_t *ffs, bp_t *pbe, int rcphone)
{
    /* DICT2PID */
    /* Get the mapping from right context phone ID to index in the
     * right context table and the bptbl->bscore_stack. */
    /* FIXME: This function gets called like 50 zillion times, either
     * it should be inlined or we should find a better way to do
     * this. */
    E_DEBUG(99,("fwdflat_search_exit_score(%d,%d)\n", bptbl_idx(ffs->bptbl, pbe), rcphone));
    assert(pbe->valid);
    if (pbe->last2_phone == -1) {
        /* No right context for single phone predecessor words. */
        E_DEBUG(99,("last2_phone = %d s_idx = %d bscore = %d\n", -1,
                    pbe->s_idx, ffs->bptbl->bscore_stack[pbe->s_idx]));
        assert(ffs->bptbl->bscore_stack[pbe->s_idx] != WORST_SCORE);
        return ffs->bptbl->bscore_stack[pbe->s_idx];
    }
    else {
        xwdssid_t *rssid;
        /* Find the index for the last diphone of the previous word +
         * the first phone of the current word. */
        rssid = dict2pid_rssid(ps_search_dict2pid(ffs),
                               pbe->last_phone, pbe->last2_phone);
        E_DEBUG(99,("last2_phone = %d s_idx = %d rc = %d n_rc = %d bscore = %d\n",
                    pbe->last2_phone, pbe->s_idx, rssid->cimap[rcphone],
                    rssid->n_ssid,
                    ffs->bptbl->bscore_stack[pbe->s_idx + rssid->cimap[rcphone]]));
        /* This may be WORST_SCORE, which means that there was no exit
         * with rcphone as right context. */
        return ffs->bptbl->bscore_stack[pbe->s_idx + rssid->cimap[rcphone]];
    }
}

static void
fwdflat_search_bp2itor(ps_seg_t *seg, int bp)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *)seg->search;
    bp_t *be, *pbe;

    be = bptbl_ent(ffs->bptbl, bp);
    pbe = bptbl_ent(ffs->bptbl, be->bp);
    seg->word = dict_wordstr(ps_search_dict(ffs), be->wid);
    seg->ef = be->frame;
    seg->sf = pbe ? pbe->frame + 1 : 0;
    seg->prob = 0; /* Bogus value... */
    /* Compute acoustic and LM scores for this segment. */
    if (pbe == NULL) {
        seg->ascr = be->score;
        seg->lscr = 0;
        seg->lback = 0;
    }
    else {
        int32 start_score;

        /* Find ending path score of previous word. */
        start_score = fwdflat_search_exit_score(ffs, pbe,
                                                dict_first_phone(ps_search_dict(ffs),
                                                                 be->wid));
        assert(start_score BETTER_THAN WORST_SCORE);
        if (be->wid == ps_search_silence_wid(ffs)) {
            /* FIXME: Nasty action at a distance here to deal with the
             * silence length limiting stuff in fwdflat_search_fwdflat.c */
            if (dict_first_phone(ps_search_dict(ffs), be->wid)
                == ps_search_acmod(ffs)->mdef->sil)
                seg->lscr = 0;
            else
                seg->lscr = ffs->silpen;
        }
        else if (dict_filler_word(ps_search_dict(ffs), be->wid)) {
            seg->lscr = ffs->fillpen;
        }
        else {
            seg->lscr = ngram_tg_score(ffs->lmset,
                                       be->real_wid,
                                       pbe->real_wid,
                                       pbe->prev_real_wid, &seg->lback)>>SENSCR_SHIFT;
            seg->lscr = (int32)(seg->lscr * seg->lwf);
        }
        seg->ascr = be->score - start_score - seg->lscr;
    }
}

static void
ngram_bp_seg_free(ps_seg_t *seg)
{
    bptbl_seg_t *itor = (bptbl_seg_t *)seg;
    
    ckd_free(itor->bpidx);
    ckd_free(itor);
}

static ps_seg_t *
ngram_bp_seg_next(ps_seg_t *seg)
{
    bptbl_seg_t *itor = (bptbl_seg_t *)seg;

    if (++itor->cur == itor->n_bpidx) {
        ngram_bp_seg_free(seg);
        return NULL;
    }

    fwdflat_search_bp2itor(seg, itor->bpidx[itor->cur]);
    return seg;
}

static ps_segfuncs_t ngram_bp_segfuncs = {
    /* seg_next */ ngram_bp_seg_next,
    /* seg_free */ ngram_bp_seg_free
};

static ps_seg_t *
fwdflat_search_bp_iter(fwdflat_search_t *ffs, int bpidx, float32 lwf)
{
    bptbl_seg_t *itor;
    int bp, cur;

    /* Calling this an "iterator" is a bit of a misnomer since we have
     * to get the entire backtrace in order to produce it.  On the
     * other hand, all we actually need is the bptbl IDs, and we can
     * allocate a fixed-size array of them. */
    itor = ckd_calloc(1, sizeof(*itor));
    itor->base.vt = &ngram_bp_segfuncs;
    itor->base.search = ps_search_base(ffs);
    itor->base.lwf = lwf;
    itor->n_bpidx = 0;
    bp = bpidx;
    while (bp != NO_BP) {
        bp_t *be = bptbl_ent(ffs->bptbl, bp);
        bp = be->bp;
        ++itor->n_bpidx;
    }
    if (itor->n_bpidx == 0) {
        ckd_free(itor);
        return NULL;
    }
    itor->bpidx = ckd_calloc(itor->n_bpidx, sizeof(*itor->bpidx));
    cur = itor->n_bpidx - 1;
    bp = bpidx;
    while (bp != NO_BP) {
        bp_t *be = bptbl_ent(ffs->bptbl, bp);
        itor->bpidx[cur] = bp;
        bp = be->bp;
        --cur;
    }

    /* Fill in relevant fields for first element. */
    fwdflat_search_bp2itor((ps_seg_t *)itor, itor->bpidx[0]);

    return (ps_seg_t *)itor;
}

static ps_seg_t *
fwdflat_search_seg_iter(ps_search_t *search, int32 *out_score)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *)search;
    int32 bpidx;

    bpidx = fwdflat_search_find_exit(ffs, -1, out_score);
    return fwdflat_search_bp_iter(ffs, bpidx, 1.0);

    return NULL;
}
