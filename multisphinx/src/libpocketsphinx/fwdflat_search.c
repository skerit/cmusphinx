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
 * @file fwdflat_search.c Flat lexicon search.
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
static int fwdflat_search_decode(ps_search_t *base);
static int fwdflat_search_finish(ps_search_t *base);
static int fwdflat_search_free(ps_search_t *base);
static char const *fwdflat_search_hyp(ps_search_t *base, int32 *out_score);
static int32 fwdflat_search_prob(ps_search_t *base);
static ps_seg_t *fwdflat_search_seg_iter(ps_search_t *base, int32 *out_score);

static ps_searchfuncs_t fwdflat_funcs = {
    /* name: */   "fwdflat",
    /* free: */   fwdflat_search_free,
    /* decode: */ fwdflat_search_decode,
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
                    dict_t *dict, dict2pid_t *d2p,
                    bptbl_t *input_bptbl)
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
    ffs->word_chan = ckd_calloc(ps_search_n_words(ffs),
                                sizeof(*ffs->word_chan));
    E_INFO("Allocated %d KiB for word HMMs\n",
           (int)ps_search_n_words(ffs) * sizeof(*ffs->word_chan) / 1024);
    ffs->word_active = bitvec_alloc(ps_search_n_words(ffs));
    ffs->word_idx = ckd_calloc(ps_search_n_words(ffs),
                               sizeof(*ffs->word_idx));
    ffs->input_words = ckd_calloc(ps_search_n_words(ffs), sizeof(*ffs->input_words));
    E_INFO("Allocated %d KiB for active word flags\n",
           (int)ps_search_n_words(ffs) * sizeof(*ffs->input_words) / 1024);
    ffs->input_arcs = fwdflat_arc_buffer_init();

    ffs->bptbl = bptbl_init(d2p, cmd_ln_int32_r(config, "-latsize"), 256);
    if (input_bptbl)
        ffs->input_bptbl = bptbl_retain(input_bptbl);

    /* Allocate active word list array */
    ffs->active_word_list = ckd_calloc_2d(2, ps_search_n_words(ffs),
                                          sizeof(**ffs->active_word_list));
    E_INFO("Allocated %d KiB for active word list\n",
           (ps_search_n_words(ffs) * sizeof(**ffs->active_word_list)
            + ps_search_n_words(ffs) * sizeof(*ffs->active_word_list)) / 1024);

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

static int
fwdflat_search_free(ps_search_t *base)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *)base;

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
    fwdflat_arc_buffer_reset(ffs->input_arcs);
    ffs->oldest_bp = -1;
    ffs->next_idx = 0;
    for (i = 0; i < ps_search_n_words(ffs); i++)
        ffs->word_idx[i] = NO_BP;

    /* Start search with <s> */
    rhmm = (first_node_t *) ffs->word_chan[ps_search_start_wid(ffs)];
    hmm_enter(&rhmm->hmm, 0, NO_BP, 0);
    ffs->active_word_list[0][0] = ps_search_start_wid(ffs);
    ffs->n_active_word[0] = 1;

    memset(ffs->input_words, 0, ps_search_n_words(ffs) * sizeof(*ffs->input_words));
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
    ffs->oldest_bp = bptbl_end_idx(ffs->bptbl);

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
    assert(ffs->oldest_bp < bptbl_end_idx(ffs->bptbl));
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
        bptbl_set_rcscore(ffs->bptbl, bpe, rc, score);
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
        rcss = bptbl_rcscores(ffs->bptbl, ent);
        if (ent->last2_phone == -1)
            rssid = NULL;
        else
            rssid = dict2pid_rssid(d2p, ent->last_phone, ent->last2_phone);

        /* Transition to all successor words. */
        for (w = 0; w < ps_search_n_words(ffs); w++) {
            int32 n_used;

            if (!ngram_model_set_known_wid(ffs->lmset, dict_basewid(dict,w)))
                continue;
            if (ffs->input_bptbl != NULL
                && ffs->input_words[w] < ffs->min_ef_width)
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
fwdflat_search_one_frame(fwdflat_search_t *ffs, int frame_idx)
{
    acmod_t *acmod = ps_search_acmod(ffs);
    int16 const *senscr;
    int32 nf, i, j;
    int32 *nawl;
    int fi;

    printf("Searching frame %d\n", frame_idx);
    /* Activate our HMMs for the current frame if need be. */
    if (!acmod->compallsen)
        compute_fwdflat_sen_active(ffs, frame_idx);

    /* Compute GMM scores for the current frame. */
    senscr = acmod_score(acmod, frame_idx);
    ffs->st.n_senone_active_utt += acmod->n_senone_active;

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

    /* Release the frame just searched. */
    acmod_release(acmod, frame_idx);

    /* Return the number of frames processed. */
    return 1;
}

static void
fwdflat_dump_active_words(fwdflat_search_t *ffs, int sf, int ef)
{
    int i, j;

    E_INFO("Active words in %d:%d:\n", sf, ef);
    for (i = 0; i < ps_search_n_words(ffs); ++i) {
        if (ffs->input_words[i] >= ffs->min_ef_width) {
            E_INFO_NOFN("%s (%d exits)\n",
                        dict_wordstr(ps_search_dict(ffs), i),
                        ffs->input_words[i]);
            ++j;
        }
    }
    E_INFO("%d active words\n", j);
}

static int
fwdflat_search_expand_arcs(fwdflat_search_t *ffs, int sf, int ef)
{
    fwdflat_arc_t *arc_start, *arc_end, *arc;
    dict_t *dict = ps_search_dict(ffs);

    arc_start = fwdflat_arc_buffer_iter(ffs->input_arcs, sf);
    arc_end = fwdflat_arc_buffer_iter(ffs->input_arcs, ef);
    memset(ffs->input_words, 0, ps_search_n_words(ffs) * sizeof(*ffs->input_words));
    for (arc = arc_start; arc != arc_end;
         arc = fwdflat_arc_next(ffs->input_arcs, arc)) {
        /* As before, only count things to which transitions can be
         * made in the language model. */
        if (!(dict_filler_word(dict, arc->wid) || arc->wid == dict->startwid))
            ++ffs->input_words[arc->wid];
    }
    if (0)
        fwdflat_dump_active_words(ffs, sf, ef);
    return 0;
}

static int
fwdflat_search_decode_2ndpass(fwdflat_search_t *ffs, acmod_t *acmod)
{
    int next_sf; /**< First frame pointed to by active bps. */
    int frame_idx, final;

    fwdflat_search_start(ps_search_base(ffs));
    frame_idx = 0;
    final = FALSE;
    /* Wait for things to get retired from the input bptbl. */
    while (!final) {
        int k, timeout, end_win;

        /* next_sf is the first starting frame that is still
         * referenced by active backpointers. */
        if (bptbl_wait(ffs->input_bptbl, -1) < 0)
            break;
        next_sf = bptbl_active_sf(ffs->input_bptbl);
        /* Extend the arc buffer the appropriate number of frames. */
        if (fwdflat_arc_buffer_extend(ffs->input_arcs, next_sf) > 0) {
            E_INFO("oldest_bp %d next_idx %d\n",
                   ffs->input_bptbl->oldest_bp, ffs->next_idx);
            /* Add the next chunk of bps to the arc buffer. */
            ffs->next_idx = fwdflat_arc_buffer_add_bps
                (ffs->input_arcs, ffs->input_bptbl,
                 ffs->next_idx, bptbl_retired_idx(ffs->input_bptbl));
            fwdflat_arc_buffer_commit(ffs->input_arcs);
            /* Release bps we won't need anymore. */
            E_INFO("oldest_bp %d next_idx %d\n",
                   ffs->input_bptbl->oldest_bp, ffs->next_idx);
            bptbl_release(ffs->input_bptbl, ffs->input_bptbl->oldest_bp);
        }
        /* We do something different depending on whether the input
         * bptbl has been finalized or not.  If it has, we run out the
         * clock and finish, otherwise we only search forward as far
         * as the arc buffer goes. */
        end_win = frame_idx + ffs->max_sf_win;
        final = (bptbl_active_frame(ffs->input_bptbl)
                 == bptbl_frame_idx(ffs->input_bptbl));
        /* If we don't have enough of a window then keep waiting. */
        if (!final /* Don't care how big it is if we're final. */
            && fwdflat_arc_buffer_iter(ffs->input_arcs, end_win - 1) == NULL)
            continue;
        /* We are going to use the window so truncate it. */
        if (end_win > bptbl_frame_idx(ffs->input_bptbl))
            end_win = bptbl_frame_idx(ffs->input_bptbl);

        /* Whether we are final or not determines whether we wait for
         * the acmod or not. */
        if (final)
            timeout = -1; /* Run until end of utterance. */
        else
            timeout = 0;  /* Don't wait for results, we will block on
                           * input_bptbl instead. */
        while ((frame_idx = acmod_wait(acmod, timeout)) >= 0) {
            E_INFO("Searching frame %d window end %d\n",
                   frame_idx, end_win);
            fwdflat_search_expand_arcs(ffs, frame_idx, end_win);
            if ((k = fwdflat_search_one_frame(ffs, frame_idx)) <= 0)
                break;
            frame_idx += k;
            fwdflat_arc_buffer_release(ffs->input_arcs, frame_idx);
        }
    }
    fwdflat_search_finish(ps_search_base(ffs));
    return frame_idx;
}

static int
fwdflat_search_decode(ps_search_t *base)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *)base;
    acmod_t *acmod = ps_search_acmod(base);
    int frame_idx, nfr, k;

    E_INFO("WTF decoding\n");
    if (ffs->input_bptbl)
        return fwdflat_search_decode_2ndpass(ffs, ps_search_acmod(base));

    nfr = 0;
    fwdflat_search_start(base);
    while ((frame_idx = acmod_wait(acmod, -1)) >= 0) {
        if ((k = fwdflat_search_one_frame(ffs, frame_idx)) <= 0)
            break;
        nfr += k;
    }
    fwdflat_search_finish(base);
    bptbl_dump(ffs->bptbl);
    /* This means we were canceled.  FIXME: Not clear whether calling
     * fwdflat_search_finish() was necessary in this case. */
    if (k < 0)
        return k;
    return nfr;
}

static int
fwdflat_search_finish(ps_search_t *base)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *)base;
    acmod_t *acmod = ps_search_acmod(ffs);
    int cf;
    
    bitvec_clear_all(ffs->word_active, ps_search_n_words(ffs));

    /* This is the number of frames of input. */
    cf = acmod->output_frame;

    /* Finalize the backpointer table. */
    bptbl_finalize(ffs->bptbl);

    /* Print out some statistics. */
    if (cf > 0) {
        E_INFO("%8d words recognized (%d/fr)\n",
               bptbl_end_idx(ffs->bptbl), (bptbl_end_idx(ffs->bptbl) + (cf >> 1)) / (cf + 1));
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
    E_INFO("Allocated %d arcs and %d start frames\n",
           garray_alloc_size(ffs->input_arcs->arcs),
           garray_alloc_size(ffs->input_arcs->sf_idx));

    return 0;
}

static int32
fwdflat_search_prob(ps_search_t *base)
{
    /* FIXME: Going to estimate this from partial results in the future. */
    return 0;
}

static char const *
fwdflat_search_hyp(ps_search_t *base, int32 *out_score)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *)base;

    ckd_free(base->hyp_str);
    base->hyp_str = bptbl_hyp(ffs->bptbl, out_score, ps_search_finish_wid(ffs));
    return base->hyp_str;
}

static ps_seg_t *
fwdflat_search_seg_iter(ps_search_t *base, int32 *out_score)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *)base;
    return bptbl_seg_iter(ffs->bptbl, out_score, ps_search_finish_wid(ffs));
}

fwdflat_arc_buffer_t *
fwdflat_arc_buffer_init(void)
{
    fwdflat_arc_buffer_t *fab;

    fab = ckd_calloc(1, sizeof(*fab));
    fab->refcount = 1;
    fab->arcs = garray_init(0, sizeof(fwdflat_arc_t));
    fab->sf_idx = garray_init(0, sizeof(int));

    return fab;
}

fwdflat_arc_buffer_t *
fwdflat_arc_buffer_retain(fwdflat_arc_buffer_t *fab)
{
    ++fab->refcount;
    return fab;
}

int
fwdflat_arc_buffer_free(fwdflat_arc_buffer_t *fab)
{
    if (fab == NULL)
        return 0;
    if (--fab->refcount > 0)
        return fab->refcount;

    garray_free(fab->sf_idx);
    garray_free(fab->arcs);
    ckd_free(fab);
    return 0;
}

void
fwdflat_arc_buffer_dump(fwdflat_arc_buffer_t *fab)
{
    size_t i, n_arcs;

    n_arcs = garray_next_idx(fab->arcs);
    E_INFO("Arc buffer %p: %d arcs:\n", fab, n_arcs);
    for (i = garray_base(fab->arcs); i < n_arcs; ++i) {
        fwdflat_arc_t *arc = garray_ptr(fab->arcs, fwdflat_arc_t, i);
        E_INFO_NOFN("wid %d sf %d ef %d\n",
                    arc->wid, arc->sf, arc->ef);
    }
}

int
fwdflat_arc_buffer_extend(fwdflat_arc_buffer_t *fab, int next_sf)
{
    if (next_sf == fab->next_sf)
        return 0;
    garray_expand_to(fab->sf_idx, next_sf);
    fab->next_sf = next_sf;
    garray_clear(fab->sf_idx, fab->active_sf, fab->next_sf - fab->active_sf);
    return next_sf - fab->active_sf;
}

bpidx_t
fwdflat_arc_buffer_add_bps(fwdflat_arc_buffer_t *fab,
                           bptbl_t *bptbl, bpidx_t start,
                           bpidx_t end)
{
    bpidx_t idx, next_idx;
    int n_arcs;

    n_arcs = 0;
    next_idx = -1;
    for (idx = start; idx < end; ++idx) {
        fwdflat_arc_t arc;
        bp_t *ent;

        /* Convert it to an arc. */
        ent = bptbl_ent(bptbl, idx);
        arc.wid = ent->wid;
        arc.sf = bptbl_sf(bptbl, idx);
        arc.ef = ent->frame;
        /* If it's inside the appropriate frame span, add it. */
        if (arc.sf >= fab->active_sf && arc.sf < fab->next_sf) {
            garray_append(fab->arcs, &arc);
            /* Increment the frame counter for its start frame. */
            ++garray_ent(fab->sf_idx, int, arc.sf);
            ++n_arcs;
        }
        else {
            if (arc.sf >= fab->active_sf && next_idx == -1)
                next_idx = idx;
        }
    }

    E_INFO("Added bps from frame %d to %d, index %d to %d\n",
           fab->active_sf, fab->next_sf,
           start, end);
    if (next_idx == -1)
        next_idx = end;
    return next_idx;
}

int
fwdflat_arc_buffer_commit(fwdflat_arc_buffer_t *fab)
{
    size_t n_arcs;
    int n_active_fr, i, prev_count;
    garray_t *active_arc;
    garray_t *active_sf;
    int *sf;

    /* Save frame and arc counts. */
    n_active_fr = fab->next_sf - fab->active_sf;
    n_arcs = garray_next_idx(fab->arcs) - fab->active_arc;

    /* Nothing to do... */
    if (n_active_fr == 0) {
        assert(n_arcs == 0);
        return 0;
    }

#if 0
    E_INFO("sf_idx before (%d:%d):", fab->active_sf, fab->next_sf);
    for (i = fab->active_sf; i < fab->next_sf; ++i)
        E_INFOCONT(" %d", garray_ent(fab->sf_idx, int, i));
    E_INFOCONT("\n");
#endif

    /* Sum forward frame counters to create arc indices. */
    sf = garray_ptr(fab->sf_idx, int, fab->active_sf);
    prev_count = sf[0];
    sf[0] = fab->active_arc;
    for (i = 1; i < n_active_fr; ++i)  {
        int tmp = sf[i];
        sf[i] = sf[i-1] + prev_count;
        prev_count = tmp;
    }

    if (n_arcs > 0) {
        /* Permute incoming arcs to match frame counters */
        active_sf = garray_slice(fab->sf_idx, fab->active_sf, n_active_fr);
        active_arc = garray_slice(fab->arcs, fab->active_arc, n_arcs);

        for (i = 0; i < n_arcs; ++i) {
            fwdflat_arc_t *arc = garray_ptr(active_arc, fwdflat_arc_t, i);
            int *pos = garray_ptr(active_sf, int, arc->sf - fab->active_sf);
            /* Copy it into place. */
            garray_ent(fab->arcs, fwdflat_arc_t, *pos) = *arc;
            /* Increment local frame counter. */
            *pos += 1;
        }
#if 0
        E_INFO("sf_idx after (%d:%d):", fab->active_sf, fab->next_sf);
        for (i = fab->active_sf; i < fab->next_sf; ++i) {
            int idx = garray_ent(fab->sf_idx, int, i);
            E_INFOCONT(" %d=%d:%d",
                       idx, garray_ent(fab->arcs, fwdflat_arc_t, idx).sf,
                       garray_ent(fab->arcs, fwdflat_arc_t, idx).ef);
        }
        E_INFOCONT("\n");
#endif

        garray_free(active_sf);
        garray_free(active_arc);
    }

    /* Update frame and arc pointers. */
    fab->active_sf += n_active_fr;
    fab->active_arc += n_arcs;
#if 0
    E_INFO("active_sf => %d active_arc => %d\n",
           fab->active_sf, fab->active_arc);
    E_INFO("active_arc-1 -> %d -> %d,%d,%d\n",
           fab->active_arc-1,
           garray_ent(fab->arcs, fwdflat_arc_t, fab->active_arc-1).wid,
           garray_ent(fab->arcs, fwdflat_arc_t, fab->active_arc-1).sf,
           garray_ent(fab->arcs, fwdflat_arc_t, fab->active_arc-1).ef
        );
#endif
    return n_arcs;
}

fwdflat_arc_t *
fwdflat_arc_buffer_iter(fwdflat_arc_buffer_t *fab, int sf)
{
    int idx;
    if (sf < garray_base(fab->sf_idx) || sf >= fab->active_sf)
        return NULL;
    idx = garray_ent(fab->sf_idx, int, sf);
    if (idx >= fab->active_arc)
        return NULL;
    return garray_ptr(fab->arcs, fwdflat_arc_t, idx);

}

fwdflat_arc_t *
fwdflat_arc_next(fwdflat_arc_buffer_t *fab, fwdflat_arc_t *ab)
{
    ab += 1;
    if (ab >= garray_ptr(fab->arcs, fwdflat_arc_t, fab->active_arc))
        return NULL;
    return ab;
}

fwdflat_arc_t *
fwdflat_arc_buffer_wait(fwdflat_arc_buffer_t *fab, int sf)
{
    /* FIXME: Implement this... */
    return NULL;
}

int
fwdflat_arc_buffer_release(fwdflat_arc_buffer_t *fab, int first_sf)
{
    int next_first_arc;
    if (first_sf == garray_base(fab->sf_idx))
        return 0;

    /* Get the new first arc. */
    next_first_arc = garray_ent(fab->sf_idx, int, first_sf);
    /* Shift back start frames and arcs. */
#if 0
    E_INFO("Shifting back %d entries first_sf -> %d -> %d\n",
           first_sf - garray_base(fab->sf_idx),
           first_sf, garray_ent(fab->sf_idx, int, first_sf));
#endif
    garray_shift_from(fab->sf_idx, first_sf);
    garray_set_base(fab->sf_idx, first_sf);
#if 0
    E_INFO("Shifting back %d entries first_arc -> %d -> %d,%d,%d\n",
           next_first_arc - garray_base(fab->arcs), next_first_arc,
           garray_ent(fab->arcs, fwdflat_arc_t, next_first_arc).wid,
           garray_ent(fab->arcs, fwdflat_arc_t, next_first_arc).sf,
           garray_ent(fab->arcs, fwdflat_arc_t, next_first_arc).ef
        );
#endif
    garray_shift_from(fab->arcs, next_first_arc);
    garray_set_base(fab->arcs, next_first_arc);

    return 0;
}

void
fwdflat_arc_buffer_reset(fwdflat_arc_buffer_t *fab)
{
    fab->active_sf = fab->next_sf = 0;
    fab->active_arc = 0;
    garray_reset(fab->arcs);
    garray_reset(fab->sf_idx);
}
