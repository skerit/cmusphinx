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
#include <sphinxbase/pio.h>
#include <sphinxbase/err.h>

/* Local headers. */
#include "search_internal.h"
#include "fwdflat_search.h"
#include "tmat.h"
#include "hmm.h"

/* Turn this on to dump channels for debugging */
#define __CHAN_DUMP__		0
#if __CHAN_DUMP__
#define chan_v_eval(chan) hmm_dump_vit_eval(&(chan)->hmm, stderr)
#else
#define chan_v_eval(chan) hmm_vit_eval(&(chan)->hmm)
#endif

#if 0
#undef E_DEBUG
#define E_DEBUG(level,x) E_INFO x
#undef E_DEBUGCONT
#define E_DEBUGCONT(level,x) E_INFOCONT x
#endif

/**
 * Phone HMM data type.
 *
 * Not the first HMM for words, which multiplex HMMs based on
 * different left contexts.  This structure is used both in the
 * dynamic HMM tree structure and in the per-word last-phone right
 * context fanout.
 */
typedef struct internal_node_s {
    hmm_t hmm; /**< Basic HMM structure.  This *must*
     be first in the structure because
     internal_node_t and first_node_t are
     sometimes used interchangeably */
    struct internal_node_s *next;/**< first descendant of this channel; or, in the
     case of the last phone of a word, the next
     alternative right context channel */
    int16 ciphone; /**< ciphone for this node */
    int16 rc_id; /**< right-context id for last phone of words */
} internal_node_t;

/**
 * Lexical tree node data type for the first phone (first) of each word HMM.
 *
 * Each state may have a different parent static HMM.  Most fields are
 * similar to those in internal_node_t.
 */
typedef struct first_node_s {
    hmm_t hmm; /**< Basic HMM structure.  This *must* be first in
     the structure because internal_node_t and first_node_t are
     sometimes used interchangeably. */
    internal_node_t *next; /**< first descendant of this channel */

    int16 ciphone; /**< first ciphone of this node; all words firsted at this
     node begin with this ciphone */
    int16 ci2phone; /**< second ciphone of this node; one first HMM for each
     unique right context */
} first_node_t;

/**
 * Various statistics for profiling.
 */
typedef struct fwdflat_stats_s {
    int32 n_fwdflat_chan;
    int32 n_fwdflat_words;
    int32 n_fwdflat_word_transition;
    int32 n_senone_active_utt;
} fwdflat_stats_t;

/**
 * Word loop-based forward search.
 */
typedef struct fwdflat_search_s {
    search_t base;
    ngram_model_t *lmset;
    hmm_context_t *hmmctx;

    listelem_alloc_t *chan_alloc; /**< For internal_node_t */
    listelem_alloc_t *root_chan_alloc; /**< For first_node_t */

    /**
     * Backpointer table (temporary storage for active word arcs).
     */
    bptbl_t *bptbl;
    /* FIXME: This may be confusing with bptbl->oldest_bp (but maybe not) */
    int32 oldest_bp; /**< Oldest bptable entry active in decoding graph. */

    int32 *word_idx; /**< BPTable index for any word in current frame;
     cleared before each frame */
    int32 *rcss; /**< Temporary storage for right context scores. */

    /**
     * Vocabulary expansion map.
     */
    vocab_map_t *vmap;

    /**
     * Cumulative vocabulary for current utterance.
     */
    bitvec_t *utt_vocab;
    garray_t *word_list;
    int max_sf_win; /**< Window size for word entries. */

    bitvec_t *expand_words;
    int32 *expand_word_list;
    int n_expand_word;

    /**
     * First HMMs for multiple-phone words.
     */
    first_node_t **word_chan;
    bitvec_t *word_active; /**< array of active flags for all words. */

    /**
     * Array of active words for current and next frame.
     *
     * Similarly to active_chan_list, active_word_list[f mod 2] = list
     * of word ids for which active channels exist in word_chan in
     * frame f.
     */
    int32 **active_word_list;
    int32 n_active_word[2]; /**< Number entries in active_word_list */

    int32 best_score; /**< Best Viterbi path score. */
    int32 renormalized; /**< renormalized? (FIXME: allow multiple renorms) */

    fwdflat_stats_t st; /**< Various statistics for profiling. */

    /**
     * Best word exit in the most recent frame (used for partial results).
     */
    bpidx_t best_exit;
    /**
     `     * Final word in best path (used for partial results).
     */
    int32 best_exit_wid;

    /* A children's treasury of beam widths. */
    int32 fwdflatbeam;
    int32 fwdflatwbeam;
    int32 fillpen;
    int32 silpen;
    int32 pip;
    int32 lw;
} fwdflat_search_t;

static search_t *
fwdflat_search_init(search_t *other, cmd_ln_t *config, acmod_t *acmod,
        dict2pid_t *d2p);
static int fwdflat_search_start(search_t *base);
static int fwdflat_search_decode(search_t *base);
static int fwdflat_search_finish(search_t *base);
static int fwdflat_search_free(search_t *base);
static char const *fwdflat_search_hyp(search_t *base, int32 *out_score);
static int32 fwdflat_search_prob(search_t *base);
static seg_iter_t *fwdflat_search_seg_iter(search_t *base, int32 *out_score);
static bptbl_t *fwdflat_search_bptbl(search_t *base);
static ngram_model_t *fwdflat_search_lmset(search_t *base);

static searchfuncs_t fwdflat_funcs = {
/* name: */"fwdflat",
/* init: */fwdflat_search_init,
/* free: */fwdflat_search_free,
/* decode: */fwdflat_search_decode,
/* hyp: */fwdflat_search_hyp,
/* prob: */fwdflat_search_prob,
/* seg_iter: */fwdflat_search_seg_iter,
/* bptbl: */fwdflat_search_bptbl,
/* lmset: */fwdflat_search_lmset, };

static void build_fwdflat_word_chan(fwdflat_search_t *ffs, int32 wid);

static void fwdflat_search_calc_beams(fwdflat_search_t *ffs)
{
    cmd_ln_t *config;
    acmod_t *acmod;

    config = search_config((search_t *)ffs);
    acmod = search_acmod(ffs);

    /* Log beam widths. */
    ffs->fwdflatbeam = logmath_log(acmod->lmath,
            cmd_ln_float64_r(config, "-fwdflatbeam")) >> SENSCR_SHIFT;
    ffs->fwdflatwbeam = logmath_log(acmod->lmath,
            cmd_ln_float64_r(config, "-fwdflatwbeam")) >> SENSCR_SHIFT;

    /* Other things. */
    ffs->pip = logmath_log(acmod->lmath, cmd_ln_float32_r(config, "-pip"))
            >> SENSCR_SHIFT;
    ffs->silpen = logmath_log(acmod->lmath,
            cmd_ln_float32_r(config, "-silprob")) >> SENSCR_SHIFT;
    ffs->fillpen = logmath_log(acmod->lmath,
            cmd_ln_float32_r(config, "-fillprob")) >> SENSCR_SHIFT;
    ffs->max_sf_win = cmd_ln_int32_r(search_config((search_t *)ffs), "-fwdflatsfwin");
}

static void fwdflat_search_update_widmap(fwdflat_search_t *ffs)
{
    const char **words;
    int32 i, n_words;

    /* It's okay to include fillers since they won't be in the LM */
    n_words = search_n_words(ffs);
    words = ckd_calloc(n_words, sizeof(*words));
    /* This will include alternates, again, that's okay since they aren't in the LM */
    for (i = 0; i < n_words; ++i)
        words[i] = (const char *) dict_wordstr(search_dict(ffs), i);
    ngram_model_set_map_words(ffs->lmset, words, n_words);
    ckd_free(words);
}

searchfuncs_t const *
fwdflat_search_query(void)
{
    return &fwdflat_funcs;
}

static search_t *
fwdflat_search_init(search_t *other, cmd_ln_t *config, acmod_t *acmod,
        dict2pid_t *d2p)
{
    fwdflat_search_t *ffs;
    const char *path;

    ffs = ckd_calloc(1, sizeof(*ffs));
    search_base_init(&ffs->base, &fwdflat_funcs, config, acmod, d2p);
    ffs->hmmctx = hmm_context_init(bin_mdef_n_emit_state(acmod->mdef),
            acmod->tmat->tp, NULL, acmod->mdef->sseq);
    if (ffs->hmmctx == NULL) {
        search_free(search_base(ffs));
        return NULL;
    }
    ffs->chan_alloc = listelem_alloc_init(sizeof(internal_node_t));
    ffs->root_chan_alloc = listelem_alloc_init(sizeof(first_node_t));
    ffs->rcss = ckd_calloc(bin_mdef_n_ciphone(acmod->mdef),
            sizeof(*ffs->rcss));

    /* Calculate various beam widths and such. */
    fwdflat_search_calc_beams(ffs);

    /* Allocate a billion different tables for stuff. */
    ffs->word_chan = ckd_calloc(search_n_words(ffs),
            sizeof(*ffs->word_chan));
    E_INFO("Allocated %d KiB for word HMMs\n",
    (int)search_n_words(ffs) * sizeof(*ffs->word_chan) / 1024);
    ffs->word_active = bitvec_alloc(search_n_words(ffs));
    ffs->word_idx = ckd_calloc(search_n_words(ffs),
            sizeof(*ffs->word_idx));
    ffs->word_list = garray_init(0, sizeof(int32));
    ffs->utt_vocab = bitvec_alloc(search_n_words(ffs));
    ffs->expand_words = bitvec_alloc(search_n_words(ffs));
    /* FIXME: Make this a garray_t */
    ffs->expand_word_list = ckd_calloc(search_n_words(ffs),
            sizeof(*ffs->expand_word_list));
    ffs->bptbl = bptbl_init("fwdflat", d2p, cmd_ln_int32_r(config, "-latsize"), 256);

    /* Allocate active word list array */
    ffs->active_word_list = ckd_calloc_2d(2, search_n_words(ffs),
            sizeof(**ffs->active_word_list));
    E_INFO("Allocated %d KiB for active word list\n",
            (search_n_words(ffs) * sizeof(**ffs->active_word_list)
                    + search_n_words(ffs) * sizeof(*ffs->active_word_list)) / 1024);

    /* Load language model(s) */
    if (other) {
        ffs->lmset = ngram_model_retain(search_lmset(other));
    }
    else {
        if ((path = cmd_ln_str_r(config, "-lmctl"))) {
            ffs->lmset = ngram_model_set_read(config, path, acmod->lmath);
            if (ffs->lmset == NULL)
            {
                E_ERROR("Failed to read language model control file: %s\n",
                        path);
                goto error_out;
            }
            /* Set the default language model if needed. */
            if ((path = cmd_ln_str_r(config, "-lmname")))
            {
                ngram_model_set_select(ffs->lmset, path);
            }
        }
        else if ((path = cmd_ln_str_r(config, "-lm")))
        {
            static const char *name = "default";
            ngram_model_t *lm;

            lm = ngram_model_read(config, path, NGRAM_AUTO, acmod->lmath);
            if (lm == NULL)
            {
                E_ERROR("Failed to read language model file: %s\n", path);
                goto error_out;
            }
            ffs->lmset = ngram_model_set_init(config,
                    &lm, (char **)&name,
                    NULL, 1);
            if (ffs->lmset == NULL)
            {
                E_ERROR("Failed to initialize language model set\n");
                goto error_out;
            }
        }
        if (ffs->lmset != NULL
                && ngram_wid(ffs->lmset, S3_FINISH_WORD) == ngram_unknown_wid(ffs->lmset))
        {
            E_ERROR("Language model/set does not contain </s>, recognition will fail\n");
            goto error_out;
        }
    }
    /* Calculate extra language model weight */
    ffs->lw = (int32)(cmd_ln_float32_r(config, "-fwdflatlw")
            / cmd_ln_float32_r(config, "-lw") * 32768);
    E_INFO("Second pass language weight %f => %d\n",
            (float)ffs->lw / 32768, ffs->lw);

    /* Load a vocabulary map if needed. */
    if ((path = cmd_ln_str_r(config, "-vm")))
    {
        FILE *fh;
        int32 ispipe;

        ffs->vmap = vocab_map_init(search_dict(ffs));
        if ((fh = fopen_comp(path, "r", &ispipe)) == NULL)
        {
            E_ERROR_SYSTEM("Failed to open vocabulary map file\n");
            goto error_out;
        }
        if (vocab_map_read(ffs->vmap, fh) < 0)
        {
            E_ERROR("Failed to read vocabulary map file\n");
            goto error_out;
        }
        fclose(fh);
    }

    /* Create word mappings. */
    fwdflat_search_update_widmap(ffs);
    return (search_t *)ffs;

    error_out:
    search_free((search_t *)ffs);
    return NULL;
}

static void fwdflat_search_free_word_chan(fwdflat_search_t *ffs, int32 w)
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

static void destroy_fwdflat_chan(fwdflat_search_t *ffs)
{
    int32 i;

    for (i = 0; i < garray_next_idx(ffs->word_list); ++i) {
        int32 wid = garray_ent(ffs->word_list, int32, i);
        assert(ffs->word_chan[wid] != NULL);
        fwdflat_search_free_word_chan(ffs, wid);
        ffs->word_chan[wid] = NULL;
    }
    garray_reset(ffs->word_list);
    bitvec_clear_all(ffs->utt_vocab, search_n_words(ffs));
}

static int fwdflat_search_free(search_t *base)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *) base;

    destroy_fwdflat_chan(ffs);

    hmm_context_free(ffs->hmmctx);
    listelem_alloc_free(ffs->chan_alloc);
    listelem_alloc_free(ffs->root_chan_alloc);
    ngram_model_free(ffs->lmset);

    ckd_free(ffs->word_idx);
    ckd_free(ffs->word_chan);
    bitvec_free(ffs->word_active);
    bitvec_free(ffs->expand_words);
    ckd_free(ffs->expand_word_list);
    bptbl_free(ffs->bptbl);
    ckd_free_2d(ffs->active_word_list);
    ckd_free(ffs->rcss);

    garray_free(ffs->word_list);
    bitvec_free(ffs->utt_vocab);

    return 0;
}

static internal_node_t *
fwdflat_search_alloc_all_rc(fwdflat_search_t *ffs, int32 w)
{
    internal_node_t *fhmm, *hmm;
    xwdssid_t *rssid;
    int32 i, tmatid, ciphone;

    /* DICT2PID */
    /* Get pointer to array of triphones for final diphone. */assert(!dict_is_single_phone(search_dict(ffs), w));
    ciphone = dict_last_phone(search_dict(ffs),w);
    rssid = dict2pid_rssid(search_dict2pid(ffs),
            ciphone,
            dict_second_last_phone(search_dict(ffs),w));
    tmatid = bin_mdef_pid2tmatid(search_acmod(ffs)->mdef, ciphone);
    fhmm = hmm = listelem_malloc(ffs->chan_alloc);
    hmm->rc_id = 0;
    hmm->next = NULL;
    hmm->ciphone = dict_last_phone(search_dict(ffs),w);
    hmm_init(ffs->hmmctx, &hmm->hmm, FALSE, rssid->ssid[0], tmatid);
    E_DEBUG(3,("allocated rc_id 0 ssid %d ciphone %d lc %d word %s\n",
                    rssid->ssid[0], hmm->ciphone,
                    dict_second_last_phone(search_dict(ffs),w),
                    dict_wordstr(search_dict(ffs),w)));
    for (i = 1; i < rssid->n_ssid; ++i) {
        if ((hmm->next == NULL) || (hmm_nonmpx_ssid(&hmm->next->hmm)
                != rssid->ssid[i])) {
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
                            dict_second_last_phone(search_dict(ffs),w),
                            dict_wordstr(search_dict(ffs),w)));
        }
        else
            hmm = hmm->next;
    }

    return fhmm;
}

static void build_fwdflat_word_chan(fwdflat_search_t *ffs, int32 wid)
{
    int32 p;
    first_node_t *rhmm;
    internal_node_t *hmm, *prevhmm;
    dict_t *dict;
    dict2pid_t *d2p;

    if (ffs->word_chan[wid] != NULL)
        return;

    dict = search_dict(ffs);
    d2p = search_dict2pid(ffs);

    /* Multiplex root HMM for first phone (one root per word, flat
     * lexicon). */
    rhmm = listelem_malloc(ffs->root_chan_alloc);
    if (dict_is_single_phone(dict, wid)) {
        rhmm->ciphone = dict_first_phone(dict, wid);
        rhmm->ci2phone = bin_mdef_silphone(search_acmod(ffs)->mdef);
    }
    else {
        rhmm->ciphone = dict_first_phone(dict, wid);
        rhmm->ci2phone = dict_second_phone(dict, wid);
    }
    rhmm->next = NULL;
    hmm_init(ffs->hmmctx, &rhmm->hmm, TRUE,
            bin_mdef_pid2ssid(search_acmod(ffs)->mdef, rhmm->ciphone),
            bin_mdef_pid2tmatid(search_acmod(ffs)->mdef, rhmm->ciphone));

    /* HMMs for word-internal phones */
    prevhmm = NULL;
    for (p = 1; p < dict_pronlen(dict, wid) - 1; p++) {
        hmm = listelem_malloc(ffs->chan_alloc);
        hmm->ciphone = dict_pron(dict, wid, p);
        hmm->rc_id = -1;
        hmm->next = NULL;
        hmm_init(ffs->hmmctx, &hmm->hmm, FALSE, dict2pid_internal(d2p, wid, p),
                bin_mdef_pid2tmatid(search_acmod(ffs)->mdef,
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
    bitvec_set(ffs->utt_vocab, wid);
    garray_append(ffs->word_list, &wid);
}

static int fwdflat_search_start(search_t *base)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *) base;
    first_node_t *rhmm;
    int i;

    bptbl_reset(ffs->bptbl);
    ffs->oldest_bp = -1;
    for (i = 0; i < search_n_words(ffs); i++)
        ffs->word_idx[i] = NO_BP;

    /* Reset output arc buffer. */
    if (search_output_arcs(ffs))
        arc_buffer_producer_start_utt(search_output_arcs(ffs), base->uttid);

    /* Create word HMM for start, end, and silence words. */
    for (i = search_start_wid(ffs); i < search_n_words(ffs); ++i)
        build_fwdflat_word_chan(ffs, i);

    /* Start search with <s> */
    rhmm = (first_node_t *) ffs->word_chan[search_start_wid(ffs)];
    hmm_enter(&rhmm->hmm, 0, NO_BP, 0);
    ffs->active_word_list[0][0] = search_start_wid(ffs);
    ffs->n_active_word[0] = 1;

    bitvec_clear_all(ffs->expand_words, search_n_words(ffs));
    ffs->best_score = 0;
    ffs->renormalized = FALSE;

    ffs->bptbl->n_frame = 0;
    ffs->st.n_fwdflat_chan = 0;
    ffs->st.n_fwdflat_words = 0;
    ffs->st.n_fwdflat_word_transition = 0;
    ffs->st.n_senone_active_utt = 0;

    search_call_event(base, SEARCH_START_UTT, base->acmod->output_frame);

    return 0;
}

static int32 update_oldest_bp(fwdflat_search_t *ffs, hmm_t *hmm)
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

static void compute_fwdflat_sen_active(fwdflat_search_t *ffs, int frame_idx)
{
    int32 i, w;
    int32 *awl;
    first_node_t *rhmm;
    internal_node_t *hmm;

    acmod_clear_active(search_acmod(ffs));
    ffs->oldest_bp = bptbl_end_idx(ffs->bptbl);

    i = ffs->n_active_word[frame_idx & 0x1];
    awl = ffs->active_word_list[frame_idx & 0x1];

    for (w = *(awl++); i > 0; --i, w = *(awl++)) {
        rhmm = (first_node_t *) ffs->word_chan[w];
        if (hmm_frame(&rhmm->hmm) == frame_idx) {
            acmod_activate_hmm(search_acmod(ffs), &rhmm->hmm);
            update_oldest_bp(ffs, &rhmm->hmm);
        }

        for (hmm = rhmm->next; hmm; hmm = hmm->next) {
            if (hmm_frame(&hmm->hmm) == frame_idx) {
                acmod_activate_hmm(search_acmod(ffs), &hmm->hmm);
                update_oldest_bp(ffs, &hmm->hmm);
            }
        }
    }
    assert(ffs->oldest_bp < bptbl_end_idx(ffs->bptbl));
}

static void fwdflat_eval_chan(fwdflat_search_t *ffs, int frame_idx)
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
            if ((score BETTER_THAN bestscore) && (w != search_finish_wid(ffs)))
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

static void fwdflat_search_save_bp(fwdflat_search_t *ffs, int frame_idx,
        int32 w, int32 score, int32 path, int32 rc)
{
    int32 bp;

    /* Look for an existing exit for this word in this frame. */
    bp = ffs->word_idx[w];
    if (bp != NO_BP) {
        bp_t bpe;

        bptbl_get_bp(ffs->bptbl, bp, &bpe);
        /* Keep only the best scoring one (this is a potential source
         * of search errors...) */
        if (bpe.score WORSE_THAN score) {
            bptbl_update_bp(ffs->bptbl, bp, rc, path, score);
        }
        /* But do keep track of scores for all right contexts, since
         * we need them to determine the starting path scores for any
         * successors of this word exit. */
        bptbl_set_rcscore(ffs->bptbl, bp, rc, score);
    }
    else {
        bpidx_t bpidx = bptbl_enter(ffs->bptbl, w, path, score, rc);
        ffs->word_idx[w] = bpidx;
    }
}

static void fwdflat_prune_chan(fwdflat_search_t *ffs, int frame_idx)
{
    int32 i, cf, nf, w, pip, newscore, thresh, wordthresh;
    int32 *awl;
    first_node_t *rhmm;
    internal_node_t *hmm, *nexthmm;

    cf = frame_idx;
    nf = cf + 1;
    i = ffs->n_active_word[cf & 0x1];
    awl = ffs->active_word_list[cf & 0x1];
    bitvec_clear_all(ffs->word_active, search_n_words(ffs));

    thresh = ffs->best_score + ffs->fwdflatbeam;
    wordthresh = ffs->best_score + ffs->fwdflatwbeam;
    pip = ffs->pip;

    /* Scan all active words. */
    for (w = *(awl++); i > 0; --i, w = *(awl++)) {
        rhmm = (first_node_t *) ffs->word_chan[w];
        /* Propagate active root channels */
        if (hmm_frame(&rhmm->hmm) == cf && hmm_bestscore(&rhmm->hmm)
                BETTER_THAN thresh) {
            hmm_frame(&rhmm->hmm) = nf;
            bitvec_set(ffs->word_active, w);

            /* Transitions out of root channel */
            newscore = hmm_out_score(&rhmm->hmm);
            if (rhmm->next) {
                assert(!dict_is_single_phone(search_dict(ffs), w));
                newscore += pip;
                if (newscore BETTER_THAN thresh) {
                    hmm = rhmm->next;
                    /* Enter all right context phones */
                    if (hmm->rc_id >= 0) {
                        for (; hmm; hmm = hmm->next) {
                            if ((hmm_frame(&hmm->hmm) < cf) || (newscore
                                    BETTER_THAN hmm_in_score(&hmm->hmm))) {
                                hmm_enter(&hmm->hmm, newscore,
                                        hmm_out_history(&rhmm->hmm), nf);
                            }
                        }
                    }
                    /* Just a normal word internal phone */
                    else {
                        if ((hmm_frame(&hmm->hmm) < cf) || (newscore
                                BETTER_THAN hmm_in_score(&hmm->hmm))) {
                            hmm_enter(&hmm->hmm, newscore,
                                    hmm_out_history(&rhmm->hmm), nf);
                        }
                    }
                }
            }
            else {
                assert(dict_is_single_phone(search_dict(ffs), w));
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
                                            || (newscore
                                                    BETTER_THAN hmm_in_score(&nexthmm->hmm))) {
                                        hmm_enter(&nexthmm->hmm, newscore,
                                                hmm_out_history(&hmm->hmm), nf);
                                    }
                                }
                            }
                            /* Enter single word-internal phone. */
                            else {
                                if ((hmm_frame(&nexthmm->hmm) < cf)
                                        || (newscore
                                                BETTER_THAN hmm_in_score(&nexthmm->hmm))) {
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
                                    hmm_out_history(&hmm->hmm), hmm->rc_id);
                        }
                    }
                }
                /* Zero out inactive HMMs. */
                else if (hmm_frame(&hmm->hmm) != nf) {
                    hmm_clear(&hmm->hmm);
                }
            }
        }
    }
}

static void fwdflat_word_transition(fwdflat_search_t *ffs, int frame_idx)
{
    int32 cf, nf, b, thresh, pip, i, w, newscore;
    int32 best_silrc_score = 0, best_silrc_bp = 0; /* FIXME: good defaults? */
    first_node_t *rhmm;
    int32 *awl;
    dict_t *dict = search_dict(ffs);
    dict2pid_t *d2p = search_dict2pid(ffs);

    cf = frame_idx;
    nf = cf + 1;
    thresh = ffs->best_score + ffs->fwdflatbeam;
    pip = ffs->pip;
    best_silrc_score = WORST_SCORE;

    /* Scan words exited in current frame */
    for (b = bptbl_ef_idx(ffs->bptbl, cf); b < bptbl_ef_idx(ffs->bptbl, cf + 1); b++) {
        xwdssid_t *rssid;
        int32 silscore;
        bp_t ent;

        bptbl_get_bp(ffs->bptbl, b, &ent);
        ffs->word_idx[ent.wid] = NO_BP;

        if (ent.wid == search_finish_wid(ffs))
            continue;

        /* Get the mapping from right context phone ID to index in the
         * right context table and the bptbl->bscore_stack. */
        bptbl_get_rcscores(ffs->bptbl, b, ffs->rcss);
        if (ent.last2_phone == -1)
            rssid = NULL;
        else
            rssid = dict2pid_rssid(d2p, ent.last_phone, ent.last2_phone);

        /* Transition to all successor words. */
        for (i = 0; i < ffs->n_expand_word; ++i) {
            int32 lmscore, n_used;

            w = ffs->expand_word_list[i];
            /* Get the exit score we recorded in save_bwd_ptr(), or
             * something approximating it. */
            if (rssid)
                newscore = ffs->rcss[rssid->cimap[dict_first_phone(dict, w)]];
            else
                newscore = ffs->rcss[0];
            if (newscore == WORST_SCORE)
                continue;
            lmscore = ngram_tg_score(ffs->lmset, dict_basewid(dict, w),
                    ent.real_wid, ent.prev_real_wid, &n_used) >> SENSCR_SHIFT;
            lmscore = lmscore * ffs->lw / 32768;
            newscore += lmscore;
            newscore += pip >> SENSCR_SHIFT;

            /* Enter the next word */
            if (newscore BETTER_THAN thresh) {
                rhmm = (first_node_t *) ffs->word_chan[w];
                if ((hmm_frame(&rhmm->hmm) < cf) || (newscore
                        BETTER_THAN hmm_in_score(&rhmm->hmm))) {
                    hmm_enter(&rhmm->hmm, newscore, b, nf);
                    /* DICT2PID: This is where mpx ssids get introduced. */
                    /* Look up the ssid to use when entering this mpx triphone. */hmm_mpx_ssid(&rhmm->hmm, 0)
                            = dict2pid_ldiph_lc(d2p, rhmm->ciphone, rhmm->ci2phone,
                                    dict_last_phone(dict, ent.wid));
                    assert(IS_S3SSID(hmm_mpx_ssid(&rhmm->hmm, 0)));
                    bitvec_set(ffs->word_active, w);
                }
            }
        }

        /* Get the best exit into silence. */
        if (rssid)
            silscore = ffs->rcss[rssid->cimap[search_acmod(ffs)->mdef->sil]];
        else
            silscore = ffs->rcss[0];
        if (silscore BETTER_THAN best_silrc_score) {
            best_silrc_score = silscore;
            best_silrc_bp = b;
        }
    }

    /* Transition to <sil> */
    newscore = best_silrc_score + ffs->silpen + pip;
    if ((newscore BETTER_THAN thresh) && (newscore BETTER_THAN WORST_SCORE)) {
        w = search_silence_wid(ffs);
        rhmm = (first_node_t *) ffs->word_chan[w];
        if ((hmm_frame(&rhmm->hmm) < cf) || (newscore
                BETTER_THAN hmm_in_score(&rhmm->hmm))) {
            hmm_enter(&rhmm->hmm, newscore, best_silrc_bp, nf);
            bitvec_set(ffs->word_active, w);
        }
    }
    /* Transition to noise words (FIXME: Depends on dictionary organization) */
    newscore = best_silrc_score + ffs->fillpen + pip;
    if ((newscore BETTER_THAN thresh) && (newscore BETTER_THAN WORST_SCORE)) {
        for (w = search_silence_wid(ffs) + 1; w < search_n_words(ffs); w++) {
            rhmm = (first_node_t *) ffs->word_chan[w];
            /* Noise words that aren't a single phone will have NULL here. */
            if (rhmm == NULL)
                continue;
            if ((hmm_frame(&rhmm->hmm) < cf) || (newscore
                    BETTER_THAN hmm_in_score(&rhmm->hmm))) {
                hmm_enter(&rhmm->hmm, newscore, best_silrc_bp, nf);
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
            hmm_clear(&rhmm->hmm);
        }
    }
}

static void fwdflat_renormalize_scores(fwdflat_search_t *ffs, int frame_idx,
        int32 norm)
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

static int fwdflat_create_expand_word_list(fwdflat_search_t *ffs)
{
    int32 i, j;

    for (i = 0, j = 0; i < garray_next_idx(ffs->word_list); ++i) {
        int32 wid = garray_ent(ffs->word_list, int32, i);
        if (bitvec_is_set(ffs->expand_words, wid))
            ffs->expand_word_list[j++] = wid;
    }
    ffs->n_expand_word = j;
    ffs->st.n_fwdflat_word_transition += ffs->n_expand_word;
    return 0;
}

#ifdef __GNUC__
#define ATTRIBUTE_UNUSED __attribute__((unused))
#else
#define ATTRIBUTE_UNUSED
#endif

static void ATTRIBUTE_UNUSED
fwdflat_dump_expand_words(fwdflat_search_t *ffs, int sf)
{
    int i;

    E_INFO("Frame %d word list:", sf);
    for (i = 0; i < ffs->n_expand_word; ++i)
    {
        E_INFOCONT(" %s",
                dict_wordstr(search_dict(ffs), ffs->expand_word_list[i]));
    }
    E_INFOCONT("\n");
}

static int fwdflat_create_active_word_list(fwdflat_search_t *ffs, int nf)
{
    int32 *nawl, i, j;

    nawl = ffs->active_word_list[nf & 0x1];
    for (i = j = 0; i < garray_next_idx(ffs->word_list); ++i) {
        int32 wid = garray_ent(ffs->word_list, int32, i);
        if (bitvec_is_set(ffs->word_active, wid)) {
            *(nawl++) = wid;
            ++j;
        }
    }
    ffs->n_active_word[nf & 0x1] = j;

    return 0;
}

static int update_partials(fwdflat_search_t *ffs, int frame_idx)
{
    ffs->best_exit = bptbl_find_exit(ffs->bptbl, -1);
    if (ffs->best_exit != NO_BP) {
        int32
                bbp =
                        dict_basewid(search_dict(ffs), bptbl_ent(ffs->bptbl, ffs->best_exit)->wid);
        if (bbp != ffs->best_exit_wid) {
            ffs->best_exit_wid = bbp;
            search_call_event(search_base(ffs), SEARCH_PARTIAL_RESULT,
                    frame_idx);
        }
    }
    return ffs->best_exit;
}

static int fwdflat_search_one_frame(fwdflat_search_t *ffs, int frame_idx)
{
    acmod_t *acmod = search_acmod(ffs);
    int16 const *senscr;
    int fi;

    E_DEBUG(2,("Searching frame %d\n", frame_idx));
    /* Activate our HMMs for the current frame if need be. */
    if (!acmod->compallsen)
        compute_fwdflat_sen_active(ffs, frame_idx);

    /* Compute GMM scores for the current frame. */
    senscr = acmod_score(acmod, frame_idx);
    ffs->st.n_senone_active_utt += acmod->n_senone_active;

    /* Mark backpointer table for current frame. */
    fi = bptbl_push_frame(ffs->bptbl, ffs->oldest_bp);
    assert(fi == frame_idx);

    /* Forward retired backpointers to the arc buffer. */
    if (search_output_arcs(ffs))
        arc_buffer_producer_sweep(search_output_arcs(ffs),
                /* FIXME: make configurable */FALSE);

    /* If the best score is equal to or worse than WORST_SCORE,
     * recognition has failed, don't bother to keep trying. */
    if (ffs->best_score == WORST_SCORE || ffs->best_score
            WORSE_THAN WORST_SCORE)
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
    fwdflat_create_expand_word_list(ffs);
    /* fwdflat_dump_expand_words(ffs, frame_idx); */
    fwdflat_word_transition(ffs, frame_idx);
    fwdflat_create_active_word_list(ffs, frame_idx + 1);

    /* Update partial results. */
    update_partials(ffs, frame_idx);

    /* Release the frame just searched. */
    acmod_consumer_release(acmod, frame_idx);

    /* Return the number of frames processed. */
    return 1;
}

static void fwdflat_search_add_expand_word(fwdflat_search_t *ffs, int32 wid)
{
    dict_t *dict = search_dict(ffs);

    if (!bitvec_is_set(ffs->expand_words, wid)) {
        /* Test this after the bitvec to avoid looking up the same
         * N-Gram repeatedly (unfortunately for unknown words that
         * will happen anyway, which may be a problem). */
        if (!ngram_model_set_known_wid(ffs->lmset, dict_basewid(dict, wid)))
            return;
        bitvec_set(ffs->expand_words, wid);
        build_fwdflat_word_chan(ffs, wid);
    }
}

static int fwdflat_search_expand_arcs(fwdflat_search_t *ffs, int sf, int ef)
{
    arc_t *arc_start, *arc_end, *arc;

    arc_start = arc_buffer_iter(search_input_arcs(ffs), sf);
    arc_end = arc_buffer_iter(search_input_arcs(ffs), ef);
    E_DEBUG(2,("Expanding %ld arcs in %d:%d\n",
                    arc_end > arc_start ? arc_end - arc_start : 0, sf, ef));
    bitvec_clear_all(ffs->expand_words, search_n_words(ffs));
    for (arc = arc_start; arc != arc_end; arc = arc_buffer_iter_next(
            search_input_arcs(ffs), arc)) {
        /* Expand things in the vocabulary map, if we have one. */
        if (ffs->vmap) {
            int32 const *wids;
            int32 nwids;
            if ((wids = vocab_map_unmap(ffs->vmap, arc->wid, &nwids)) != NULL) {
                int32 i;
                for (i = 0; i < nwids; ++i)
                    fwdflat_search_add_expand_word(ffs, wids[i]);
                continue;
            }
        }
        /* Otherwise, just expand the word normally. */
        fwdflat_search_add_expand_word(ffs, arc->wid);
    }
    return 0;
}

static int fwdflat_search_decode(search_t *base)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *) base;
    acmod_t *acmod = search_acmod(base);
    int frame_idx;

    ptmr_start(&ffs->base.t);
    frame_idx = 0;
    /* Wait for the arc buffer. */E_INFO("fwdflat: waiting for acmod start\n");
    if (acmod_consumer_start_utt(base->acmod, -1) < 0)
    {
        if (base->output_arcs)
        arc_buffer_producer_shutdown(base->output_arcs);
        return -1;
    }
    base->uttid = base->acmod->uttid;
    E_INFO("waiting for arc buffer start\n");
    if (arc_buffer_consumer_start_utt(search_input_arcs(base), -1) < 0)
        return -1;
    fwdflat_search_start(base);
    while (!acmod_eou(base->acmod))
    {
        int end_win;

        /* Stop timing and wait for the arc buffer. */
        ptmr_stop(&ffs->base.t);
        if (arc_buffer_consumer_wait(search_input_arcs(ffs), -1) < 0)
        goto canceled;

        /* Figure out the last frame we need. */
        end_win = frame_idx + ffs->max_sf_win;
        /* Decode as many frames as possible. */
        while (arc_buffer_eou(search_input_arcs(ffs))
                || arc_buffer_iter(search_input_arcs(ffs),
                        end_win - 1) != NULL)
        {
            int start_win, k;

            /* Note that if search_input_arcs(ffs)->final changes state
             * between the tests above and now, there will be no ill
             * effect, since we will never block on acmod if we
             * believe it to be false.  A full analysis by cases
             * follows:
             *
             * 1) final = TRUE on exiting arc_buffer_wait(): it
             *    will not be reset until the next utterance, we wait
             *    on acmod below exclusively until utterance end.
             *
             * 2) final = FALSE on exiting arc_buffer_wait():
             *
             *    2a) final = FALSE in while() above: arc buffer
             *        either does or does not contain new frames, we
             *        do not block either way.  If arc buffer becomes
             *        final in the meantime, arc_buffer_wait will not
             *        block above (FIXME: This relies on there only
             *        being one consumer thread)
             *
             *    2b) final = TRUE in while() above: it will not be
             *        reset until the next utterance and we wait on
             *        acmod below exclusively
             */
            if (arc_buffer_eou(search_input_arcs(ffs)))
            {
                int nfx;
                /* We are no longer waiting on the arc buffer so it is
                 * okay to wait as long as necessary for acmod
                 * frames. */
                if ((nfx = acmod_consumer_wait(acmod, -1)) < 0)
                {
                    if (acmod_eou(acmod))
                    break;
                    else
                    goto canceled;
                }
            }
            else
            {
                int nfx;
                /* Don't wait on the acoustic model as that could
                 * cause deadlock.  If this times out it will return
                 * -1, so there is no danger of accidentally searching
                 * the same frame twice. */
                if ((nfx = acmod_consumer_wait(acmod, 0)) < 0)
                {
                    if (acmod_eou(acmod))
                    break;
                    else
                    goto canceled;
                }
            }
            ptmr_start(&ffs->base.t);

            /* Lock the arc buffer while we expand arcs. */
            arc_buffer_lock(search_input_arcs(ffs));
            end_win = frame_idx + ffs->max_sf_win;
            start_win = frame_idx - ffs->max_sf_win;
            if (start_win < 0) start_win = 0;
            fwdflat_search_expand_arcs(ffs, start_win, end_win);
            arc_buffer_unlock(search_input_arcs(ffs));

            /* Now do our search. */
            if ((k = fwdflat_search_one_frame(ffs, frame_idx)) <= 0)
            break;
            frame_idx += k;
            arc_buffer_consumer_release(search_input_arcs(ffs), start_win);
            ptmr_stop(&ffs->base.t);
        }
    }
    arc_buffer_consumer_end_utt(search_input_arcs(ffs));
    ptmr_start(&ffs->base.t);
    fwdflat_search_finish(search_base(ffs));
    ptmr_stop(&ffs->base.t);
    return frame_idx;

    canceled:
    if (base->output_arcs)
    arc_buffer_producer_shutdown(base->output_arcs);
    return -1;
}

static int fwdflat_search_finish(search_t *base)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *) base;
    acmod_t *acmod = search_acmod(ffs);
    int cf;

    /* This is the number of frames of input. */
    cf = acmod->output_frame;

    /* Final result for this thread. */
    search_call_event(base, SEARCH_FINAL_RESULT, cf);

    /* Finalize the backpointer table. */
    bptbl_finalize(ffs->bptbl);

    /* Finalize the output arc buffer and wait for consumer. */
    if (search_output_arcs(ffs))
        arc_buffer_producer_end_utt(search_output_arcs(ffs), FALSE);

    /* Finalize the input acmod (signals producer) */
    acmod_consumer_end_utt(base->acmod);
    base->total_frames += base->acmod->output_frame;

    search_call_event(base, SEARCH_END_UTT, base->acmod->output_frame);

    /* Print out some statistics. */
    if (cf > 0) {
        E_INFO("%8d words recognized in %d frames (%d/fr)\n",
        bptbl_end_idx(ffs->bptbl), cf + 1,
        (bptbl_end_idx(ffs->bptbl) + (cf >> 1)) / (cf + 1));
        E_INFO("%8d senones evaluated (%d/fr)\n", ffs->st.n_senone_active_utt,
                (ffs->st.n_senone_active_utt + (cf >> 1)) / (cf + 1));
        E_INFO("%8d channels searched (%d/fr)\n",
                ffs->st.n_fwdflat_chan, ffs->st.n_fwdflat_chan / (cf + 1));
        E_INFO("%8d words searched (%d/fr)\n",
                ffs->st.n_fwdflat_words, ffs->st.n_fwdflat_words / (cf + 1));
        E_INFO("%8d word transitions (%d/fr)\n",
                ffs->st.n_fwdflat_word_transition,
                ffs->st.n_fwdflat_word_transition / (cf + 1));
        E_INFO("time %f wall %.2f xRT\n",
                base->t.t_elapsed,
                base->t.t_elapsed / base->acmod->output_frame
                * cmd_ln_int32_r(base->config, "-frate"));
        E_INFO("utterance vocabulary had %d words\n",
                garray_next_idx(ffs->word_list));
    }

    /* bptbl_dump(ffs->bptbl); */

    /* Reset the utterance vocabulary. */
    destroy_fwdflat_chan(ffs);

    return 0;
}

static int32 fwdflat_search_prob(search_t *base)
{
    /* FIXME: Going to estimate this from partial results in the future. */
    return 0;
}

static char const *
fwdflat_search_hyp(search_t *base, int32 *out_score)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *)base;

    ckd_free(base->hyp_str);
    if (bptbl_is_final(ffs->bptbl)) {
        base->hyp_str = bptbl_hyp(ffs->bptbl, out_score, base->finish_wid);
    }
    else {
        *out_score = ffs->best_score;
        base->hyp_str = bptbl_backtrace(ffs->bptbl, ffs->best_exit);
    }
    return base->hyp_str;
}

static seg_iter_t *
fwdflat_search_seg_iter(search_t *base, int32 *out_score)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *) base;
    if (bptbl_is_final(ffs->bptbl))
        return bptbl_seg_iter(ffs->bptbl, out_score, search_finish_wid(ffs));
    else {
        *out_score = ffs->best_score;
        return bptbl_seg_backtrace(ffs->bptbl, ffs->best_exit);
    }
}

vocab_map_t *
fwdflat_search_set_vocab_map(search_t *search, vocab_map_t *vm)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *) search;
    ffs->vmap = vm;
    return vm;
}

static bptbl_t *
fwdflat_search_bptbl(search_t *base)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *) base;
    return ffs->bptbl;
}

static ngram_model_t *
fwdflat_search_lmset(search_t *base)
{
    fwdflat_search_t *ffs = (fwdflat_search_t *) base;
    return ffs->lmset;
}
