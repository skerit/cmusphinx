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
 * @file fwdtree_search.c Lexicon tree search.
 */

/* System headers. */
#include <string.h>
#include <assert.h>

/* SphinxBase headers. */
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/listelem_alloc.h>
#include <sphinxbase/err.h>

/* Local headers. */
#include "tmat.h"
#include "hmm.h"
#include "search_internal.h"
#include "fwdtree_search.h"
#include "arc_buffer.h"

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
 * Lexical tree node data type.
 *
 * Not the first HMM for words, which multiplex HMMs based on
 * different left contexts.  This structure is used both in the
 * dynamic HMM tree structure and in the per-word last-phone right
 * context fanout.
 */
typedef struct nonroot_node_s {
    hmm_t hmm;                  /**< Basic HMM structure.  This *must*
                                   be first in the structure because
                                   nonroot_node_t and root_node_t are
                                   sometimes used interchangeably */
    struct nonroot_node_s *next;/**< first descendant of this channel; or, in the
				   case of the last phone of a word, the next
				   alternative right context channel */
    struct nonroot_node_s *alt;	/**< sibling; i.e., next descendant of parent HMM */

    int32    ciphone;		/**< ciphone for this node */
    union {
	int32 penult_phn_wid;	/**< list of words whose last phone follows this one;
				   this field indicates the first of the list; the
				   rest must be built up in a separate array.  Used
				   only within HMM tree.  -1 if none */
	int32 rc_id;		/**< right-context id for last phone of words */
    } info;
} nonroot_node_t;

/**
 * Lexical tree node data type for the first phone (root) of each dynamic HMM tree
 * structure.
 *
 * Each state may have a different parent static HMM.  Most fields are
 * similar to those in nonroot_node_t.
 */
typedef struct root_node_s {
    hmm_t hmm;                  /**< Basic HMM structure.  This *must* be first in
                                   the structure because nonroot_node_t and root_node_t are
                                   sometimes used interchangeably. */
    nonroot_node_t *next;	/**< first descendant of this channel */

    int32  penult_phn_wid;
    int32  this_phn_wid;	/**< list of words consisting of this single phone;
				   actually the first of the list, like penult_phn_wid;
				   -1 if none */
    int16    ciphone;		/**< first ciphone of this node; all words rooted at this
				   node begin with this ciphone */
    int16    ci2phone;		/**< second ciphone of this node; one root HMM for each
                                   unique right context */
} root_node_t;

/**
 * Candidate words for entering their last phones.
 *
 * Cleared and rebuilt in each frame.  NOTE: candidates can only be
 * multi-phone, real dictionary words.
 */
typedef struct lastphn_cand_s {
    int32 wid;           /**< Word ID. */
    int32 score;         /**< Path score. */
    int32 bp;            /**< Backpointer ID of path. */
    int32 next;          /**< Next candidate starting at the same frame. */
} lastphn_cand_t;

/**
 * Cached precomputed component of path score for words entering last
 * frame.
 *
 * Since the same instance of a word (i.e., <word,start-frame>)
 * reaches its last phone several times, we can compute its best BP
 * and LM transition score info just the first time and cache it for
 * future occurrences.
 */
typedef struct {
    int32 sf;                   /**< Start frame */
    int32 dscr;                 /**< Delta-score upon entering last
                                 * phone, i.e. path score for just
                                 * this word's frames */
    int32 bp;                   /**< Best BP */
} last_ltrans_t;

/** Size of increments for reallocating cand_sf */
#define CAND_SF_ALLOCSIZE	32
/**
 * Frame, index pair for external sort of last_phn_cand.
 */
typedef struct {
    int32 bp_ef;        /**< Exit frame for previous word. */
    int32 cand;         /**< Index into last_phn_cand. */
} cand_sf_t;

/**
 * Structure for reorganizing the BP table entries in the current
 * frame according to distinct right context ci-phones.
 *
 * Each entry contains the best BP entry for a given right context.
 * Each successor word will pick up the correct entry based on its
 * first ci-phone.
 */
typedef struct bestbp_rc_s {
    int32 score;
    int32 path;     /**< BP table index corresponding to this entry */
    int32 lc;       /**< right most ci-phone of above BP entry word */
} bestbp_rc_t;

/**
 * Various statistics for profiling.
 */
typedef struct fwdtree_stats_s {
    int32 n_phone_eval;
    int32 n_root_chan_eval;
    int32 n_nonroot_chan_eval;
    int32 n_last_chan_eval;
    int32 n_word_lastchan_eval;
    int32 n_lastphn_cand_utt;
    int32 n_senone_active_utt;
} fwdtree_stats_t;

/**
 * Approximate tree-based forward search.
 */
typedef struct fwdtree_search_s {
    search_t base;
    ngram_model_t *lmset;
    hmm_context_t *hmmctx;

    listelem_alloc_t *chan_alloc;      /**< For nonroot_node_t */
    listelem_alloc_t *root_chan_alloc; /**< For root_node_t */

    /**
     * Backpointer table (temporary storage for active word arcs).
     */
    bptbl_t *bptbl;
    int32 oldest_bp; /**< Oldest bptable entry active in decoding graph. */

    int32 *word_idx;     /**< BPTable index for any word in current frame;
                            cleared before each frame */
    int32 *rcss; /**< Temporary storage for right context scores. */

    /**
     * Search structure of HMM instances.
     *
     * The word triphone sequences (HMM instances) are transformed
     * into tree structures, one tree per unique left triphone in the
     * entire dictionary (actually diphone, since its left context
     * varies dyamically during the search process).  The entire set
     * of trees of channels is allocated once and for all during
     * initialization (since dynamic management of active CHANs is
     * time consuming), with one exception: the last phones of words,
     * that need multiple right context modelling, are not maintained
     * in this static structure since there are too many of them and
     * few are active at any time.  Instead they are maintained as
     * linked lists of CHANs, one list per word, and each CHAN in this
     * set is allocated only on demand and freed if inactive.
     */
    root_node_t *root_chan;  /**< Roots of search tree. */
    int32 n_root_chan_alloc; /**< Number of root_chan allocated */
    int32 n_root_chan;       /**< Number of valid root_chan */
    int32 n_nonroot_chan;    /**< Number of valid non-root channels */
    int32 max_nonroot_chan;  /**< Maximum possible number of non-root channels */
    root_node_t *rhmm_1ph;   /**< Root HMMs for single-phone words */

    /**
     * Channels associated with a given word (only used for right
     * contexts and single-phone words).  WARNING: For single-phone
     * words this actually contains pointers to root_node_t, which are
     * allocated using root_chan_alloc.  This is a suboptimal state of
     * affairs.
     */
    nonroot_node_t **word_chan;
    bitvec_t *word_active;      /**< array of active flags for all words. */

    /**
     * Each node in the HMM tree structure may point to a set of words
     * whose last phone would follow that node in the tree structure
     * (but is not included in the tree structure for reasons
     * explained above).  The channel node points to one word in this
     * set of words.  The remaining words are linked through
     * homophone_set[].
     * 
     * Single-phone words are not represented in the HMM tree; they
     * are kept in word_chan.
     *
     * Specifically, homophone_set[w] = wid of next word in the same
     * set as w.
     */
    int32 *homophone_set;
    int32 *single_phone_wid; /**< list of single-phone word ids */
    int32 n_1ph_words;       /**< Number single phone words in dict (total) */
    int32 n_1ph_LMwords;     /**< Number single phone dict words also in LM;
                                these come first in single_phone_wid */
    /**
     * Array of active channels for current and next frame.
     *
     * In any frame, only some HMM tree nodes are active.
     * active_chan_list[f mod 2] = list of nonroot channels in the HMM
     * tree active in frame f.
     */
    nonroot_node_t ***active_chan_list;
    int32 n_active_chan[2];  /**< Number entries in active_chan_list */
    /**
     * Array of active multi-phone words for current and next frame.
     *
     * Similarly to active_chan_list, active_word_list[f mod 2] = list
     * of word ids for which active channels exist in word_chan in
     * frame f.
     *
     * Statically allocated single-phone words are always active and
     * should not appear in this list.
     */
    int32 **active_word_list;
    int32 n_active_word[2]; /**< Number entries in active_word_list */

    /**
     * Array of of words entering their last phone in the current frame.
     */
    lastphn_cand_t *lastphn_cand;
    int32 n_lastphn_cand; /**< Number of entries in lastphn_cand. */

    /**
     * Cached partial path scores for words entering last phone.
     */
    last_ltrans_t *last_ltrans;

    /**
     * Array of pointers into lastphn_cand, sorted by start frame.
     *
     * This is only used in last_phone_transition() and is cleared at
     * the beginning of each call.  The list of candidates is chained
     * via the cand and next fields in this and lastphn_cand,
     * respectively.
     */
    cand_sf_t *cand_sf;
    int32 cand_sf_alloc; /**< Number of cand_sf_t allocated. */

    /**
     * Best word exiting in this frame with each right context.
     */
    bestbp_rc_t *bestbp_rc;

    /**
     * Best word exit in the most recent frame (used for partial results).
     */
    bpidx_t best_exit;
    /**
     * Final word in best path (used for partial results).
     */
    int32 best_exit_wid;

    uint16 *zeroPermTab; /**< Null right context table (contains zeros) */

    int32 best_score; /**< Best Viterbi path score. */
    int32 last_phone_best_score; /**< Best Viterbi path score for last phone. */
    int32 renormalized; /**< renormalized? (FIXME: allow multiple renorms) */

    fwdtree_stats_t st; /**< Various statistics for profiling. */

    /* A children's treasury of beam widths. */
    int32 beam;
    int32 dynamic_beam;
    int32 pbeam;
    int32 wbeam;
    int32 lpbeam;
    int32 lponlybeam;
    int32 fillpen;
    int32 silpen;
    int32 wip;
    int32 nwpen;
    int32 pip;
    int32 maxwpf;
    int32 maxhmmpf;

    /** Maximum number of frames a silence word is allowed to persist. */
    int max_silence;
} fwdtree_search_t;

static int fwdtree_search_decode(search_t *base);
static int fwdtree_search_free(search_t *base);
static char const *fwdtree_search_hyp(search_t *base, int32 *out_score);
static int32 fwdtree_search_prob(search_t *base);
static seg_iter_t *fwdtree_search_seg_iter(search_t *base, int32 *out_score);
static bptbl_t *fwdtree_search_bptbl(search_t *base);
static ngram_model_t *fwdtree_search_lmset(search_t *base);
static search_t *fwdtree_search_init(search_t *other, cmd_ln_t *config,
        acmod_t *acmod, dict2pid_t *d2p);


static searchfuncs_t fwdtree_funcs = {
    /* name: */   "fwdtree",
    /* init: */   fwdtree_search_init,
    /* free: */   fwdtree_search_free,
    /* decode: */ fwdtree_search_decode,
    /* hyp: */      fwdtree_search_hyp,
    /* prob: */     fwdtree_search_prob,
    /* seg_iter: */ fwdtree_search_seg_iter,
    /* bptbl: */  fwdtree_search_bptbl,
    /* lmset: */ fwdtree_search_lmset
};

static void fwdtree_search_free_all_rc(fwdtree_search_t *fts, int32 w);
static void fwdtree_search_save_bp(fwdtree_search_t *fts, int frame_idx,
                                   int32 w, int32 score, int32 path, int32 rc);
static void fwdtree_search_alloc_all_rc(fwdtree_search_t *fts, int32 w);
static int32 fwdtree_search_exit_score(fwdtree_search_t *fts, bpidx_t bp,
                                       int last_phone, int last2_phone, int rcphone);
static void fwdtree_search_update_widmap(fwdtree_search_t *fts);
static void fwdtree_search_calc_beams(fwdtree_search_t *fts);
static void init_search_tree(fwdtree_search_t *fts);
static void init_nonroot_chan(fwdtree_search_t *fts, nonroot_node_t * hmm,
                              int32 ph, int32 ci, int32 tmatid);
static void create_search_tree(fwdtree_search_t *fts);
static void reinit_search_subtree(fwdtree_search_t *fts, nonroot_node_t * hmm);
static void reinit_search_tree(fwdtree_search_t *fts);
static void deinit_search_tree(fwdtree_search_t *fts);
static void deactivate_channels(fwdtree_search_t *fts, int frame_idx);
static void word_transition(fwdtree_search_t *fts, int frame_idx);
static void prune_channels(fwdtree_search_t *fts, int frame_idx);
static void prune_word_chan(fwdtree_search_t *fts, int frame_idx);
static void last_phone_transition(fwdtree_search_t *fts, int frame_idx);
static void prune_nonroot_chan(fwdtree_search_t *fts, int frame_idx);
static void prune_root_chan(fwdtree_search_t *fts, int frame_idx);
static int32 evaluate_channels(fwdtree_search_t *fts,
                               int16 const *senone_scores, int frame_idx);
static int32 eval_word_chan(fwdtree_search_t *fts, int frame_idx);
static int32 eval_nonroot_chan(fwdtree_search_t *fts, int frame_idx);
static int32 eval_root_chan(fwdtree_search_t *fts, int frame_idx);
static void renormalize_scores(fwdtree_search_t *fts,
                               int frame_idx, int32 norm);
static void compute_sen_active(fwdtree_search_t *fts, int frame_idx);

searchfuncs_t const *
fwdtree_search_query(void)
{
    return &fwdtree_funcs;
}

static search_t *
fwdtree_search_init(search_t *other, cmd_ln_t *config, acmod_t *acmod, dict2pid_t *d2p)
{
    fwdtree_search_t *fts;
    const char *path;

    fts = ckd_calloc(1, sizeof(*fts));
    search_base_init(&fts->base, &fwdtree_funcs, config, acmod, d2p);
    fts->hmmctx = hmm_context_init(bin_mdef_n_emit_state(acmod->mdef),
                                   acmod->tmat->tp,
                                   NULL, acmod->mdef->sseq);
    if (fts->hmmctx == NULL) {
        search_free(search_base(fts));
        return NULL;
    }
    fts->chan_alloc = listelem_alloc_init(sizeof(nonroot_node_t));
    fts->root_chan_alloc = listelem_alloc_init(sizeof(root_node_t));
    fts->rcss = ckd_calloc(bin_mdef_n_ciphone(acmod->mdef),
                           sizeof(*fts->rcss));

    /* Calculate various beam widths and such. */
    fwdtree_search_calc_beams(fts);

    /* Allocate a billion different tables for stuff. */
    fts->word_chan = ckd_calloc(dict_size(d2p->dict),
                                sizeof(*fts->word_chan));
    fts->word_active = bitvec_alloc(dict_size(d2p->dict));
    fts->last_ltrans = ckd_calloc(dict_size(d2p->dict),
                                  sizeof(*fts->last_ltrans));

    fts->bptbl = bptbl_init("fwdtree",
                            d2p, cmd_ln_int32_r(config, "-latsize"), 256);
    fts->word_idx = ckd_calloc(dict_size(d2p->dict), sizeof(*fts->word_idx));

    /* Allocate active word list array */
    fts->active_word_list = ckd_calloc_2d(2, dict_size(d2p->dict),
                                          sizeof(**fts->active_word_list));

    /* Load language model(s) */
    if ((path = cmd_ln_str_r(config, "-lmctl"))) {
        fts->lmset = ngram_model_set_read(config, path, acmod->lmath);
        if (fts->lmset == NULL) {
            E_ERROR("Failed to read language model control file: %s\n",
                    path);
            goto error_out;
        }
        /* Set the default language model if needed. */
        if ((path = cmd_ln_str_r(config, "-lmname"))) {
            ngram_model_set_select(fts->lmset, path);
        }
    }
    /* First try -fwdtreelm if it exists. */
    else if ((path = cmd_ln_str_r(config, "-fwdtreelm"))
             || (path = cmd_ln_str_r(config, "-lm"))) {
        static const char *name = "default";
        ngram_model_t *lm;

        lm = ngram_model_read(config, path, NGRAM_AUTO, acmod->lmath);
        if (lm == NULL) {
            E_ERROR("Failed to read language model file: %s\n", path);
            goto error_out;
        }
        fts->lmset = ngram_model_set_init(config,
                                          &lm, (char **)&name,
                                          NULL, 1);
        if (fts->lmset == NULL) {
            E_ERROR("Failed to initialize language model set\n");
            goto error_out;
        }
    }
    if (fts->lmset != NULL
        && ngram_wid(fts->lmset, S3_FINISH_WORD) == ngram_unknown_wid(fts->lmset)) {
        E_ERROR("Language model/set does not contain </s>, recognition will fail\n");
        goto error_out;
    }

    /* Create word mappifts. */
    fwdtree_search_update_widmap(fts);

    /* Allocate bestbp_rc, lastphn_cand, last_ltrans */
    fts->bestbp_rc = ckd_calloc(bin_mdef_n_ciphone(search_acmod(fts)->mdef),
                                sizeof(*fts->bestbp_rc));
    fts->lastphn_cand = ckd_calloc(search_n_words(fts),
                                   sizeof(*fts->lastphn_cand));
    init_search_tree(fts);
    create_search_tree(fts);

    return (search_t *)fts;

error_out:
    search_free((search_t *)fts);
    return NULL;
}

static void
fwdtree_search_update_widmap(fwdtree_search_t *fts)
{
    const char **words;
    int32 i, n_words;

    /* It's okay to include fillers since they won't be in the LM */
    n_words = search_n_words(fts);
    words = ckd_calloc(n_words, sizeof(*words));
    /* This will include alternates, again, that's okay since they aren't in the LM */
    for (i = 0; i < n_words; ++i)
        words[i] = (const char *)dict_wordstr(search_dict(fts), i);
    ngram_model_set_map_words(fts->lmset, words, n_words);
    ckd_free(words);
}

static void
fwdtree_search_calc_beams(fwdtree_search_t *fts)
{
    cmd_ln_t *config;
    acmod_t *acmod;

    config = search_config((search_t *)fts);
    acmod = search_acmod(fts);

    /* Log beam widths. */
    fts->beam = logmath_log(acmod->lmath,
                            cmd_ln_float64_r(config, "-beam")) >> SENSCR_SHIFT;
    fts->wbeam = logmath_log(acmod->lmath,
                             cmd_ln_float64_r(config, "-wbeam")) >> SENSCR_SHIFT;
    fts->pbeam = logmath_log(acmod->lmath,
                             cmd_ln_float64_r(config, "-pbeam")) >> SENSCR_SHIFT;
    fts->lpbeam = logmath_log(acmod->lmath,
                              cmd_ln_float64_r(config, "-lpbeam")) >> SENSCR_SHIFT;
    fts->lponlybeam = logmath_log(acmod->lmath, 
                                  cmd_ln_float64_r(config, "-lponlybeam")) >> SENSCR_SHIFT;

    /* Absolute pruning parameters. */
    fts->maxwpf = cmd_ln_int32_r(config, "-maxwpf");
    fts->maxhmmpf = cmd_ln_int32_r(config, "-maxhmmpf");
    fts->max_silence = cmd_ln_int32_r(config, "-maxsilfr");

    /* Various penalties which may or may not be useful. */
    fts->wip = logmath_log(acmod->lmath,
                           cmd_ln_float32_r(config, "-wip")) >> SENSCR_SHIFT;
    fts->nwpen = logmath_log(acmod->lmath,
                             cmd_ln_float32_r(config, "-nwpen")) >> SENSCR_SHIFT;
    fts->pip = logmath_log(acmod->lmath,
                           cmd_ln_float32_r(config, "-pip")) >> SENSCR_SHIFT;
    fts->silpen = logmath_log(acmod->lmath,
                              cmd_ln_float32_r(config, "-silprob")) >> SENSCR_SHIFT;
    fts->fillpen = logmath_log(acmod->lmath,
                               cmd_ln_float32_r(config, "-fillprob")) >> SENSCR_SHIFT;
}

/*
 * Allocate that part of the search channel tree structure that is independent of the
 * LM in use.
 */
static void
init_search_tree(fwdtree_search_t *fts)
{
    int32 w, ndiph, i, n_words, n_ci;
    dict_t *dict = search_dict(fts);
    bitvec_t *dimap;

    n_words = search_n_words(fts);
    fts->homophone_set = ckd_calloc(n_words, sizeof(*fts->homophone_set));

    /* Find #single phone words, and #unique first diphones (#root channels) in dict. */
    ndiph = 0;
    fts->n_1ph_words = 0;
    n_ci = bin_mdef_n_ciphone(search_acmod(fts)->mdef);
    /* Allocate a bitvector with flags for each possible diphone. */
    dimap = bitvec_alloc(n_ci * n_ci);
    for (w = 0; w < n_words; w++) {
        if (!dict_real_word(dict, w))
            continue;
        if (dict_is_single_phone(dict, w))
            ++fts->n_1ph_words;
        else {
            int ph0, ph1;
            ph0 = dict_first_phone(dict, w);
            ph1 = dict_second_phone(dict, w);
            /* Increment ndiph the first time we see a diphone. */
            if (bitvec_is_clear(dimap, ph0 * n_ci + ph1)) {
                bitvec_set(dimap, ph0 * n_ci + ph1);
                ++ndiph;
            }
        }
    }
    E_INFO("%d unique initial diphones\n", ndiph);
    bitvec_free(dimap);

    /* Add remaining dict words (</s>, <s>, <sil>, noise words) to single-phone words */
    fts->n_1ph_words += dict_num_fillers(dict) + 2;
    fts->n_root_chan_alloc = ndiph + 1;
    /* Verify that these are all *actually* single-phone words,
     * otherwise really bad thifts will happen to us. */
    for (w = 0; w < n_words; ++w) {
        if (dict_real_word(dict, w))
            continue;
        if (!dict_is_single_phone(dict, w)) {
            E_WARN("Filler word %d = %s has more than one phone, ignoring it.\n",
                   w, dict_wordstr(dict, w));
            --fts->n_1ph_words;
        }
    }

    /* Allocate and initialize root channels */
    fts->root_chan =
        ckd_calloc(fts->n_root_chan_alloc, sizeof(*fts->root_chan));
    for (i = 0; i < fts->n_root_chan_alloc; i++) {
        hmm_init(fts->hmmctx, &fts->root_chan[i].hmm, TRUE, -1, -1);
        fts->root_chan[i].penult_phn_wid = -1;
        fts->root_chan[i].next = NULL;
    }

    /* Permanently allocate and initialize channels for single-phone
     * words (1/word). */
    fts->rhmm_1ph = ckd_calloc(fts->n_1ph_words, sizeof(*fts->rhmm_1ph));
    i = 0;
    for (w = 0; w < n_words; w++) {
        if (!dict_is_single_phone(dict, w))
            continue;
        /* Use SIL as right context for these. */
        fts->rhmm_1ph[i].ci2phone = bin_mdef_silphone(search_acmod(fts)->mdef);
        fts->rhmm_1ph[i].ciphone = dict_first_phone(dict, w);
        hmm_init(fts->hmmctx, &fts->rhmm_1ph[i].hmm, TRUE,
                 bin_mdef_pid2ssid(search_acmod(fts)->mdef, fts->rhmm_1ph[i].ciphone),
                 bin_mdef_pid2tmatid(search_acmod(fts)->mdef, fts->rhmm_1ph[i].ciphone));
        fts->rhmm_1ph[i].next = NULL;

        fts->word_chan[w] = (nonroot_node_t *) &(fts->rhmm_1ph[i]);
        i++;
    }

    fts->single_phone_wid = ckd_calloc(fts->n_1ph_words,
                                       sizeof(*fts->single_phone_wid));
    E_INFO("%d root, %d non-root channels, %d single-phone words\n",
           fts->n_root_chan, fts->n_nonroot_chan, fts->n_1ph_words);
}

/*
 * One-time initialization of internal channels in HMM tree.
 */
static void
init_nonroot_chan(fwdtree_search_t *fts, nonroot_node_t * hmm, int32 ph, 
                  int32 ci, int32 tmatid)
{
    hmm->next = NULL;
    hmm->alt = NULL;
    hmm->info.penult_phn_wid = -1;
    hmm->ciphone = ci;
    hmm_init(fts->hmmctx, &hmm->hmm, FALSE, ph, tmatid);
}

/*
 * Allocate and initialize search channel-tree structure.
 * At this point, all the root-channels have been allocated and partly initialized
 * (as per init_search_tree()), and channels for all the single-phone words have been
 * allocated and initialized.  None of the interior channels of search-trees have
 * been allocated.
 * This routine may be called on every utterance, after reinit_search_tree() clears
 * the search tree created for the previous utterance.  Meant for reconfiguring the
 * search tree to suit the currently active LM.
 */
static void
create_search_tree(fwdtree_search_t *fts)
{
    nonroot_node_t *hmm;
    root_node_t *rhmm;
    int32 w, i, j, p, ph, tmatid;
    int32 n_words;
    dict_t *dict = search_dict(fts);
    dict2pid_t *d2p = search_dict2pid(fts);

    n_words = search_n_words(fts);

    E_INFO("Creating search tree\n");

    for (w = 0; w < n_words; w++)
        fts->homophone_set[w] = -1;

    E_INFO("before: %d root, %d non-root channels, %d single-phone words\n",
           fts->n_root_chan, fts->n_nonroot_chan, fts->n_1ph_words);

    fts->n_1ph_LMwords = 0;
    fts->n_root_chan = 0;
    fts->n_nonroot_chan = 0;

    for (w = 0; w < n_words; w++) {
        int ciphone, ci2phone;

        /* Ignore dictionary words not in LM */
        if (!ngram_model_set_known_wid(fts->lmset, dict_basewid(dict, w)))
            continue;

        /* Handle single-phone words individually; not in channel tree */
        if (dict_is_single_phone(dict, w)) {
            E_DEBUG(1,("single_phone_wid[%d] = %s\n",
                       fts->n_1ph_LMwords, dict_wordstr(dict, w)));
            fts->single_phone_wid[fts->n_1ph_LMwords++] = w;
            continue;
        }

        /* Find a root channel matching the initial diphone, or
         * allocate one if not found. */
        ciphone = dict_first_phone(dict, w);
        ci2phone = dict_second_phone(dict, w);
        for (i = 0; i < fts->n_root_chan; ++i) {
            if (fts->root_chan[i].ciphone == ciphone
                && fts->root_chan[i].ci2phone == ci2phone)
                break;
        }
        if (i == fts->n_root_chan) {
            rhmm = &(fts->root_chan[fts->n_root_chan]);
            rhmm->hmm.tmatid = bin_mdef_pid2tmatid(search_acmod(fts)->mdef, ciphone);
            /* Begin with CI phone?  Not sure this makes a difference... */
            hmm_mpx_ssid(&rhmm->hmm, 0) =
                bin_mdef_pid2ssid(search_acmod(fts)->mdef, ciphone);
            rhmm->ciphone = ciphone;
            rhmm->ci2phone = ci2phone;
            fts->n_root_chan++;
        }
        else
            rhmm = &(fts->root_chan[i]);

        /* Now, rhmm = root channel for w.  Go on to remaining phones */
        if (dict_pronlen(dict, w) == 2) {
            /* Next phone is the last; not kept in tree; add w to penult_phn_wid set */
            if ((j = rhmm->penult_phn_wid) < 0)
                rhmm->penult_phn_wid = w;
            else {
                for (; fts->homophone_set[j] >= 0; j = fts->homophone_set[j]);
                fts->homophone_set[j] = w;
            }
        }
        else {
            /* Add remaining phones, except the last, to tree */
            ph = dict2pid_internal(d2p, w, 1);
            tmatid = bin_mdef_pid2tmatid(search_acmod(fts)->mdef, dict_pron(dict, w, 1));
            hmm = rhmm->next;
            if (hmm == NULL) {
                rhmm->next = hmm = listelem_malloc(fts->chan_alloc);
                init_nonroot_chan(fts, hmm, ph, dict_pron(dict, w, 1), tmatid);
                fts->n_nonroot_chan++;
            }
            else {
                nonroot_node_t *prev_hmm = NULL;

                for (; hmm && (hmm_nonmpx_ssid(&hmm->hmm) != ph); hmm = hmm->alt)
                    prev_hmm = hmm;
                if (!hmm) {     /* thanks, rkm! */
                    prev_hmm->alt = hmm = listelem_malloc(fts->chan_alloc);
                    init_nonroot_chan(fts, hmm, ph, dict_pron(dict, w, 1), tmatid);
                    fts->n_nonroot_chan++;
                }
            }
            E_DEBUG(3,("phone %s = %d\n",
                       bin_mdef_ciphone_str(search_acmod(fts)->mdef,
                                            dict_second_phone(dict, w)), ph));
            for (p = 2; p < dict_pronlen(dict, w) - 1; p++) {
                ph = dict2pid_internal(d2p, w, p);
                tmatid = bin_mdef_pid2tmatid(search_acmod(fts)->mdef, dict_pron(dict, w, p));
                if (!hmm->next) {
                    hmm->next = listelem_malloc(fts->chan_alloc);
                    hmm = hmm->next;
                    init_nonroot_chan(fts, hmm, ph, dict_pron(dict, w, p), tmatid);
                    fts->n_nonroot_chan++;
                }
                else {
                    nonroot_node_t *prev_hmm = NULL;

                    for (hmm = hmm->next; hmm && (hmm_nonmpx_ssid(&hmm->hmm) != ph);
                         hmm = hmm->alt)
                        prev_hmm = hmm;
                    if (!hmm) { /* thanks, rkm! */
                        prev_hmm->alt = hmm = listelem_malloc(fts->chan_alloc);
                        init_nonroot_chan(fts, hmm, ph, dict_pron(dict, w, p), tmatid);
                        fts->n_nonroot_chan++;
                    }
                }
                E_DEBUG(3,("phone %s = %d\n",
                           bin_mdef_ciphone_str(search_acmod(fts)->mdef,
                                                dict_pron(dict, w, p)), ph));
            }

            /* All but last phone of w in tree; add w to hmm->info.penult_phn_wid set */
            if ((j = hmm->info.penult_phn_wid) < 0)
                hmm->info.penult_phn_wid = w;
            else {
                for (; fts->homophone_set[j] >= 0; j = fts->homophone_set[j]);
                fts->homophone_set[j] = w;
            }
        }
    }

    fts->n_1ph_words = fts->n_1ph_LMwords;

    /* Add filler words to the array of 1ph words. */
    for (w = 0; w < n_words; ++w) {
        /* Skip anything that doesn't actually have a single phone. */
        if (!dict_is_single_phone(dict, w))
            continue;
        /* Also skip "real words" and thifts that are in the LM. */
        if (dict_real_word(dict, w))
            continue;
        if (ngram_model_set_known_wid(fts->lmset, dict_basewid(dict, w)))
            continue;
        E_DEBUG(1,("single_phone_wid[%d] = %s\n",
                   fts->n_1ph_words, dict_wordstr(dict, w)));
        fts->single_phone_wid[fts->n_1ph_words++] = w;
    }

    if (fts->n_nonroot_chan >= fts->max_nonroot_chan) {
        /* Give some room for channels for new words added dynamically at run time */
        fts->max_nonroot_chan = fts->n_nonroot_chan + 128;
        E_INFO("after: max nonroot chan increased to %d\n", fts->max_nonroot_chan);

        /* Free old active channel list array if any and allocate new one */
        if (fts->active_chan_list)
            ckd_free_2d(fts->active_chan_list);
        fts->active_chan_list = ckd_calloc_2d(2, fts->max_nonroot_chan,
                                              sizeof(**fts->active_chan_list));
    }

    E_INFO("after: %d root, %d non-root channels, %d single-phone words\n",
           fts->n_root_chan, fts->n_nonroot_chan, fts->n_1ph_words);
}

static void
reinit_search_subtree(fwdtree_search_t *fts, nonroot_node_t * hmm)
{
    nonroot_node_t *child, *sibling;

    /* First free all children under hmm */
    for (child = hmm->next; child; child = sibling) {
        sibling = child->alt;
        reinit_search_subtree(fts, child);
    }

    /* Now free hmm */
    hmm_deinit(&hmm->hmm);
    listelem_free(fts->chan_alloc, hmm);
}

/*
 * Delete search tree by freeing all interior channels within search tree and
 * restoring root channel state to the init state (i.e., just after init_search_tree()).
 */
static void
reinit_search_tree(fwdtree_search_t *fts)
{
    int32 i;
    nonroot_node_t *hmm, *sibling;

    for (i = 0; i < fts->n_root_chan; i++) {
        hmm = fts->root_chan[i].next;

        while (hmm) {
            sibling = hmm->alt;
            reinit_search_subtree(fts, hmm);
            hmm = sibling;
        }

        fts->root_chan[i].penult_phn_wid = -1;
        fts->root_chan[i].next = NULL;
    }
    fts->n_nonroot_chan = 0;
}

static void
deinit_search_tree(fwdtree_search_t *fts)
{
    int i, w, n_words;

    n_words = search_n_words(fts);
    for (i = 0; i < fts->n_root_chan_alloc; i++) {
        hmm_deinit(&fts->root_chan[i].hmm);
    }
    if (fts->rhmm_1ph) {
        for (i = w = 0; w < n_words; ++w) {
            if (!dict_is_single_phone(search_dict(fts), w))
                continue;
            hmm_deinit(&fts->rhmm_1ph[i].hmm);
            ++i;
        }
        ckd_free(fts->rhmm_1ph);
        fts->rhmm_1ph = NULL;
    }
    fts->n_root_chan = 0;
    fts->n_root_chan_alloc = 0;
    ckd_free(fts->root_chan);
    fts->root_chan = NULL;
    ckd_free(fts->single_phone_wid);
    fts->single_phone_wid = NULL;
    ckd_free(fts->homophone_set);
    fts->homophone_set = NULL;
}

int
fwdtree_search_free(search_t *base)
{
    fwdtree_search_t *fts = (fwdtree_search_t *)base;
    int cf;

    cf = acmod_frame(search_acmod(fts));

    /* Reset non-root channels. */
    reinit_search_tree(fts);
    /* Free the search tree. */
    deinit_search_tree(fts);
    /* Free other stuff. */
    fts->max_nonroot_chan = 0;
    ckd_free_2d(fts->active_chan_list);
    fts->active_chan_list = NULL;
    ckd_free(fts->cand_sf);
    fts->cand_sf = NULL;
    ckd_free(fts->bestbp_rc);
    fts->bestbp_rc = NULL;
    ckd_free(fts->lastphn_cand);
    fts->lastphn_cand = NULL;

    hmm_context_free(fts->hmmctx);
    listelem_alloc_free(fts->chan_alloc);
    listelem_alloc_free(fts->root_chan_alloc);
    ngram_model_free(fts->lmset);

    ckd_free(fts->rcss);
    ckd_free(fts->word_idx);
    ckd_free(fts->word_chan);
    bitvec_free(fts->word_active);
    bptbl_free(fts->bptbl);
    ckd_free_2d(fts->active_word_list);
    ckd_free(fts->last_ltrans);

    return 0;
}

static int
fwdtree_search_start(search_t *base)
{
    fwdtree_search_t *fts = (fwdtree_search_t *)base;
    int32 i, w, n_words;
    root_node_t *rhmm;

    n_words = search_n_words(fts);

    /* Reset trigram cache */
    ngram_model_flush(fts->lmset);

    /* Reset utterance statistics. */
    memset(&fts->st, 0, sizeof(fts->st));

    /* Reset backpointer table. */
    bptbl_reset(fts->bptbl);
    fts->oldest_bp = NO_BP;

    /* Reset output arc buffer. */
    if (search_output_arcs(fts))
        arc_buffer_producer_start_utt(search_output_arcs(fts),
                                      base->uttid);

    /* Reset word lattice. */
    for (i = 0; i < n_words; ++i)
        fts->word_idx[i] = NO_BP;

    /* Reset active HMM and word lists. */
    fts->n_active_chan[0] = fts->n_active_chan[1] = 0;
    fts->n_active_word[0] = fts->n_active_word[1] = 0;

    /* Reset scores. */
    fts->best_score = 0;
    fts->renormalized = 0;

    /* Reset other stuff. */
    for (i = 0; i < n_words; i++)
        fts->last_ltrans[i].sf = -1;

    /* Clear the hypothesis string. */
    ckd_free(base->hyp_str);
    base->hyp_str = NULL;

    /* Reset the permanently allocated single-phone words, since they
     * may have junk left over in them from FWDFLAT. */
    for (i = 0; i < fts->n_1ph_words; i++) {
        w = fts->single_phone_wid[i];
        rhmm = (root_node_t *) fts->word_chan[w];
        hmm_clear(&rhmm->hmm);
    }

    /* Start search with <s>; word_chan[<s>] is permanently allocated */
    rhmm = (root_node_t *) fts->word_chan[dict_startwid(search_dict(fts))];
    hmm_clear(&rhmm->hmm);
    hmm_enter(&rhmm->hmm, 0, NO_BP, 0);

    search_call_event(base, SEARCH_START_UTT, base->acmod->output_frame);

    return 0;
}

/*
 * Mark the active senones for all senones belonging to channels that are active in the
 * current frame, and update the oldest backpointer entry.
 */
#define MPX_BITVEC_SET(a,h,i)                                   \
    if (hmm_mpx_ssid(h,i) != BAD_SSID) {                        \
        bitvec_set((a)->senone_active_vec, hmm_mpx_senid(h,i)); \
        if (hmm_history(h,i) < fts->oldest_bp)                  \
            fts->oldest_bp = hmm_history(h,i);                  \
    }

static void
fwdtree_activate_mpx_hmm(fwdtree_search_t *fts, hmm_t *hmm)
{
    acmod_t *acmod = search_acmod(fts);
    int i;

    switch (hmm_n_emit_state(hmm)) {
    case 5:
        MPX_BITVEC_SET(acmod, hmm, 4);
        MPX_BITVEC_SET(acmod, hmm, 3);
    case 3:
        MPX_BITVEC_SET(acmod, hmm, 2);
        MPX_BITVEC_SET(acmod, hmm, 1);
        MPX_BITVEC_SET(acmod, hmm, 0);
        break;
    default:
        for (i = 0; i < hmm_n_emit_state(hmm); ++i) {
            MPX_BITVEC_SET(acmod, hmm, i);
        }
    }
    if (hmm_out_score(hmm) BETTER_THAN WORST_SCORE)
        if (hmm_out_history(hmm) < fts->oldest_bp)
            fts->oldest_bp = hmm_out_history(hmm);
}

#define NONMPX_BITVEC_SET(a,h,i)                                        \
    bitvec_set((a)->senone_active_vec,                                  \
               hmm_nonmpx_senid(h,i));                                  \
    if (hmm_history(h,i) < fts->oldest_bp)                              \
        fts->oldest_bp = hmm_history(h,i);                              \

static void
fwdtree_activate_hmm(fwdtree_search_t *fts, hmm_t *hmm)
{
    acmod_t *acmod = search_acmod(fts);
    int i;

    switch (hmm_n_emit_state(hmm)) {
    case 5:
        NONMPX_BITVEC_SET(acmod, hmm, 4);
        NONMPX_BITVEC_SET(acmod, hmm, 3);
    case 3:
        NONMPX_BITVEC_SET(acmod, hmm, 2);
        NONMPX_BITVEC_SET(acmod, hmm, 1);
        NONMPX_BITVEC_SET(acmod, hmm, 0);
        break;
    default:
        for (i = 0; i < hmm_n_emit_state(hmm); ++i) {
            NONMPX_BITVEC_SET(acmod, hmm, i);
        }
    }
    if (hmm_out_score(hmm) BETTER_THAN WORST_SCORE)
        if (hmm_out_history(hmm) < fts->oldest_bp)
            fts->oldest_bp = hmm_out_history(hmm);
}

static void
compute_sen_active(fwdtree_search_t *fts, int frame_idx)
{
    root_node_t *rhmm;
    nonroot_node_t *hmm, **acl;
    int32 i, w, *awl;

    acmod_clear_active(search_acmod(fts));
    fts->oldest_bp = bptbl_end_idx(fts->bptbl);

    /* Flag active senones for root channels */
    for (i = fts->n_root_chan, rhmm = fts->root_chan; i > 0; --i, rhmm++) {
        if (hmm_frame(&rhmm->hmm) == frame_idx)
            fwdtree_activate_mpx_hmm(fts, &rhmm->hmm);
    }

    /* Flag active senones for nonroot channels in HMM tree */
    i = fts->n_active_chan[frame_idx & 0x1];
    acl = fts->active_chan_list[frame_idx & 0x1];
    for (hmm = *(acl++); i > 0; --i, hmm = *(acl++))
        fwdtree_activate_hmm(fts, &hmm->hmm);

    /* Flag active senones for individual word channels */
    /* These are the only ones that can actually generate new backpointers. */
    i = fts->n_active_word[frame_idx & 0x1];
    awl = fts->active_word_list[frame_idx & 0x1];
    for (w = *(awl++); i > 0; --i, w = *(awl++)) {
        for (hmm = fts->word_chan[w]; hmm; hmm = hmm->next)
            fwdtree_activate_hmm(fts, &hmm->hmm);
    }
    for (i = 0; i < fts->n_1ph_words; i++) {
        w = fts->single_phone_wid[i];
        rhmm = (root_node_t *) fts->word_chan[w];

        if (hmm_frame(&rhmm->hmm) == frame_idx)
            fwdtree_activate_mpx_hmm(fts, &rhmm->hmm);
    }
}

static void
renormalize_scores(fwdtree_search_t *fts, int frame_idx, int32 norm)
{
    root_node_t *rhmm;
    nonroot_node_t *hmm, **acl;
    int32 i, w, *awl;

    /* Renormalize root channels */
    for (i = fts->n_root_chan, rhmm = fts->root_chan; i > 0; --i, rhmm++) {
        if (hmm_frame(&rhmm->hmm) == frame_idx) {
            hmm_normalize(&rhmm->hmm, norm);
        }
    }

    /* Renormalize nonroot channels in HMM tree */
    i = fts->n_active_chan[frame_idx & 0x1];
    acl = fts->active_chan_list[frame_idx & 0x1];
    for (hmm = *(acl++); i > 0; --i, hmm = *(acl++)) {
        hmm_normalize(&hmm->hmm, norm);
    }

    /* Renormalize individual word channels */
    i = fts->n_active_word[frame_idx & 0x1];
    awl = fts->active_word_list[frame_idx & 0x1];
    for (w = *(awl++); i > 0; --i, w = *(awl++)) {
        for (hmm = fts->word_chan[w]; hmm; hmm = hmm->next) {
            hmm_normalize(&hmm->hmm, norm);
        }
    }
    for (i = 0; i < fts->n_1ph_words; i++) {
        w = fts->single_phone_wid[i];
        rhmm = (root_node_t *) fts->word_chan[w];
        if (hmm_frame(&rhmm->hmm) == frame_idx) {
            hmm_normalize(&rhmm->hmm, norm);
        }
    }

    fts->renormalized = TRUE;
}

static int32
eval_root_chan(fwdtree_search_t *fts, int frame_idx)
{
    root_node_t *rhmm;
    int32 i, bestscore;

    bestscore = WORST_SCORE;
    for (i = fts->n_root_chan, rhmm = fts->root_chan; i > 0; --i, rhmm++) {
        if (hmm_frame(&rhmm->hmm) == frame_idx) {
            int32 score = chan_v_eval(rhmm);
            if (score BETTER_THAN bestscore)
                bestscore = score;
            ++fts->st.n_root_chan_eval;
        }
    }
    return (bestscore);
}

static int32
eval_nonroot_chan(fwdtree_search_t *fts, int frame_idx)
{
    nonroot_node_t *hmm, **acl;
    int32 i, bestscore;

    i = fts->n_active_chan[frame_idx & 0x1];
    acl = fts->active_chan_list[frame_idx & 0x1];
    bestscore = WORST_SCORE;
    fts->st.n_nonroot_chan_eval += i;

    for (hmm = *(acl++); i > 0; --i, hmm = *(acl++)) {
        int32 score = chan_v_eval(hmm);
        assert(hmm_frame(&hmm->hmm) == frame_idx);
        if (score BETTER_THAN bestscore)
            bestscore = score;
    }

    return bestscore;
}

static int32
eval_word_chan(fwdtree_search_t *fts, int frame_idx)
{
    root_node_t *rhmm;
    nonroot_node_t *hmm;
    int32 i, w, bestscore, *awl, j, k;

    k = 0;
    bestscore = WORST_SCORE;
    awl = fts->active_word_list[frame_idx & 0x1];

    i = fts->n_active_word[frame_idx & 0x1];
    for (w = *(awl++); i > 0; --i, w = *(awl++)) {
        assert(bitvec_is_set(fts->word_active, w));
        bitvec_clear(fts->word_active, w);
        assert(fts->word_chan[w] != NULL);

        for (hmm = fts->word_chan[w]; hmm; hmm = hmm->next) {
            int32 score;

            assert(hmm_frame(&hmm->hmm) == frame_idx);
            score = chan_v_eval(hmm);

            if (score BETTER_THAN bestscore)
                bestscore = score;

            k++;
        }
    }

    /* Similarly for statically allocated single-phone words */
    j = 0;
    for (i = 0; i < fts->n_1ph_words; i++) {
        int32 score;

        w = fts->single_phone_wid[i];
        rhmm = (root_node_t *) fts->word_chan[w];
        if (hmm_frame(&rhmm->hmm) < frame_idx)
            continue;

#if 0
        E_INFOCONT("1ph word chan %d: %s\n", w,
                   dict_wordstr(search_dict(fts), w));
        score = hmm_dump_vit_eval(&(rhmm)->hmm, stderr);
#else
        score = chan_v_eval(rhmm);
#endif
        if (score BETTER_THAN bestscore && w != search_finish_wid(fts))
            bestscore = score;

        j++;
    }

    fts->st.n_last_chan_eval += k + j;
    fts->st.n_nonroot_chan_eval += k + j;
    fts->st.n_word_lastchan_eval +=
        fts->n_active_word[frame_idx & 0x1] + j;

    return bestscore;
}

static int32
evaluate_channels(fwdtree_search_t *fts, int16 const *senone_scores, int frame_idx)
{
    int32 bs;

    hmm_context_set_senscore(fts->hmmctx, senone_scores);
    fts->best_score = eval_root_chan(fts, frame_idx);
    if ((bs = eval_nonroot_chan(fts, frame_idx)) BETTER_THAN fts->best_score)
        fts->best_score = bs;
    if ((bs = eval_word_chan(fts, frame_idx)) BETTER_THAN fts->best_score)
        fts->best_score = bs;
    fts->last_phone_best_score = bs;

    return fts->best_score;
}

/*
 * Prune currently active root channels for next frame.  Also, perform exit
 * transitions out of them and activate successors.
 * score[] of pruned root chan set to WORST_SCORE elsewhere.
 */
static void
prune_root_chan(fwdtree_search_t *fts, int frame_idx)
{
    root_node_t *rhmm;
    nonroot_node_t *hmm;
    int32 i, nf, w;
    int32 thresh, newphone_thresh, lastphn_thresh, newphone_score;
    nonroot_node_t **nacl;              /* next active list */
    lastphn_cand_t *candp;

    nf = frame_idx + 1;
    thresh = fts->best_score + fts->dynamic_beam;
    newphone_thresh = fts->best_score + fts->pbeam;
    lastphn_thresh = fts->best_score + fts->lpbeam;
    nacl = fts->active_chan_list[nf & 0x1];

    for (i = 0, rhmm = fts->root_chan; i < fts->n_root_chan; i++, rhmm++) {
        /* First check if this channel was active in current frame */
        if (hmm_frame(&rhmm->hmm) < frame_idx)
            continue;

        if (hmm_bestscore(&rhmm->hmm) BETTER_THAN thresh) {
            hmm_frame(&rhmm->hmm) = nf;  /* rhmm will be active in next frame */
            /* transitions out of this root channel */
            /* transition to all next-level channels in the HMM tree */
            newphone_score = hmm_out_score(&rhmm->hmm) + fts->pip;
            if (newphone_score BETTER_THAN newphone_thresh) {
                for (hmm = rhmm->next; hmm; hmm = hmm->alt) {
                    if ((hmm_frame(&hmm->hmm) < frame_idx)
                        || (newphone_score BETTER_THAN hmm_in_score(&hmm->hmm))) {
                        hmm_enter(&hmm->hmm, newphone_score,
                                  hmm_out_history(&rhmm->hmm), nf);
                        *(nacl++) = hmm;
                    }
                }
            }

            /*
             * Transition to last phone of all words for which this is the
             * penultimate phone (the last phones may need multiple right contexts).
             * Remember to remove the temporary newword_penalty.
             */
            if (newphone_score BETTER_THAN lastphn_thresh) {
                for (w = rhmm->penult_phn_wid; w >= 0;
                     w = fts->homophone_set[w]) {
                    candp = fts->lastphn_cand + fts->n_lastphn_cand;
                    fts->n_lastphn_cand++;
                    candp->wid = w;
                    candp->score = newphone_score - fts->nwpen;
                    candp->bp = hmm_out_history(&rhmm->hmm);
                }
            }
        }
    }
    fts->n_active_chan[nf & 0x1] = nacl - fts->active_chan_list[nf & 0x1];
}

/*
 * Prune currently active nonroot channels in HMM tree for next frame.  Also, perform
 * exit transitions out of such channels and activate successors.
 */
static void
prune_nonroot_chan(fwdtree_search_t *fts, int frame_idx)
{
    nonroot_node_t *hmm, *nexthmm;
    int32 nf, w, i;
    int32 thresh, newphone_thresh, lastphn_thresh, newphone_score;
    nonroot_node_t **acl, **nacl;       /* active list, next active list */
    lastphn_cand_t *candp;

    nf = frame_idx + 1;

    thresh = fts->best_score + fts->dynamic_beam;
    newphone_thresh = fts->best_score + fts->pbeam;
    lastphn_thresh = fts->best_score + fts->lpbeam;

    acl = fts->active_chan_list[frame_idx & 0x1];   /* currently active HMMs in tree */
    nacl = fts->active_chan_list[nf & 0x1] + fts->n_active_chan[nf & 0x1];

    for (i = fts->n_active_chan[frame_idx & 0x1], hmm = *(acl++); i > 0;
         --i, hmm = *(acl++)) {
        assert(hmm_frame(&hmm->hmm) >= frame_idx);

        if (hmm_bestscore(&hmm->hmm) BETTER_THAN thresh) {
            /* retain this channel in next frame */
            if (hmm_frame(&hmm->hmm) != nf) {
                hmm_frame(&hmm->hmm) = nf;
                *(nacl++) = hmm;
            }

            /* transition to all next-level channel in the HMM tree */
            newphone_score = hmm_out_score(&hmm->hmm) + fts->pip;
            if (newphone_score BETTER_THAN newphone_thresh) {
                for (nexthmm = hmm->next; nexthmm; nexthmm = nexthmm->alt) {
                    if (((hmm_frame(&nexthmm->hmm) < frame_idx)
                         || (newphone_score
                             BETTER_THAN hmm_in_score(&nexthmm->hmm)))) {
                        if (hmm_frame(&nexthmm->hmm) != nf) {
                            /* Keep this HMM on the active list */
                            *(nacl++) = nexthmm;
                        }
                        hmm_enter(&nexthmm->hmm, newphone_score,
                                  hmm_out_history(&hmm->hmm), nf);
                    }
                }
            }

            /*
             * Transition to last phone of all words for which this is the
             * penultimate phone (the last phones may need multiple right contexts).
             * Remember to remove the temporary newword_penalty.
             */
            if (newphone_score BETTER_THAN lastphn_thresh) {
                for (w = hmm->info.penult_phn_wid; w >= 0;
                     w = fts->homophone_set[w]) {
                    candp = fts->lastphn_cand + fts->n_lastphn_cand;
                    fts->n_lastphn_cand++;
                    candp->wid = w;
                    candp->score =
                        newphone_score - fts->nwpen;
                    candp->bp = hmm_out_history(&hmm->hmm);
                }
            }
        }
        else if (hmm_frame(&hmm->hmm) != nf) {
            hmm_clear(&hmm->hmm);
        }
    }
    fts->n_active_chan[nf & 0x1] = nacl - fts->active_chan_list[nf & 0x1];
}

/*
 * Execute the transition into the last phone for all candidates words emerging from
 * the HMM tree.  Attach LM scores to such transitions.
 * (Executed after pruning root and non-root, but before pruning word-chan.)
 */
static void
last_phone_transition(fwdtree_search_t *fts, int frame_idx)
{
    int32 i, j, k, nf, bp, w;
    lastphn_cand_t *candp;
    int32 *nawl;
    int32 thresh;
    int32 bestscore, dscr;
    nonroot_node_t *hmm;
    int32 n_cand_sf = 0;

    nf = frame_idx + 1;
    nawl = fts->active_word_list[nf & 0x1];
    fts->st.n_lastphn_cand_utt += fts->n_lastphn_cand;

    /* For each candidate word (entering its last phone) */
    /* If best LM score and bp for candidate known use it, else sort cands by startfrm */
    for (i = 0, candp = fts->lastphn_cand; i < fts->n_lastphn_cand; i++, candp++) {
        int32 start_score;
        bp_t bpe;

        /* This can happen if recognition fails. */
        if (candp->bp == NO_BP)
            continue;
        /* Backpointer entry for it. */
        bptbl_get_bp(fts->bptbl, candp->bp, &bpe);

        /* Subtract starting score for candidate, leave it with only word score */
        start_score = fwdtree_search_exit_score
            (fts, candp->bp, bpe.last_phone, bpe.last2_phone,
             dict_first_phone(search_dict(fts), candp->wid));
        assert(start_score BETTER_THAN WORST_SCORE);
        candp->score -= start_score;
        E_DEBUG(4, ("candp->score %d\n", candp->score));

        /*
         * If this candidate not occurred in an earlier frame, prepare for finding
         * best transition score into last phone; sort by start frame.
         */
        /* i.e. if we don't have an entry in last_ltrans for this
         * <word,sf>, then create one */
        if (fts->last_ltrans[candp->wid].sf != bpe.frame + 1) {
            /* Look for an entry in cand_sf matching the backpointer
             * for this candidate. */
            for (j = 0; j < n_cand_sf; j++) {
                if (fts->cand_sf[j].bp_ef == bpe.frame)
                    break;
            }
            /* Oh, we found one, so chain onto it. */
            if (j < n_cand_sf)
                candp->next = fts->cand_sf[j].cand;
            else {
                /* Nope, let's make a new one, allocating cand_sf if necessary. */
                if (n_cand_sf >= fts->cand_sf_alloc) {
                    if (fts->cand_sf_alloc == 0) {
                        fts->cand_sf =
                            ckd_calloc(CAND_SF_ALLOCSIZE,
                                       sizeof(*fts->cand_sf));
                        fts->cand_sf_alloc = CAND_SF_ALLOCSIZE;
                    }
                    else {
                        fts->cand_sf_alloc += CAND_SF_ALLOCSIZE;
                        fts->cand_sf = ckd_realloc(fts->cand_sf,
                                                   fts->cand_sf_alloc
                                                   * sizeof(*fts->cand_sf));
                        E_INFO("cand_sf[] increased to %d entries\n",
                               fts->cand_sf_alloc);
                    }
                }

                /* Use the newly created cand_sf. */
                j = n_cand_sf++;
                candp->next = -1; /* End of the chain. */
                fts->cand_sf[j].bp_ef = bpe.frame;
            }
            /* Update it to point to this candidate. */
            fts->cand_sf[j].cand = i;

            fts->last_ltrans[candp->wid].dscr = WORST_SCORE;
            fts->last_ltrans[candp->wid].sf = bpe.frame + 1;
        }
    }

    /* Compute best LM score and bp for new cands entered in the sorted lists above */
    for (i = 0; i < n_cand_sf; i++) {
        int32 bpend;
        /* This is the last frame that contains end-sorted
         * backpointers.  Luckily by the definition of window_sf it's
         * not possible to have any candidates that fail this
         * assertion. */
        assert(fts->cand_sf[i].bp_ef >= fts->bptbl->active_fr);
        /* For the i-th unique end frame... */
        bp = bptbl_ef_idx(fts->bptbl, fts->cand_sf[i].bp_ef);
        bpend = bptbl_ef_idx(fts->bptbl, fts->cand_sf[i].bp_ef + 1);
        for (; bp < bpend; bp++) {
            bp_t bpe;
            bptbl_get_bp(fts->bptbl, bp, &bpe);
            if (!bpe.valid)
                continue;
            /* For each candidate at the start frame find bp->cand transition-score */
            for (j = fts->cand_sf[i].cand; j >= 0; j = candp->next) {
                int32 n_used;
                candp = &(fts->lastphn_cand[j]);
                dscr = 
                    fwdtree_search_exit_score
                    (fts, bp, bpe.last_phone, bpe.last2_phone,
                     dict_first_phone(search_dict(fts), candp->wid));
                if (dscr BETTER_THAN WORST_SCORE) {
                    assert(!dict_filler_word(search_dict(fts), candp->wid));
                    dscr += ngram_tg_score(fts->lmset,
                                           dict_basewid(search_dict(fts), candp->wid),
                                           bpe.real_wid,
                                           bpe.prev_real_wid, &n_used)>>SENSCR_SHIFT;
                }
                if (dscr BETTER_THAN fts->last_ltrans[candp->wid].dscr) {
                    fts->last_ltrans[candp->wid].dscr = dscr;
                    fts->last_ltrans[candp->wid].bp = bp;
                }
            }
        }
    }

    /* Update best transitions for all candidates; also update best lastphone score */
    bestscore = fts->last_phone_best_score;
    for (i = 0, candp = fts->lastphn_cand; i < fts->n_lastphn_cand; i++, candp++) {
        candp->score += fts->last_ltrans[candp->wid].dscr;
        candp->bp = fts->last_ltrans[candp->wid].bp;

        if (candp->score BETTER_THAN bestscore)
            bestscore = candp->score;
    }
    fts->last_phone_best_score = bestscore;

    /* At this pt, we know the best entry score (with LM component) for all candidates */
    thresh = bestscore + fts->lponlybeam;
    for (i = fts->n_lastphn_cand, candp = fts->lastphn_cand; i > 0; --i, candp++) {
        if (candp->score BETTER_THAN thresh) {
            w = candp->wid;

            fwdtree_search_alloc_all_rc(fts, w);

            k = 0;
            for (hmm = fts->word_chan[w]; hmm; hmm = hmm->next) {
                if ((hmm_frame(&hmm->hmm) < frame_idx)
                    || (candp->score BETTER_THAN hmm_in_score(&hmm->hmm))) {
                    assert(hmm_frame(&hmm->hmm) != nf);
                    hmm_enter(&hmm->hmm,
                              candp->score, candp->bp, nf);
                    k++;
                }
            }
            if (k > 0) {
                assert(bitvec_is_clear(fts->word_active, w));
                assert(!dict_is_single_phone(search_dict(fts), w));
                *(nawl++) = w;
                bitvec_set(fts->word_active, w);
            }
        }
    }
    fts->n_active_word[nf & 0x1] = nawl - fts->active_word_list[nf & 0x1];
}

/*
 * Prune currently active word channels for next frame.  Also, perform exit
 * transitions out of such channels and active successors.
 */
static void
prune_word_chan(fwdtree_search_t *fts, int frame_idx)
{
    root_node_t *rhmm;
    nonroot_node_t *hmm, *thmm;
    nonroot_node_t **phmmp;             /* previous HMM-pointer */
    int32 nf, w, i, k;
    int32 newword_thresh, lastphn_thresh;
    int32 *awl, *nawl;

    nf = frame_idx + 1;
    newword_thresh = fts->last_phone_best_score + fts->wbeam;
    lastphn_thresh = fts->last_phone_best_score + fts->lponlybeam;

    awl = fts->active_word_list[frame_idx & 0x1];
    nawl = fts->active_word_list[nf & 0x1] + fts->n_active_word[nf & 0x1];

    /* Dynamically allocated last channels of multi-phone words */
    for (i = fts->n_active_word[frame_idx & 0x1], w = *(awl++); i > 0;
         --i, w = *(awl++)) {
        k = 0;
        phmmp = &(fts->word_chan[w]);
        for (hmm = fts->word_chan[w]; hmm; hmm = thmm) {
            assert(hmm_frame(&hmm->hmm) >= frame_idx);

            thmm = hmm->next;
            if (hmm_bestscore(&hmm->hmm) BETTER_THAN lastphn_thresh) {
                /* retain this channel in next frame */
                hmm_frame(&hmm->hmm) = nf;
                k++;
                phmmp = &(hmm->next);

                /* Could if ((! skip_alt_frm) || (frame_idx & 0x1)) the following */
                if (hmm_out_score(&hmm->hmm) BETTER_THAN newword_thresh) {
                    /* can exit channel and recognize word */
                    fwdtree_search_save_bp(fts, frame_idx, w,
                                 hmm_out_score(&hmm->hmm),
                                 hmm_out_history(&hmm->hmm),
                                 hmm->info.rc_id);
                }
            }
            else if (hmm_frame(&hmm->hmm) == nf) {
                phmmp = &(hmm->next);
            }
            else {
                hmm_deinit(&hmm->hmm);
                listelem_free(fts->chan_alloc, hmm);
                *phmmp = thmm;
            }
        }
        if ((k > 0) && (bitvec_is_clear(fts->word_active, w))) {
            assert(!dict_is_single_phone(search_dict(fts), w));
            *(nawl++) = w;
            bitvec_set(fts->word_active, w);
        }
    }
    fts->n_active_word[nf & 0x1] = nawl - fts->active_word_list[nf & 0x1];

    /*
     * Prune permanently allocated single-phone channels.
     * NOTES: score[] of pruned channels set to WORST_SCORE elsewhere.
     */
    for (i = 0; i < fts->n_1ph_words; i++) {
        w = fts->single_phone_wid[i];
        rhmm = (root_node_t *) fts->word_chan[w];
        E_DEBUG(3,("Single phone word %s frame %d score %d thresh %d outscore %d nwthresh %d\n",
                   dict_wordstr(search_dict(fts),w),
                   hmm_frame(&rhmm->hmm), hmm_bestscore(&rhmm->hmm),
                   lastphn_thresh, hmm_out_score(&rhmm->hmm), newword_thresh));
        if (hmm_frame(&rhmm->hmm) < frame_idx)
            continue;
        if (hmm_bestscore(&rhmm->hmm) BETTER_THAN lastphn_thresh) {
            hmm_frame(&rhmm->hmm) = nf;

            /* Could if ((! skip_alt_frm) || (frame_idx & 0x1)) the following */
            if (hmm_out_score(&rhmm->hmm) BETTER_THAN newword_thresh) {
                E_DEBUG(4,("Exiting single phone word %s with %d > %d, %d\n",
                           dict_wordstr(search_dict(fts),w),
                           hmm_out_score(&rhmm->hmm),
                           lastphn_thresh, newword_thresh));
                fwdtree_search_save_bp(fts, frame_idx, w,
                             hmm_out_score(&rhmm->hmm),
                             hmm_out_history(&rhmm->hmm), 0);
            }
        }
    }
}

static void
prune_channels(fwdtree_search_t *fts, int frame_idx)
{
    /* Clear last phone candidate list. */
    fts->n_lastphn_cand = 0;
    /* Set the dynamic beam based on maxhmmpf here. */
    fts->dynamic_beam = fts->beam;
    if (fts->maxhmmpf != -1
        && fts->st.n_root_chan_eval + fts->st.n_nonroot_chan_eval > fts->maxhmmpf) {
        /* Build a histogram to approximately prune them. */
        int32 bins[256], bw, nhmms, i;
        root_node_t *rhmm;
        nonroot_node_t **acl, *hmm;

        /* Bins go from zero (best score) to edge of beam. */
        bw = -fts->beam / 256;
        memset(bins, 0, sizeof(bins));
        /* For each active root channel. */
        for (i = 0, rhmm = fts->root_chan; i < fts->n_root_chan; i++, rhmm++) {
            int32 b;

            /* Put it in a bin according to its bestscore. */
            b = (fts->best_score - hmm_bestscore(&rhmm->hmm)) / bw;
            if (b >= 256)
                b = 255;
            ++bins[b];
        }
        /* For each active non-root channel. */
        acl = fts->active_chan_list[frame_idx & 0x1];       /* currently active HMMs in tree */
        for (i = fts->n_active_chan[frame_idx & 0x1], hmm = *(acl++);
             i > 0; --i, hmm = *(acl++)) {
            int32 b;

            /* Put it in a bin according to its bestscore. */
            b = (fts->best_score - hmm_bestscore(&hmm->hmm)) / bw;
            if (b >= 256)
                b = 255;
            ++bins[b];
        }
        /* Walk down the bins to find the new beam. */
        for (i = nhmms = 0; i < 256; ++i) {
            nhmms += bins[i];
            if (nhmms > fts->maxhmmpf)
                break;
        }
        fts->dynamic_beam = -(i * bw);
    }

    prune_root_chan(fts, frame_idx);
    prune_nonroot_chan(fts, frame_idx);
    last_phone_transition(fts, frame_idx);
    prune_word_chan(fts, frame_idx);
}

static int
fwdtree_maxwpf(fwdtree_search_t *fts, int frame_idx)
{
    int32 n, maxwpf;
    int32 bestscr, worstscr;
    bpidx_t start, end;
    bpidx_t bp, bestbp, worstbp;
    bptbl_t *bptbl;

    /* Don't prune if no pruing. */
    maxwpf = fts->maxwpf;
    bptbl = fts->bptbl;
    if (maxwpf == -1)
        return 0;

    /* Get bps. */
    start = bptbl_ef_idx(bptbl, frame_idx);
    end = bptbl_end_idx(bptbl);
    /* Allow only one filler word exit (the best) per frame */
    bestscr = (int32) 0x80000000;
    bestbp = NO_BP;
    n = 0;
    for (bp = start; bp != end; ++bp) {
        bp_t bpe;
        bptbl_get_bp(bptbl, bp, &bpe);
        if (dict_filler_word(bptbl->d2p->dict, bpe.wid)) {
            if (bpe.score BETTER_THAN bestscr) {
                bestscr = bpe.score;
                bestbp = bp;
            }
            /* FIXME: Should have a shortcut to invalidate. */
            bpe.valid = FALSE; /* Flag to indicate invalidation */
            bptbl_set_bp(bptbl, bp, &bpe);
            fts->word_idx[bpe.wid] = NO_BP;
            n++;                /* No. of filler words */
        }
    }
    /* Restore bestbpe to valid state */
    if (bestbp != NO_BP) {
        bp_t bestbpe;
        bptbl_get_bp(bptbl, bestbp, &bestbpe);
        fts->word_idx[bestbpe.wid] = bestbp;
        bestbpe.valid = TRUE;
        bptbl_set_bp(bptbl, bestbp, &bestbpe);
        --n;
    }

    /* Allow up to maxwpf best entries to survive; mark the remaining
     * with valid = 0 */
    /* No. of entries after limiting fillers */
    n = bptbl_ef_count(bptbl, frame_idx) - n;
    for (; n > maxwpf; --n) {
        bp_t worstbpe;
        /* Find worst BPTable entry */
        worstscr = (int32) 0x7fffffff;
        worstbp = NO_BP;
        for (bp = start; bp != end; ++bp) {
            bp_t bpe;
            bptbl_get_bp(bptbl, bp, &bpe);
            if (bpe.valid && (bpe.score WORSE_THAN worstscr)) {
                worstscr = bpe.score;
                worstbp = bp;
            }
        }
        assert(worstbp != NO_BP);
        bptbl_get_bp(bptbl, worstbp, &worstbpe);
        fts->word_idx[worstbpe.wid] = NO_BP;
        worstbpe.valid = FALSE;
        bptbl_set_bp(bptbl, worstbp, &worstbpe);
    }

    return n;
}

static void
word_transition(fwdtree_search_t *fts, int frame_idx)
{
    int32 i, k, w, nf;
    int32 rc;
    int32 thresh, newscore;
    bpidx_t bp;
    root_node_t *rhmm;
    struct bestbp_rc_s *bestbp_rc_ptr;
    dict_t *dict = search_dict(fts);
    dict2pid_t *d2p = search_dict2pid(fts);

    /*
     * Transition to start of new word instances (HMM tree roots); but only if words
     * other than </s> finished here.
     * But, first, find the best starting score for each possible right context phone.
     */
    for (i = bin_mdef_n_ciphone(search_acmod(fts)->mdef) - 1; i >= 0; --i)
        fts->bestbp_rc[i].score = WORST_SCORE;
    k = 0;
    /* Ugh, this is complicated.  Scan all word exits for this frame
     * (they have already been created by prune_word_chan()). */
    for (bp = bptbl_ef_idx(fts->bptbl, frame_idx);
         bp < bptbl_ef_idx(fts->bptbl, frame_idx + 1); bp++) {
        bp_t bpe;
        bptbl_get_bp(fts->bptbl, bp, &bpe);
        fts->word_idx[bpe.wid] = NO_BP;

        /* No transitions from the finish word for obvious reasons. */
        if (bpe.wid == search_finish_wid(fts))
            continue;
        k++;

        /* DICT2PID */
        /* Array of HMM scores corresponding to all the possible right
         * context expansions of the final phone.  It's likely that a
         * lot of these are going to be missing, actually. */
        bptbl_get_rcscores(fts->bptbl, bp, fts->rcss);
        if (bpe.last2_phone == -1) {
            /* No right context expansion. */
            for (rc = 0;
                 rc < bin_mdef_n_ciphone(search_acmod(fts)->mdef); ++rc) {
                if (fts->rcss[0] BETTER_THAN fts->bestbp_rc[rc].score) {
                    E_DEBUG(4,("bestbp_rc[0] = %d lc %d\n",
                               fts->rcss[0], bpe.last_phone));
                    fts->bestbp_rc[rc].score = fts->rcss[0];
                    fts->bestbp_rc[rc].path = bp;
                    fts->bestbp_rc[rc].lc = bpe.last_phone;
                }
            }
        }
        else {
            xwdssid_t *rssid = dict2pid_rssid(d2p, bpe.last_phone, bpe.last2_phone);
            for (rc = 0; rc < bin_mdef_n_ciphone(search_acmod(fts)->mdef); ++rc) {
                if (fts->rcss[rssid->cimap[rc]] BETTER_THAN fts->bestbp_rc[rc].score) {
                    E_DEBUG(4,("bestbp_rc[%d] = %d lc %d\n",
                               rc, fts->rcss[rssid->cimap[rc]], bpe.last_phone));
                    fts->bestbp_rc[rc].score = fts->rcss[rssid->cimap[rc]];
                    fts->bestbp_rc[rc].path = bp;
                    fts->bestbp_rc[rc].lc = bpe.last_phone;
                }
            }
        }
    }
    if (k == 0)
        return;

    nf = frame_idx + 1;
    thresh = fts->best_score + fts->dynamic_beam;
    /*
     * Hypothesize successors to words finished in this frame.
     * Main dictionary, multi-phone words transition to HMM-trees roots.
     */
    for (i = 0, rhmm = fts->root_chan; i < fts->n_root_chan; ++i, ++rhmm) {
        bestbp_rc_ptr = &(fts->bestbp_rc[rhmm->ciphone]);

        newscore = bestbp_rc_ptr->score + fts->nwpen + fts->pip;
        if (newscore BETTER_THAN thresh) {
            if ((hmm_frame(&rhmm->hmm) < frame_idx)
                || (newscore BETTER_THAN hmm_in_score(&rhmm->hmm))) {
                hmm_enter(&rhmm->hmm, newscore,
                          bestbp_rc_ptr->path, nf);
                /* DICT2PID: Another place where mpx ssids are entered. */
                /* Look up the ssid to use when entering this mpx triphone. */
                hmm_mpx_ssid(&rhmm->hmm, 0) =
                    dict2pid_ldiph_lc(d2p, rhmm->ciphone, rhmm->ci2phone, bestbp_rc_ptr->lc);
                assert(hmm_mpx_ssid(&rhmm->hmm, 0) != BAD_SSID);
            }
        }
    }

    /*
     * Single phone words; no right context for these.  Cannot use
     * bestbp_rc as LM scores have to be included.  First find best
     * transition to these words, including language model score, and
     * store it in last_ltrans (which is otherwise not used here).
     */
    for (i = 0; i < fts->n_1ph_LMwords; i++) {
        w = fts->single_phone_wid[i];
        fts->last_ltrans[w].dscr = (int32) 0x80000000;
    }
    for (bp = bptbl_ef_idx(fts->bptbl, frame_idx);
         bp < bptbl_ef_idx(fts->bptbl, frame_idx + 1); bp++) {
        bp_t bpe;

        bptbl_get_bp(fts->bptbl, bp, &bpe);
        if (!bpe.valid)
            continue;

        for (i = 0; i < fts->n_1ph_LMwords; i++) {
            int32 n_used;
            w = fts->single_phone_wid[i];
            newscore = fwdtree_search_exit_score
                (fts, bp, bpe.last_phone, bpe.last2_phone,
                 dict_first_phone(dict, w));
            E_DEBUG(4, ("initial newscore for %s: %d\n",
                        dict_wordstr(dict, w), newscore));
            if (newscore != WORST_SCORE)
                newscore += ngram_tg_score(fts->lmset,
                                           dict_basewid(dict, w),
                                           bpe.real_wid,
                                           bpe.prev_real_wid, &n_used)>>SENSCR_SHIFT;

            /* FIXME: Not sure how WORST_SCORE could be better, but it
             * apparently happens (uh, maybe because WORST_SCORE >
             * (int32)0x80000000?). */
            if (newscore BETTER_THAN fts->last_ltrans[w].dscr) {
                fts->last_ltrans[w].dscr = newscore;
                fts->last_ltrans[w].bp = bp;
            }
        }
    }

    /* Now transition to in-LM single phone words */
    for (i = 0; i < fts->n_1ph_LMwords; i++) {
        w = fts->single_phone_wid[i];
        /* Never transition into the start word (for one thing, it is
           a non-event in the language model.) */
        if (w == dict_startwid(search_dict(fts)))
            continue;
        rhmm = (root_node_t *) fts->word_chan[w];
        newscore = fts->last_ltrans[w].dscr + fts->pip;
        if (newscore BETTER_THAN thresh) {
            bp_t bpe;
            bptbl_get_bp(fts->bptbl, fts->last_ltrans[w].bp, &bpe);
            if ((hmm_frame(&rhmm->hmm) < frame_idx)
                || (newscore BETTER_THAN hmm_in_score(&rhmm->hmm))) {
                hmm_enter(&rhmm->hmm,
                          newscore, fts->last_ltrans[w].bp, nf);
                /* DICT2PID: another place where mpx ssids are entered. */
                /* Look up the ssid to use when entering this mpx triphone. */
                hmm_mpx_ssid(&rhmm->hmm, 0) =
                    dict2pid_ldiph_lc(d2p, rhmm->ciphone, rhmm->ci2phone,
                                      dict_last_phone(dict, bpe.wid));
                assert(hmm_mpx_ssid(&rhmm->hmm, 0) != BAD_SSID);
            }
        }
    }

    /* Remaining words: <sil>, noise words.  No mpx for these! */
    w = search_silence_wid(fts);
    rhmm = (root_node_t *) fts->word_chan[w];
    bestbp_rc_ptr = &(fts->bestbp_rc[search_acmod(fts)->mdef->sil]);
    newscore = bestbp_rc_ptr->score + fts->silpen + fts->pip;
    if (newscore BETTER_THAN thresh) {
        if ((hmm_frame(&rhmm->hmm) < frame_idx)
            || (newscore BETTER_THAN hmm_in_score(&rhmm->hmm))) {
            hmm_enter(&rhmm->hmm,
                      newscore, bestbp_rc_ptr->path, nf);
        }
    }
    for (w = dict_filler_start(dict); w <= dict_filler_end(dict); w++) {
        if (w == search_silence_wid(fts))
            continue;
        /* Never transition into the start word (for one thing, it is
           a non-event in the language model.) */
        if (w == dict_startwid(search_dict(fts)))
            continue;
        rhmm = (root_node_t *) fts->word_chan[w];
        /* If this was not actually a single-phone word, rhmm will be NULL. */
        if (rhmm == NULL)
            continue;
        newscore = bestbp_rc_ptr->score + fts->fillpen + fts->pip;
        if (newscore BETTER_THAN thresh) {
            if ((hmm_frame(&rhmm->hmm) < frame_idx)
                || (newscore BETTER_THAN hmm_in_score(&rhmm->hmm))) {
                hmm_enter(&rhmm->hmm,
                          newscore, bestbp_rc_ptr->path, nf);
            }
        }
    }
}

static void
deactivate_channels(fwdtree_search_t *fts, int frame_idx)
{
    root_node_t *rhmm;
    int i;

    /* Clear score[] of pruned root channels */
    for (i = fts->n_root_chan, rhmm = fts->root_chan; i > 0; --i, rhmm++) {
        if (hmm_frame(&rhmm->hmm) == frame_idx) {
            hmm_clear(&rhmm->hmm);
        }
    }
    /* Clear score[] of pruned single-phone channels */
    for (i = 0; i < fts->n_1ph_words; i++) {
        int32 w = fts->single_phone_wid[i];
        rhmm = (root_node_t *) fts->word_chan[w];
        if (hmm_frame(&rhmm->hmm) == frame_idx) {
            hmm_clear(&rhmm->hmm);
        }
    }
}

static int
update_partials(fwdtree_search_t *fts, int frame_idx)
{
    fts->best_exit = bptbl_find_exit(fts->bptbl, -1);
    if (fts->best_exit != NO_BP) {
        int32 bbp = dict_basewid(search_dict(fts), bptbl_ent(fts->bptbl, fts->best_exit)->wid);
        if (bbp != fts->best_exit_wid) {
            fts->best_exit_wid = bbp;
            search_call_event(search_base(fts), SEARCH_PARTIAL_RESULT, frame_idx);
        }
    }
    return fts->best_exit;
}

static int
fwdtree_search_one_frame(fwdtree_search_t *fts)
{
    acmod_t *acmod = search_acmod(fts);
    int16 const *senscr;
    int fi, frame_idx;

    if ((frame_idx = acmod_consumer_wait(acmod, -1)) < 0) {
        /* Normal end of utterance... */
        if (acmod_eou(acmod))
            return 0;
        /* Or something bad (i.e. cancellation)? */
        return -1;
    }
    E_DEBUG(2,("Searching frame %d\n", frame_idx));
    ptmr_start(&fts->base.t);
    /* Activate our HMMs for the current frame if need be. */
    if (!acmod->compallsen)
        compute_sen_active(fts, frame_idx);

    /* Compute GMM scores for the current frame. */
    if ((senscr = acmod_score(acmod, frame_idx)) == NULL) {
        ptmr_stop(&fts->base.t);
        return 0;
    }
    fts->st.n_senone_active_utt += acmod->n_senone_active;

    /* Mark backpointer table for current frame. */
    fi = bptbl_push_frame(fts->bptbl, fts->oldest_bp);
    assert(fi == frame_idx);

    /* Forward retired backpointers to the arc buffer. */
    if (search_output_arcs(fts))
        arc_buffer_producer_sweep(search_output_arcs(fts),
                                  /* FIXME: should be configurable */
                                  FALSE);

    /* If the best score is equal to or worse than WORST_SCORE,
     * recognition has failed, don't bother to keep trying. */
    if (fts->best_score == WORST_SCORE
        || fts->best_score WORSE_THAN WORST_SCORE) {
        ptmr_stop(&fts->base.t);
        return 0;
    }
    /* Renormalize if necessary */
    if (fts->best_score + (2 * fts->beam) WORSE_THAN WORST_SCORE) {
        E_INFO("Renormalizing Scores at frame %d, best score %d\n",
               frame_idx, fts->best_score);
        renormalize_scores(fts, frame_idx, fts->best_score);
    }

    /* Evaluate HMMs */
    evaluate_channels(fts, senscr, frame_idx);
    /* Prune HMMs and do phone transitions. */
    prune_channels(fts, frame_idx);
    /* Prune and commit new word exits. */
    fwdtree_maxwpf(fts, frame_idx);
    bptbl_commit(fts->bptbl);
    /* Do word transitions. */
    word_transition(fts, frame_idx);
    /* Deactivate pruned HMMs. */
    deactivate_channels(fts, frame_idx);

    /* Release the frame just searched. */
    acmod_consumer_release(acmod, frame_idx);

    /* Update partial results. */
    update_partials(fts, frame_idx);

    /* Return the number of frames processed. */
    ptmr_stop(&fts->base.t);
    return 1;
}

static int
fwdtree_search_finish(search_t *base)
{
    fwdtree_search_t *fts = (fwdtree_search_t *)base;
    int32 i, w, cf, *awl;
    root_node_t *rhmm;
    nonroot_node_t *hmm, **acl;

    ptmr_start(&fts->base.t);

    /* This is the number of frames processed. */
    cf = acmod_frame(search_acmod(fts));

    /* Final result for this thread. */
    search_call_event(base, SEARCH_FINAL_RESULT, cf);

    /* Finalize the backpointer table. */
    bptbl_finalize(fts->bptbl);

    /* Finalize the output arc buffer. */
    if (search_output_arcs(fts))
        arc_buffer_producer_end_utt(search_output_arcs(fts), FALSE);

    /* Deactivate channels lined up for the next frame */
    /* First, root channels of HMM tree */
    for (i = fts->n_root_chan, rhmm = fts->root_chan; i > 0; --i, rhmm++) {
        hmm_clear(&rhmm->hmm);
    }

    /* nonroot channels of HMM tree */
    i = fts->n_active_chan[cf & 0x1];
    acl = fts->active_chan_list[cf & 0x1];
    for (hmm = *(acl++); i > 0; --i, hmm = *(acl++)) {
        hmm_clear(&hmm->hmm);
    }

    /* word channels */
    i = fts->n_active_word[cf & 0x1];
    awl = fts->active_word_list[cf & 0x1];
    for (w = *(awl++); i > 0; --i, w = *(awl++)) {
        /* Don't accidentally free single-phone words! */
        if (dict_is_single_phone(search_dict(fts), w))
            continue;
        bitvec_clear(fts->word_active, w);
        if (fts->word_chan[w] == NULL)
            continue;
        fwdtree_search_free_all_rc(fts, w);
    }
    ptmr_stop(&fts->base.t);

    /* Finalize the input acmod (signals producer) */
    acmod_consumer_end_utt(base->acmod);
    base->total_frames += base->acmod->output_frame;

    search_call_event(base, SEARCH_END_UTT, base->acmod->output_frame);

    /*
     * The previous search code did a postprocessing of the
     * backpointer table here, but we will postpone this until it is
     * absolutely necessary, i.e. when generating a word graph.
     * Likewise we don't actually have to decide what the exit word is
     * until somebody requests a backtrace.
     */

    /* Print out some statistics. */
    if (cf > 0) {
        E_INFO("%8d words recognized in %d frames (%d/fr)\n",
               bptbl_end_idx(fts->bptbl), cf + 1,
               (bptbl_end_idx(fts->bptbl) + (cf >> 1)) / (cf + 1));
        E_INFO("%8d senones evaluated (%d/fr)\n", fts->st.n_senone_active_utt,
               (fts->st.n_senone_active_utt + (cf >> 1)) / (cf + 1));
        E_INFO("%8d channels searched (%d/fr), %d 1st, %d last\n",
               fts->st.n_root_chan_eval + fts->st.n_nonroot_chan_eval,
               (fts->st.n_root_chan_eval + fts->st.n_nonroot_chan_eval) / (cf + 1),
               fts->st.n_root_chan_eval, fts->st.n_last_chan_eval);
        E_INFO("%8d words for which last channels evaluated (%d/fr)\n",
               fts->st.n_word_lastchan_eval,
               fts->st.n_word_lastchan_eval / (cf + 1));
        E_INFO("%8d candidate words for entering last phone (%d/fr)\n",
               fts->st.n_lastphn_cand_utt, fts->st.n_lastphn_cand_utt / (cf + 1));
        E_INFO("time %f wall %.2f xRT\n",
               base->t.t_elapsed,
               base->t.t_elapsed / base->acmod->output_frame
               * cmd_ln_int32_r(base->config, "-frate"));
    }
    /* bptbl_dump(fts->bptbl); */

    return 0;
}

int
fwdtree_search_decode(search_t *base)
{
    fwdtree_search_t *fts = (fwdtree_search_t *)base;
    int nfr, k;

    if (acmod_consumer_start_utt(base->acmod, -1) < 0) {
        if (base->output_arcs)
            arc_buffer_producer_shutdown(base->output_arcs);
        return -1;
    }
    base->uttid = base->acmod->uttid;
    nfr = 0;
    fwdtree_search_start(base);
    while ((k = fwdtree_search_one_frame(fts)) > 0) {
        nfr += k;
    }

    /* Abnormal termination (cancellation, error) */
    if (k < 0) {
        if (base->output_arcs)
            arc_buffer_producer_shutdown(base->output_arcs);
        return k;
    }

    /* This calls arc_buffer_finalize() for us. */
    fwdtree_search_finish(base);
    return nfr;
}

static void
fwdtree_search_alloc_all_rc(fwdtree_search_t *fts, int32 w)
{
    nonroot_node_t *hmm, *thmm;
    xwdssid_t *rssid;
    int32 i, tmatid, ciphone;

    /* DICT2PID */
    /* Get pointer to array of triphones for final diphone. */
    assert(!dict_is_single_phone(search_dict(fts), w));
    ciphone = dict_last_phone(search_dict(fts),w);
    rssid = dict2pid_rssid(search_dict2pid(fts),
                           ciphone,
                           dict_second_last_phone(search_dict(fts),w));
    tmatid = bin_mdef_pid2tmatid(search_acmod(fts)->mdef, ciphone);
    hmm = fts->word_chan[w];
    if ((hmm == NULL) || (hmm_nonmpx_ssid(&hmm->hmm) != rssid->ssid[0])) {
        hmm = listelem_malloc(fts->chan_alloc);
        hmm->next = fts->word_chan[w];
        fts->word_chan[w] = hmm;

        hmm->info.rc_id = 0;
        hmm->ciphone = ciphone;
        hmm_init(fts->hmmctx, &hmm->hmm, FALSE, rssid->ssid[0], tmatid);
        E_DEBUG(3,("allocated rc_id 0 ssid %d ciphone %d lc %d word %s\n",
                   rssid->ssid[0], hmm->ciphone,
                   dict_second_last_phone(search_dict(fts),w),
                   dict_wordstr(search_dict(fts),w)));
    }
    for (i = 1; i < rssid->n_ssid; ++i) {
        if ((hmm->next == NULL) || (hmm_nonmpx_ssid(&hmm->next->hmm) != rssid->ssid[i])) {
            thmm = listelem_malloc(fts->chan_alloc);
            thmm->next = hmm->next;
            hmm->next = thmm;
            hmm = thmm;

            hmm->info.rc_id = i;
            hmm->ciphone = ciphone;
            hmm_init(fts->hmmctx, &hmm->hmm, FALSE, rssid->ssid[i], tmatid);
            E_DEBUG(3,("allocated rc_id %d ssid %d ciphone %d lc %d word %s\n",
                       i, rssid->ssid[i], hmm->ciphone,
                       dict_second_last_phone(search_dict(fts),w),
                       dict_wordstr(search_dict(fts),w)));
        }
        else
            hmm = hmm->next;
    }
}

static int32
fwdtree_search_exit_score(fwdtree_search_t *fts, bpidx_t bp,
                          int last_phone, int last2_phone, int rcphone)
{
    int rcsize;
    /* Get the right context scores for this backpointer. */
    rcsize = bptbl_get_rcscores(fts->bptbl, bp, fts->rcss);
    if (rcsize == 1) {
        /* This indicates a single phone predecessor word, which has
         * no right context scores. */
        assert(fts->rcss[0] != WORST_SCORE);
        return fts->rcss[0];
    }
    else {
        xwdssid_t *rssid;
        /* Find the index for the last diphone of the previous word +
         * the first phone of the current word. */
        rssid = dict2pid_rssid(search_dict2pid(fts),
                               last_phone, last2_phone);
        /* This may be WORST_SCORE, which means that there was no exit
         * with rcphone as right context. */
        assert(rssid->cimap[rcphone] < rcsize);
        return fts->rcss[rssid->cimap[rcphone]];
    }
}

static void
fwdtree_search_free_all_rc(fwdtree_search_t *fts, int32 w)
{
    nonroot_node_t *hmm, *thmm;

    for (hmm = fts->word_chan[w]; hmm; hmm = thmm) {
        thmm = hmm->next;
        hmm_deinit(&hmm->hmm);
        listelem_free(fts->chan_alloc, hmm);
    }
    fts->word_chan[w] = NULL;
}

static void
fwdtree_search_save_bp(fwdtree_search_t *fts, int frame_idx,
                       int32 w, int32 score, int32 path, int32 rc)
{
    int32 bp;

    /* Look for an existing exit for this word in this frame. */
    bp = fts->word_idx[w];
    if (bp != NO_BP) {
        bp_t bpe;

        assert(bp >= bptbl_ef_idx(fts->bptbl, frame_idx));
        bptbl_get_bp(fts->bptbl, bp, &bpe);
        assert(frame_idx == bpe.frame);
        /* Keep only the best scoring language model state (it's
         * actually rare that it will change because the only reason
         * there would be more than one exit for a word is due to
         * multiple right contexts for its end phone) */
        if (bpe.score WORSE_THAN score) {
            bptbl_update_bp(fts->bptbl, bp, rc, path, score);
        }
        /* But do keep track of scores for all right contexts, since
         * we need them to determine the starting path scores for any
         * successors of this word exit. */
        bptbl_set_rcscore(fts->bptbl, bp, rc, score);
    }
    else {
        bpidx_t bpidx = bptbl_enter(fts->bptbl, w, path, score, rc);
        fts->word_idx[w] = bpidx;
    }
}

static int32
fwdtree_search_prob(search_t *base)
{
    /* FIXME: Going to estimate this from partial results in the future. */
    return 0;
}

static char const *
fwdtree_search_hyp(search_t *base, int32 *out_score)
{
    fwdtree_search_t *fts = (fwdtree_search_t *)base;

    ckd_free(base->hyp_str);
    if (bptbl_is_final(fts->bptbl)) {
        base->hyp_str = bptbl_hyp(fts->bptbl, out_score, base->finish_wid);
    }
    else {
        *out_score = fts->best_score;
        base->hyp_str = bptbl_backtrace(fts->bptbl, fts->best_exit);
    }
    return base->hyp_str;
}

static seg_iter_t *
fwdtree_search_seg_iter(search_t *base, int32 *out_score)
{
    fwdtree_search_t *fts = (fwdtree_search_t *) base;
    if (bptbl_is_final(fts->bptbl))
        return bptbl_seg_iter(fts->bptbl, out_score, search_finish_wid(fts));
    else {
        *out_score = fts->best_score;
        return bptbl_seg_backtrace(fts->bptbl, fts->best_exit);
    }
}

static bptbl_t *
fwdtree_search_bptbl(search_t *base)
{
    fwdtree_search_t *fts = (fwdtree_search_t *)base;
    return fts->bptbl;
}

static ngram_model_t *
fwdtree_search_lmset(search_t *base)
{
    fwdtree_search_t *fts = (fwdtree_search_t *)base;
    return fts->lmset;
}
