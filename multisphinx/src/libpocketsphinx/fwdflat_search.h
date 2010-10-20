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
 * @file fwdflat_search.h Flat lexicon based Viterbi search.
 */

#ifndef __FWDFLAT_SEARCH_H__
#define __FWDFLAT_SEARCH_H__

/* SphinxBase headers. */

/* Local headers. */
#include "bptbl.h"
#include "hmm.h"

/**
 * Phone HMM data type.
 *
 * Not the first HMM for words, which multiplex HMMs based on
 * different left contexts.  This structure is used both in the
 * dynamic HMM tree structure and in the per-word last-phone right
 * context fanout.
 */
typedef struct internal_node_s {
    hmm_t hmm;                  /**< Basic HMM structure.  This *must*
                                   be first in the structure because
                                   internal_node_t and first_node_t are
                                   sometimes used interchangeably */
    struct internal_node_s *next;/**< first descendant of this channel; or, in the
				   case of the last phone of a word, the next
				   alternative right context channel */
    int16    ciphone;		/**< ciphone for this node */
    int16    rc_id;		/**< right-context id for last phone of words */
} internal_node_t;

/**
 * Lexical tree node data type for the first phone (first) of each word HMM.
 *
 * Each state may have a different parent static HMM.  Most fields are
 * similar to those in internal_node_t.
 */
typedef struct first_node_s {
    hmm_t hmm;                  /**< Basic HMM structure.  This *must* be first in
                                   the structure because internal_node_t and first_node_t are
                                   sometimes used interchangeably. */
    internal_node_t *next;	/**< first descendant of this channel */

    int16    ciphone;		/**< first ciphone of this node; all words firsted at this
				   node begin with this ciphone */
    int16    ci2phone;		/**< second ciphone of this node; one first HMM for each
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
    ps_search_t base;
    ngram_model_t *lmset;
    hmm_context_t *hmmctx;

    listelem_alloc_t *chan_alloc;      /**< For internal_node_t */
    listelem_alloc_t *root_chan_alloc; /**< For first_node_t */

    /**
     * Backpointer table (temporary storage for active word arcs).
     */
    bptbl_t *bptbl;
    int32 oldest_bp; /**< Oldest bptable entry active in decoding graph. */

    int32 *word_idx;     /**< BPTable index for any word in current frame;
                            cleared before each frame */

    /**
     * Input backpointer table (determines successor words).
     */
    bptbl_t *input_bptbl;
    int16 min_ef_width;   /**< Minimum number of active endpoints. */
    int16 max_sf_win;     /**< Lookahead window for active words */

    uint8 *input_words;     /**< Active words (end frame counts) */
    int16 input_first_sf; /**< First frame of current lookahead window. */
    int16 input_last_sf;  /**< Last start frame found in the input. */

    bpidx_t input_first_bp;
    bpidx_t input_last_bp;

    /**
     * First HMMs for multiple-phone words.
     */
    first_node_t **word_chan;
    bitvec_t *word_active;    /**< array of active flags for all words. */

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

    /* A children's treasury of beam widths. */
    int32 fwdflatbeam;
    int32 fwdflatwbeam;
    int32 fillpen;
    int32 silpen;
    int32 pip;

    /** Are we done? */
    int done;
} fwdflat_search_t;

/**
 * Initialize fwdflat search.
 */
ps_search_t *fwdflat_search_init(cmd_ln_t *config, acmod_t *acmod,
                                 dict_t *dict, dict2pid_t *d2p,
                                 bptbl_t *input_bptbl);

#endif /* __FWDFLAT_SEARCH_H__ */
