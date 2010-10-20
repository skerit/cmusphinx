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
 * @file fwdtree_search.h Lexicon tree based Viterbi search.
 */

#ifndef __FWDTREE_SEARCH_H__
#define __FWDTREE_SEARCH_H__

/* SphinxBase headers. */

/* Local headers. */
#include "bptbl.h"

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
    ps_search_t base;
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

/**
 * Initialize fwdtree search.
 */
ps_search_t *fwdtree_search_init(cmd_ln_t *config, acmod_t *acmod,
                                 dict_t *dict, dict2pid_t *d2p);

/**
 * Get the output backpointer table.
 */
bptbl_t *fwdtree_search_bptbl(ps_search_t *base);

/**
 * Get the language model set.
 */
ngram_model_t *fwdtree_search_lmset(ps_search_t *base);

#endif /* __FWDTREE_SEARCH_H__ */
