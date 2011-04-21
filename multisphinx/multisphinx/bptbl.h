/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008-2010 Carnegie Mellon University.  All rights
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
 * @file bptbl.h Forward search lattice for N-Gram search.
 */

#ifndef __BPTBL_H__
#define __BPTBL_H__

/* SphinxBase headers. */
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/logmath.h>
#include <sphinxbase/ngram_model.h>
#include <sphinxbase/listelem_alloc.h>
#include <sphinxbase/err.h>
#include <sphinxbase/garray.h>
#include <sphinxbase/sbthread.h>
#include <sphinxbase/profile.h>

#include <multisphinx/dict2pid.h>
#include <multisphinx/bptbl.h>

#define NO_BP		-1

/**
 * Numeric backpointer ID
 */
typedef int32 bpidx_t;

/**
 * Back pointer table entry (32 bytes)
 */
typedef struct bp_s {
    int16    frame;		/**< start or end frame */
    uint8    valid;		/**< For absolute pruning */
    uint8    refcnt;            /**< Reference count (number of successors) */
    int32    wid;		/**< Word index */
    int32    bp;		/**< Back Pointer */
    int32    score;		/**< Score (best among all right contexts) */
    int32    s_idx;		/**< Start of BScoreStack for various right contexts*/
    int32    real_wid;		/**< wid of this or latest predecessor real word */
    int32    prev_real_wid;	/**< real word predecessor of real_wid */
    int16    last_phone;        /**< last phone of this word */
    int16    last2_phone;       /**< next-to-last phone of this word */
} bp_t;

/**
 * Back pointer table
 */
typedef struct bptbl_s {
    int refcount;        /**< Reference count. */
    char *name;          /**< Name for diagnostic output. */
    dict2pid_t *d2p;     /**< Tied state mapping. */

    ptmr_t t_bptbl;      /**< Performance timer. */

    garray_t *retired;   /**< Retired backpointer entries. */
    garray_t *ent;       /**< Active backpointer entries. */
    garray_t *rc;        /**< Right context scores for word exits. */

    int32 n_frame;       /**< Number of frames searched. */
    /**
     * First frame containing active backpointers.  Also the first
     * frame which is indexed by end frame.  ef_idx indices should be
     * added to this to get actual frame indices.
     */
    int32 active_fr;
    /**
     * Oldest retired backpointer referenced by active backpointers.
     * This bp's endpoint plus one is the first start frame of active
     * backpointers.  Although no new backpointers will be generated
     * with starting frames before active_fr, and thus all
     * backpointers before active_fr can be retired, existing active
     * backpointers will still have starting frames before active_fr.
     * However, no active backpointer has a starting frame before this.
     */
    int32 oldest_bp;
    int32 dest_s_idx;        /**< bscorestack index corresponding to
                              * first_invert_bp (which is invalid by
                              * definition) */

    garray_t *permute;      /**< Current permutation of entries (indexed
                            * by ent - first_invert_bp). */
    garray_t *ef_idx;    /**< First BPTable entry exiting in each frame */

    /* This is indexed by frame - active_fr */
    int32 n_frame_alloc; /**< Number of frames allocated for frame-based arrays. */
    bitvec_t *valid_fr;  /**< Set of accessible frames (used in gc) */
} bptbl_t;


/**
 * Data type used to score right-context score deltas.
 */
typedef uint16 rcdelta_t;

/**
 * Indicates no such right-context exists.
 */
#define NO_RC ((rcdelta_t)-1)

/**
 * Create a new bptbl.
 */
bptbl_t *bptbl_init(char const *name,
                    dict2pid_t *d2p, int n_alloc, int n_frame_alloc);

/**
 * Retain a pointer to a bptbl.
 */
bptbl_t *bptbl_retain(bptbl_t *bpt);

/**
 * Release a bptbl.
 */
int bptbl_free(bptbl_t *bpt);

/**
 * Dump contents of a bptbl for debugging.
 */
void bptbl_dump(bptbl_t *bptbl);

/**
 * Clear the backpointer table.
 */
void bptbl_reset(bptbl_t *bptbl);

/**
 * Check if the backpointer table is final.
 */
int bptbl_is_final(bptbl_t *bptbl);

/**
 * Finalize the backpointer table.
 *
 * Garbage collects and retires all active backpointers.
 */
int bptbl_finalize(bptbl_t *bptbl);

/**
 * Find the best early exit from an utterance, exp
 */
bpidx_t bptbl_find_exit(bptbl_t *bptbl, int32 wid);

/**
 * Record the current frame's index in the backpointer table.
 *
 * @param oldest_bp Index of the oldest backpointer still active in
 *                  the search graph.
 * @return the new frame index
 */
int bptbl_push_frame(bptbl_t *bptbl, bpidx_t oldest_bp);

/**
 * Release retired backpointers before a given index.
 */
int bptbl_release(bptbl_t *bptbl, bpidx_t first_idx);

/**
 * Add a backpointer to the table.
 */
bpidx_t bptbl_enter(bptbl_t *bptbl, int32 w, int32 path,
                    int32 score, int rc);

/**
 * Commit all valid backpointers from the current frame.
 */
int bptbl_commit(bptbl_t *bptbl);

/**
 * Get the index of the current frame from the backpointer table.
 */
int bptbl_frame_idx(bptbl_t *bptbl);

/**
 * Obtain the index of the first word exit for a given frame.
 */
bpidx_t bptbl_ef_idx(bptbl_t *bptbl, int frame_idx);

/**
 * Obtain the first active backpointer index.
 */
bpidx_t bptbl_active_idx(bptbl_t *bptbl);

/**
 * Obtain the first active frame's index.
 */
int bptbl_active_frame(bptbl_t *bptbl);

/**
 * Obtain the index of the end of the backpointer table (the last index plus one).
 */
bpidx_t bptbl_end_idx(bptbl_t *bptbl);

/**
 * Obtain the index of the end of the retired portion of the table (last index plus one).
 */
bpidx_t bptbl_retired_idx(bptbl_t *bptbl);

/**
 * Obtain the index of the first active start frame.
 *
 * Although no new backpointers will be generated with starting frames
 * before bptbl_active_fr(), and thus all backpointers before
 * bptbl_active_fr() are retired, existing active backpointers will
 * still have starting frames before bptbl_active_fr().  However, no
 * active backpointer has a starting frame before this one.
 */
int bptbl_active_sf(bptbl_t *bptbl);

/**
 * Get a backpointer entry.
 */
int bptbl_get_bp(bptbl_t *bptbl, bpidx_t bpidx, bp_t *out_ent);

/**
 * Get pointer to a backpointer entry.
 */
bp_t *bptbl_ent(bptbl_t *bptbl, bpidx_t bpidx);

/**
 * Update a backpointer entry.
 */
int bptbl_set_bp(bptbl_t *bptbl, bpidx_t bpidx, bp_t const *ent);

/**
 * Get the start frame for a given backpointer index.
 */
int bptbl_sf(bptbl_t *bptbl, bpidx_t bpidx);

/**
 * Get the number of word exits in a given frame.
 */
int bptbl_ef_count(bptbl_t *bptbl, int frame_idx);

/**
 * Set a right context score for a given backpointer entry.
 */
void bptbl_set_rcscore(bptbl_t *bptbl, bpidx_t bpidx, int rc, int32 score);

/**
 * Get the array of right context scores for a backpointer entry.
 *
 * @return Number of right context scores.
 */
int bptbl_get_rcscores(bptbl_t *bptbl, bpidx_t bpidx, int32 *out_rcscores);

/**
 * Get the array of right context deltasfor a backpointer entry.
 *
 * @return Number of right context scores.
 */
int bptbl_get_rcdeltas(bptbl_t *bptbl, bpidx_t bpidx, rcdelta_t *out_rcdeltas);

/**
 * Update best score and language model state for a backpointer table entry.
 */
void bptbl_update_bp(bptbl_t *bptbl, int32 bp, int rc,
                     bpidx_t new_path, int32 new_score);

/**
 * Calculate approximate language model score for an entry.
 */
int32 bptbl_fake_lmscore(bptbl_t *bptbl, ngram_model_t *lm,
                         bpidx_t bp, int *out_n_used);

/**
 * Construct a hypothesis string by backtracing from an entry.
 * @param bp Entry from which to backtrace
 * @return Newly allocated hypothesis string (free with ckd_free())
 */
char *bptbl_backtrace(bptbl_t *bptbl, bpidx_t bp);

/**
 * Construct a hypothesis string from the best path in bptbl.
 *
 * @param out_score Output: Score of hypothesis.
 * @param wid End word ID (or BAD_S3WID to return the best exit
 * regardless of word ID)
 * @return Newly allocated hypothesis string (free with ckd_free()).
 */
char *bptbl_hyp(bptbl_t *bptbl, int32 *out_score, int32 finish_wid);

/**
 * Construct a segmentation iterator from the best path in bptbl.
 *
 * @param out_score Output: Score of hypothesis.
 * @param wid End word ID (or BAD_S3WID to return the best exit
 * regardless of word ID)
 * @return New segmentation iterator.
 */
struct seg_iter_s *bptbl_seg_iter(bptbl_t *bptbl, int32 *out_score, int32 finish_wid);

/**
 * Construct a segmentation iterator from the best path ending with a given bp.
 *
 * @param out_score Output: Score of hypothesis.
 * @param bp Entry from which to backtrace.
  * @return New segmentation iterator.
 */
struct seg_iter_s *bptbl_seg_backtrace(bptbl_t *bptbl, bpidx_t bp);

#endif /* __BPTBL_H__ */
