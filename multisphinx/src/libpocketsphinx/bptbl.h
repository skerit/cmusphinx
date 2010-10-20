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

/* Local headers. */
#include "pocketsphinx_internal.h"

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
    dict2pid_t *d2p;     /**< Tied state mapping. */

    bp_t *retired;       /**< Retired backpointer entries. */
    int32 n_retired;     /**< Number of retired entries. */
    int32 n_retired_alloc; /**< Allocated size of retired array. */

    bp_t *ent;           /**< Active backpointer entries. */
    int32 n_ent;         /**< First free BPTable entry. */
    int32 n_ent_alloc;   /**< Allocated size of active array. */

    int32 *bscore_stack;     /**< Array containing right context scores for word exits. */
    int32 bss_head;          /**< First free BScoreStack entry */
    int32 bscore_stack_size; /**< Number of entries allocated in bscore_stack. */

    int32 n_frame;       /**< Number of frames searched. */
    /**
     * First frame containing active backpointers.  Also the first
     * frame which is indexed by end frame.  ef_idx indices should be
     * added to this to get actual frame indices.
     */
    int32 active_fr;
    /**
     * Index of first bptbl entry to be reordered by compaction and
     * insertion sort.  Also marks the end of the region of retired
     * bptbl entries.
     */
    bpidx_t first_invert_bp;
    int32 dest_s_idx;        /**< bscorestack index corresponding to
                              * first_invert_bp (which is invalid by
                              * definition) */
    bpidx_t *permute;      /**< Current permutation of entries (indexed
                            * by ent - first_invert_bp). */
    int32 n_permute_alloc;

    /* All these are indexed by frame - active_fr */
    int32 n_frame_alloc; /**< Number of frames allocated for frame-based arrays. */
    bpidx_t *ef_idx;       /**< First BPTable entry exiting in each frame */
    bitvec_t *valid_fr;  /**< Set of accessible frames (used in gc) */
} bptbl_t;


/**
 * Segmentation "iterator" for backpointer table results.
 */
typedef struct bptbl_seg_s {
    ps_seg_t base;  /**< Base structure. */
    bptbl_t *bptbl; /**< Backpointer table. */
    int32 *bpidx;   /**< Sequence of backpointer IDs. */
    int16 n_bpidx;  /**< Number of backpointer IDs. */
    int16 cur;      /**< Current position in bpidx. */
} bptbl_seg_t;


/**
 * Create a new bptbl.
 */
bptbl_t *bptbl_init(dict2pid_t *d2p, int n_alloc, int n_frame_alloc);

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
 * Finalize the backpointer table.
 *
 * Garbage collects and retires all active backpointers.
 */
int bptbl_finalize(bptbl_t *bptbl);

/**
 * Find the best scoring exit in the final frame.
 *
 * @param wid End word ID (or BAD_S3WID to return the best exit
 * regardless of word ID)
 */
bp_t *bptbl_find_exit(bptbl_t *bptbl, int32 wid);

/**
 * Record the current frame's index in the backpointer table.
 *
 * @param oldest_bp Index of the oldest backpointer still active in
 *                  the search graph.
 * @return the new frame index
 */
int bptbl_push_frame(bptbl_t *bptbl, bpidx_t oldest_bp);

/**
 * Add a backpointer to the table.
 */
bp_t *bptbl_enter(bptbl_t *bptbl, int32 w, int32 path,
                  int32 score, int rc);

/**
 * Get the index of the current frame from the backpointer table.
 */
int bptbl_frame_idx(bptbl_t *bptbl);

/**
 * Obtain the index of the first word exit for a given frame.
 */
bpidx_t bptbl_ef_idx(bptbl_t *bptbl, int frame_idx);

/**
 * Obtain a pointer to the backpointer with a given index.
 */
bp_t *bptbl_ent(bptbl_t *bptbl, bpidx_t bpidx);

/**
 * Obtain the best predecessor to a given backpointer entry.
 */
bp_t *bptbl_prev(bptbl_t *bptbl, bp_t *ent);

/**
 * Get the start frame for a given backpointer index.
 */
int bptbl_sf(bptbl_t *bptbl, bpidx_t bpidx);

/**
 * Get the index for a given backpointer.
 */
bpidx_t bptbl_idx(bptbl_t *bptbl, bp_t *bpe);

/**
 * Get the number of word exits in a given frame.
 */
int bptbl_ef_count(bptbl_t *bptbl, int frame_idx);

/**
 * Update language model state for a backpointer table entry.
 *
 * This is called "fake" language model state since backpointers don't
 * actually have unique LM states.  In other words this is "poor man's
 * trigram" search.
 */
void bptbl_fake_lmstate(bptbl_t *bptbl, int32 bp);

/**
 * Construct a hypothesis string from the best path in bptbl.
 *
 * @param out_score Output: Score of hypothesis.
 * @param wid End word ID (or BAD_S3WID to return the best exit
 * regardless of word ID)
 * @return Newly allocated hypothesis string (free with ckd_free()).
 */
char *
bptbl_hyp(bptbl_t *bptbl, int32 *out_score, int32 finish_wid);

/**
 * Construct a segmentation iterator from the best path in bptbl.
 *
 * @param out_score Output: Score of hypothesis.
 * @param wid End word ID (or BAD_S3WID to return the best exit
 * regardless of word ID)
 * @return New segmentation iterator.
 */
ps_seg_t *bptbl_seg_iter(bptbl_t *bptbl, int32 *out_score, int32 finish_wid);

#endif /* __BPTBL_H__ */
