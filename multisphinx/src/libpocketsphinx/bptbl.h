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
 * Back pointer table entry
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

#define bp_sf(bpt,i) ((i == NO_BP || (bpt)->ent[i].bp == NO_BP) ? 0 \
                      : (bpt)->ent[(bpt)->ent[i].bp].frame + 1)

/**
 * Back pointer table
 */
typedef struct bptbl_s {
    dict2pid_t *d2p;     /**< Tied state mapping. */
    bp_t *ent;           /**< Backpointer entries. */
    int32 n_ent;         /**< First free BPTable entry. */
    int32 n_alloc;       /**< Number of entries allocated for entry-based arrays (ent, permute) */
    int32 active_sf;     /**< First frame containing active backpointers. */
    int32 *permute;      /**< Current permutation of entries (used for gc/sorting). */
    int32 first_invert_bp; /**< First reordered backpointer (used in gc) */

    int32 *bscore_stack;     /**< Array containing right context scores for word exits. */
    int32 bss_head;          /**< First free BScoreStack entry */
    int32 bscore_stack_size; /**< Number of entries allocated in bscore_stack. */

    int32 n_frame_alloc; /**< Number of frames allocated for frame-based arrays. */
    int32 n_frame;       /**< Number of frames actually present. */
    int32 *ef_idx;       /**< First BPTable entry exiting in each frame */
    int32 *sf_idx;       /**< First BPTable entry entering in each frame */
    bitvec_t *valid_fr;  /**< Set of accessible frames (used in gc) */
    ps_latnode_t **frm_wordlist;   /**< List of active words in each frame. */

    int32 *word_idx;     /**< BPTable index for any word in current frame;
                            cleared before each frame */
} bptbl_t;


/**
 * Segmentation "iterator" for backpointer table results.
 */
typedef struct bptbl_seg_s {
    ps_seg_t base;  /**< Base structure. */
    int32 *bpidx;   /**< Sequence of backpointer IDs. */
    int16 n_bpidx;  /**< Number of backpointer IDs. */
    int16 cur;      /**< Current position in bpidx. */
} bptbl_seg_t;


bptbl_t *bptbl_init(dict2pid_t *d2p, int n_alloc, int n_frame_alloc);

void bptbl_free(bptbl_t *bpt);

void dump_bptable(bptbl_t *bptbl, int start, int end);

/**
 * Record the current frame's index in the backpointer table.
 *
 * @return the current backpointer index.
 */
int bptbl_push_frame(bptbl_t *bptbl, int oldest_bp, int frame_idx);

/**
 * Add a backpointer to the table.
 */
bp_t *bptbl_enter(bptbl_t *bptbl, int32 w, int frame_idx,
                  int32 path, int32 score, int rc);

/**
 * Clear the backpointer table.
 */
void bptbl_reset(bptbl_t *bptbl);

/**
 * Cache trigram predecessors for a backpointer table entry.
 */
void bptbl_fake_lmstate(bptbl_t *bptbl, int32 bp);

#endif /* __BPTBL_H__ */
