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

#define bp_sf(bpt,i) (((bpt)->ent[i].bp == NO_BP) ? 0 \
                      : (bpt)->ent[(bpt)->ent[i].bp].frame + 1)

/**
 * Back pointer table
 */
typedef struct bptbl_s {
    bp_t *ent;       /* Forward pass lattice */
    dict_t *dict;    /* Dictionary */
    int32 n_ent;             /* First free BPTable entry */
    int32 n_alloc;
    int32 first_unsorted; /**< First entry that has not yet been sorted/gc-ed */
    int32 *bscore_stack;     /* Score stack for all possible right contexts */
    int32 bss_head;          /* First free BScoreStack entry */
    int32 bscore_stack_size;
    int32 n_frame_alloc; /**< Number of frames allocated in bp_table_idx and friends. */
    int32 n_frame;       /**< Number of frames actually present. */
    int32 *frame_idx;    /**< First BPTable entry for each frame */
    int32 *word_idx;     /**< BPTable index for any word in current frame;
                            cleared before each frame */
    ps_latnode_t **frm_wordlist;   /**< List of active words in each frame. */
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


bptbl_t *bptbl_init(dict_t *dict, int n_alloc, int n_frame_alloc);

void bptbl_free(bptbl_t *bpt);

void dump_bptable(bptbl_t *bptbl);

/**
 * Record the current frame's index in the backpointer table.
 *
 * @return the current backpointer index.
 */
int bptbl_push_frame(bptbl_t *bptbl, int frame_idx);

/**
 * Order backpointer table entries according to start frame and remove
 * invalid paths.
 */
void bptable_gc(bptbl_t *ngs, int oldest_bp, int frame_idx);


#endif /* __BPTBL_H__ */
