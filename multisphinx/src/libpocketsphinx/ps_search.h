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
 * @file ps_search.h Search algorithm class
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#ifndef __PS_SEARCH_H__
#define __PS_SEARCH_H__

/* SphinxBase headers. */
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/logmath.h>
#include <sphinxbase/sbthread.h>
#include <sphinxbase/profile.h>

/* Local headers. */
#include "pocketsphinx.h"
#include "acmod.h"
#include "dict.h"
#include "dict2pid.h"

/* Can't actually do this because of header loops or some such evil */
/* #include "arc_buffer.h" */

/**
 * Search algorithm structure.
 */
typedef struct ps_search_s ps_search_t;

/**
 * V-table for search algorithm functions, not called directly by users.
 */
typedef struct ps_searchfuncs_s {
    char const *name;

    int (*free)(ps_search_t *search);   /**< Free search-specific stuff. */
    int (*decode)(ps_search_t *search); /**< Decode an utterance. */

    char const *(*hyp)(ps_search_t *search, int32 *out_score);
    int32 (*prob)(ps_search_t *search);
    ps_seg_t *(*seg_iter)(ps_search_t *search, int32 *out_score);
} ps_searchfuncs_t;

struct arc_buffer_s;

/**
 * Base structure for search module.
 */
struct ps_search_s {
    ps_searchfuncs_t *vt;  /**< V-table of search methods. */
    sbthread_t *thr;       /**< Thread in which this search runs. */
    sbmtx_t *mtx;          /**< Lock for this search. */
    ptmr_t t;              /**< Overall performance timer for this search. */
    int32 total_frames;    /**< Total number of frames processed. */

    cmd_ln_t *config;      /**< Configuration. */
    acmod_t *acmod;        /**< Acoustic model. */
    dict_t *dict;          /**< Pronunciation dictionary. */
    dict2pid_t *d2p;       /**< Dictionary to senone mappings. */
    char *hyp_str;         /**< Current hypothesis string. */
    int32 post;            /**< Utterance posterior probability. */
    int32 n_words;         /**< Number of words known to search (may
                              be less than in the dictionary) */

    struct arc_buffer_s *output_arcs; /**< Arc buffer, used to forward search
                                       * results to next pass of search. */

    /* Magical word IDs that must exist in the dictionary: */
    int32 start_wid;       /**< Start word ID. */
    int32 silence_wid;     /**< Silence word ID. */
    int32 finish_wid;      /**< Finish word ID. */
};

/* A variety of accessors. */
#define ps_search_base(s) ((ps_search_t *)s)
#define ps_search_thread(s) ps_search_base(s)->thr
#define ps_search_config(s) ps_search_base(s)->config
#define ps_search_acmod(s) ps_search_base(s)->acmod
#define ps_search_dict(s) ps_search_base(s)->dict
#define ps_search_dict2pid(s) ps_search_base(s)->d2p
#define ps_search_post(s) ps_search_base(s)->post
#define ps_search_output_arcs(s) ps_search_base(s)->output_arcs
#define ps_search_n_words(s) ps_search_base(s)->n_words
#define ps_search_silence_wid(s) ps_search_base(s)->silence_wid
#define ps_search_start_wid(s) ps_search_base(s)->start_wid
#define ps_search_finish_wid(s) ps_search_base(s)->finish_wid

/**
 * Initialize base structure.
 */
void ps_search_init(ps_search_t *search, ps_searchfuncs_t *vt,
                    cmd_ln_t *config, acmod_t *acmod, dict_t *dict,
                    dict2pid_t *d2p);

/**
 * De-initialize base structure.
 */
void ps_search_deinit(ps_search_t *search);

/**
 * Start a search thread.
 */
sbthread_t *ps_search_run(ps_search_t *search);

/**
 * Wait for a search thread to complete.
 */
int ps_search_wait(ps_search_t *search);

/**
 * Free a search structure.
 */
int ps_search_free(ps_search_t *search);

/**
 * Get the latest hypothesis from a search.
 *
 * FIXME: This will probably go away due to hypothesis splicing
 */
char const *ps_search_hyp(ps_search_t *search, int32 *out_score);

/**
 * V-table for segmentation iterators.
 */
typedef struct ps_segfuncs_s {
    ps_seg_t *(*seg_next)(ps_seg_t *seg);
    void (*seg_free)(ps_seg_t *seg);
} ps_segfuncs_t;

/**
 * Base structure for hypothesis segmentation iterator.
 */
struct ps_seg_s {
    ps_segfuncs_t *vt;     /**< V-table of seg methods */
    ps_search_t *search;   /**< Search object from whence this came */
    char const *word;      /**< Word string (pointer into dictionary hash) */
    int16 sf;                /**< Start frame. */
    int16 ef;                /**< End frame. */
    int32 ascr;            /**< Acoustic score. */
    int32 lscr;            /**< Language model score. */
    int32 prob;            /**< Log posterior probability. */
    /* This doesn't need to be 32 bits, so once the scores above are
     * reduced to 16 bits (or less!), this will be too. */
    int32 lback;           /**< Language model backoff. */
    /* Not sure if this should be here at all. */
    float32 lwf;           /**< Language weight factor (for second-pass searches) */
};

#define ps_search_seg_next(seg) (*(seg->vt->seg_next))(seg)
#define ps_search_seg_free(s) (*(seg->vt->seg_free))(seg)

/**
 * Get the latest segmentation from a search
 *
 * FIXME: This will probably go away due to hypothesis splicing
 */
ps_seg_t *ps_search_seg_iter(ps_search_t *search, int32 *out_score);

#endif /* __PS_SEARCH_H__ */
