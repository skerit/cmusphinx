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
 * @file search_internal.h Search algorithm class (internal API)
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#ifndef __PS_SEARCH_INTERNAL_H__
#define __PS_SEARCH_INTERNAL_H__

#include <sphinxbase/logmath.h>
#include <sphinxbase/sbthread.h>
#include <sphinxbase/profile.h>
#include <sphinxbase/cmd_ln.h>

#include <multisphinx/search.h>
#include <multisphinx/acmod.h>
#include <multisphinx/bptbl.h>
#include <multisphinx/arc_buffer.h>
#include <multisphinx/dict.h>
#include <multisphinx/dict2pid.h>

/**
 * V-table for search algorithm functions, not called directly by users.
 */
struct searchfuncs_s {
    char const *name;

    search_t *(*init)(search_t *other, cmd_ln_t *config, acmod_t *acmod, dict2pid_t *d2p);
    int (*free)(search_t *search);   /**< Free search-specific stuff. */
    int (*decode)(search_t *search); /**< Decode an utterance. */

    char const *(*hyp)(search_t *search, int32 *out_score);
    int32 (*prob)(search_t *search);
    seg_iter_t *(*seg_iter)(search_t *search, int32 *out_score);

    bptbl_t *(*bptbl)(search_t *search);
    ngram_model_t *(*lmset)(search_t *search);
};

/**
 * Base structure for search module.
 */
struct search_s {
    searchfuncs_t *vt;  /**< V-table of search methods. */
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
    char *uttid;

    struct arc_buffer_s *input_arcs;  
    struct arc_buffer_s *output_arcs;

    search_cb_func cb;
    void *cb_data;

    /* Magical word IDs that must exist in the dictionary: */
    int32 start_wid;       /**< Start word ID. */
    int32 silence_wid;     /**< Silence word ID. */
    int32 finish_wid;      /**< Finish word ID. */
};

/* A variety of accessors. */
#define search_base(s) ((search_t *)s)
#define search_thread(s) search_base(s)->thr
#define search_acmod(s) search_base(s)->acmod
#define search_dict(s) search_base(s)->dict
#define search_dict2pid(s) search_base(s)->d2p
#define search_post(s) search_base(s)->post
#define search_input_arcs(s) search_base(s)->input_arcs
#define search_output_arcs(s) search_base(s)->output_arcs
#define search_n_words(s) search_base(s)->n_words
#define search_silence_wid(s) search_base(s)->silence_wid
#define search_start_wid(s) search_base(s)->start_wid
#define search_finish_wid(s) search_base(s)->finish_wid
#define search_cb(s) search_base(s)->cb;
#define search_cb_data(s) search_base(s)->cb_data;

/**
 * Initialize base structure.
 */
void search_base_init(search_t *search, searchfuncs_t *vt,
                      cmd_ln_t *config, acmod_t *acmod,
                      dict2pid_t *d2p);

/**
 * De-initialize base structure.
 */
void search_deinit(search_t *search);

/**
 * V-table for segmentation iterators.
 */
typedef struct segfuncs_s {
    seg_iter_t *(*seg_next)(seg_iter_t *seg);
    void (*seg_free)(seg_iter_t *seg);
} segfuncs_t;

/**
 * Base structure for hypothesis segmentation iterator.
 */
struct seg_iter_s {
    segfuncs_t *vt;     /**< V-table of seg methods */
    search_t *search;   /**< Search object from whence this came */
    int32 wid;          /**< Word ID. */
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

#define search_seg_next(seg) (*(seg->vt->seg_next))(seg)
#define search_seg_free(s) (*(seg->vt->seg_free))(seg)

/**
 * Call an event simply.
 */
int search_call_event(search_t *search, int event, int frame);

#endif /* __PS_SEARCH_INTERNAL_H__ */
