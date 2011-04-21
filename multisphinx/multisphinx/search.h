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
 * @file search.h Search algorithm class
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#ifndef __PS_SEARCH_H__
#define __PS_SEARCH_H__

#include <sphinxbase/sbthread.h>

#include <multisphinx/arc_buffer.h>
#include <multisphinx/bptbl.h>

/**
 * Search algorithm structure.
 */
typedef struct search_s search_t;

/**
 * Structure describing a search algorithm.
 */
typedef struct searchfuncs_s searchfuncs_t;

/**
 * Segmentation iterator.
 */
typedef struct seg_iter_s seg_iter_t;

/**
 * Event in search.
 */
typedef struct search_event_s search_event_t;
struct search_event_s {
    int16 event;
    int16 frame;
};

/**
 * Event types.
 */
enum search_event_e {
    SEARCH_START_UTT,
    SEARCH_PARTIAL_RESULT,
    SEARCH_FINAL_RESULT,
    SEARCH_END_UTT
};

/**
 * Callback for search events.
 */
typedef int (*search_cb_func)(search_t *search, search_event_t *evt, void *udata);

/**
 * Start a search thread.
 */
sbthread_t *search_run(search_t *search);

/**
 * Wait for a search thread to complete.
 */
int search_wait(search_t *search);

/**
 * Free a search structure.
 */
int search_free(search_t *search);

/**
 * Get the name of the search implementation.
 */
char const *search_name(search_t *search);

/**
 * Link one search structure to another via an arc buffer.
 */
arc_buffer_t *search_link(search_t *from, search_t *to,
                          char const *name, int keep_scores);

/**
 * Add a callback for search events.
 */
void search_set_cb(search_t *search, search_cb_func cb, void *udata);

/**
 * Get the latest hypothesis from a search.
 */
char const *search_hyp(search_t *search, int32 *out_score);

/**
 * Get the latest segmentation from a search
 */
seg_iter_t *search_seg_iter(search_t *search, int32 *out_score);

/**
 * Move iterator vorwart.
 */
seg_iter_t *seg_iter_next(seg_iter_t *itor);

/**
 * Free iterator.
 */
void seg_iter_free(seg_iter_t *itor);

/**
 * Get word from iterator.
 */
char const *seg_iter_word(seg_iter_t *itor);

/**
 * Get word ID from iterator.
 */
int32 seg_iter_wid(seg_iter_t *itor);

/**
 * Get frames from iterator.
 */
void seg_iter_times(seg_iter_t *itor, int *out_sf, int *out_ef);

/**
 * Get the backpointer table, if any, from a search.
 */
bptbl_t *search_bptbl(search_t *search);

/**
 * Get the N-Gram language model set, if any, from a search.
 */
ngram_model_t *search_lmset(search_t *search);

/**
 * Get the configuration from a search.
 */
cmd_ln_t *search_config(search_t *search);

/**
 * Get the current utterance ID from a search.
 */
char const *search_uttid(search_t *search);

#endif /* __PS_SEARCH_H__ */
