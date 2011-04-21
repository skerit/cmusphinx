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
 * Link one search structure to another via an arc buffer.
 */
arc_buffer_t *search_link(search_t *from, search_t *to,
                          char const *name, int keep_scores);

/**
 * Get the latest hypothesis from a search.
 *
 * FIXME: This will probably go away due to hypothesis splicing
 */
char const *search_hyp(search_t *search, int32 *out_score);

/**
 * Splice hypotheses from multiple searches.
 */
char const *search_splice(search_t **searches, int nsearches,
                          int32 *out_score);

/**
 * Get the backpointer table, if any, from a search.
 */
bptbl_t *search_bptbl(search_t *search);

/**
 * Get the N-Gram language model set, if any, from a search.
 */
ngram_model_t *search_lmset(search_t *search);

#endif /* __PS_SEARCH_H__ */
