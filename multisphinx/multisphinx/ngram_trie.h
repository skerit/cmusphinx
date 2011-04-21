/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2010 Carnegie Mellon University.  All rights
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
 * @file ngram_trie.h
 * @brief Mutable trie implementation of N-Gram language models
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#ifndef __NGRAM_TRIE_H__
#define __NGRAM_TRIE_H__

#include <stdarg.h>
#include <stdio.h>

#include <sphinxbase/logmath.h>

#include "dict.h"

/**
 * Trie-based N-Gram model structure.
 */
typedef struct ngram_trie_s ngram_trie_t;

/**
 * Node (individual N-Gram) in the trie.
 */
typedef struct ngram_trie_node_s ngram_trie_node_t;

/**
 * Iterator over N-Grams in a trie.
 */
typedef struct ngram_trie_iter_s ngram_trie_iter_t;

/**
 * Initialize an empty N-Gram trie.
 */
ngram_trie_t *ngram_trie_init(dict_t *d, logmath_t *lmath);

ngram_trie_t *ngram_trie_retain(ngram_trie_t *t);
int ngram_trie_free(ngram_trie_t *t);
dict_t *ngram_trie_dict(ngram_trie_t *t);
logmath_t *ngram_trie_logmath(ngram_trie_t *t);
int32 ngram_trie_zero(ngram_trie_t *t);

/**
 * Get the order of the N-Gram model.
 */
int ngram_trie_n(ngram_trie_t *t);

/**
 * Read N-Grams from an ARPA text format file.
 */
int ngram_trie_read_arpa(ngram_trie_t *t, FILE *arpafile);

/**
 * Write N-Grams to an ARPA text format file.
 */
int ngram_trie_write_arpa(ngram_trie_t *t, FILE *arpafile);

/**
 * Get the root node of the trie.
 */
ngram_trie_node_t *ngram_trie_root(ngram_trie_t *t);

/**
 * Look up an N-Gram in the trie by strings.
 */
ngram_trie_node_t *ngram_trie_ngram(ngram_trie_t *t, char const *w, ...); 

/**
 * Look up an N-Gram in the trie by numeric IDs.
 */
ngram_trie_node_t *ngram_trie_ngram_v(ngram_trie_t *t, int32 w,
				      int32 const *hist, int32 n_hist);

/**
 * Create an N-Gram (and its parents) in the trie by strings.
 */
ngram_trie_node_t *ngram_trie_ngram_init(ngram_trie_t *t, char const *w, ...); 

/**
 * Create an N-Gram (and its parents) in the trie by strings.
 */
ngram_trie_node_t *ngram_trie_ngram_init_v(ngram_trie_t *t, int32 w,
                                           int32 const *hist, int32 n_hist);

/**
 * Get model probability (with backoff) for a word with history.
 */
int32 ngram_trie_prob(ngram_trie_t *t, int *n_used, char const *w, ...);

/**
 * Get model probability (with backoff) for a word with history.
 */
int32 ngram_trie_prob_v(ngram_trie_t *t, int *n_used, int32 w,
			int32 const *hist, int32 n_hist);

/**
 * Get an iterator over all N-Grams of a given order in the trie.
 */
ngram_trie_iter_t *ngram_trie_ngrams(ngram_trie_t *t, int n);

/**
 * Get an iterator over all successors to an N-Gram.
 */
ngram_trie_iter_t *ngram_trie_successors(ngram_trie_t *t, ngram_trie_node_t *h);

/**
 * Get the raw array of successors.
 */
ngram_trie_node_t **ngram_trie_successors_unchecked(ngram_trie_t *t,
                                                    ngram_trie_node_t *h, size_t *out_nsucc);

/**
 * Move iterator forward.
 */
ngram_trie_iter_t *ngram_trie_iter_next(ngram_trie_iter_t *itor);

/**
 * Move iterator up one level.
 */
ngram_trie_iter_t *ngram_trie_iter_up(ngram_trie_iter_t *itor);

/**
 * Move iterator down one level.
 */
ngram_trie_iter_t *ngram_trie_iter_down(ngram_trie_iter_t *itor);

/**
 * Get the node pointed to by an iterator.
 */
ngram_trie_node_t *ngram_trie_iter_get(ngram_trie_iter_t *itor);

/**
 * Get the parent of the node pointed to by an iterator.
 */
ngram_trie_node_t *ngram_trie_iter_get_parent(ngram_trie_iter_t *itor);

/**
 * Get the word ID from a node.
 */
int32 ngram_trie_node_word(ngram_trie_t *t, ngram_trie_node_t *node);

/**
 * Set the word ID in a node.
 */
void ngram_trie_node_set_word(ngram_trie_t *t, ngram_trie_node_t *node,
                              int32 wid);

/**
 * Get the probability and backoff weight values from a node.
 */
void ngram_trie_node_params(ngram_trie_t *t,
                            ngram_trie_node_t *node,
                            int32 *out_log_prob,
                            int32 *out_log_bowt);

/**
 * Set the probability and backoff weight values in a node.
 */
void ngram_trie_node_set_params(ngram_trie_t *t,
                                ngram_trie_node_t *node,
                                int32 log_prob,
                                int32 log_bowt);

/**
 * Get the "raw" (unscaled) probability and backoff weight values from
 * a node.
 */
void ngram_trie_node_params_raw(ngram_trie_t *t,
                                ngram_trie_node_t *node,
                                int16 *out_log_prob,
                                int16 *out_log_bowt);

/**
 * Set the "raw" (unscaled) probability and backoff weight values in a node.
 */
void ngram_trie_node_set_params_raw(ngram_trie_t *t,
                                    ngram_trie_node_t *node,
                                    int16 log_prob,
                                    int16 log_bowt);

/**
 * Look up a successor to an N-Gram.
 */
ngram_trie_node_t *ngram_trie_successor(ngram_trie_t *t, ngram_trie_node_t *h, int32 w);

/**
 * Delete a successor to an N-Gram
 */
int ngram_trie_delete_successor(ngram_trie_t *t, ngram_trie_node_t *h, int32 w);

/**
 * Add a new successor to an N-Gram
 */
ngram_trie_node_t *ngram_trie_add_successor(ngram_trie_t *t,
                                            ngram_trie_node_t *h, int32 w);

/**
 * Create a new node with no history or successors.
 */
ngram_trie_node_t *ngram_trie_node_alloc(ngram_trie_t *t);

/**
 * Add a new successor to an N-Gram
 */
int ngram_trie_add_successor_ngram(ngram_trie_t *t, 
                                   ngram_trie_node_t *h,
				   ngram_trie_node_t *w);

/**
 * Change the head word of a successor N-Gram.
 */
int ngram_trie_rename_successor(ngram_trie_t *t,
                                ngram_trie_node_t *h,
                                ngram_trie_node_t *w,
                                int32 new_wid);

/**
 * Get the backoff N-gram (one fewer element of history) for an N-Gram.
 */
ngram_trie_node_t *ngram_trie_backoff(ngram_trie_t *t,
				      ngram_trie_node_t *ng);

/**
 * Get the model probability (with backoff) for a successor.
 */
int32 ngram_trie_successor_prob(ngram_trie_t *t,
				ngram_trie_node_t *h, int32 w);

/**
 * Calculate the backoff weight for an N-Gram.
 */
int32 ngram_trie_calc_bowt(ngram_trie_t *t, ngram_trie_node_t *h);

/**
 * Validate successor distribution for an N-Gram.
 */
int32 ngram_trie_node_validate(ngram_trie_t *t, ngram_trie_node_t *h);

/**
 * Get the order of a node.
 */
int ngram_trie_node_n(ngram_trie_t *t, ngram_trie_node_t *ng);

/**
 * Get the word history from a node.
 */
int32 ngram_trie_node_get_word_hist(ngram_trie_t *t,
                                    ngram_trie_node_t *ng,
                                    int32 *out_hist);

/**
 * Print out the words from an N-gram.
 */
int ngram_trie_node_print(ngram_trie_t *t,
                          ngram_trie_node_t *ng,
                          FILE *outfh);

/**
 * Update the N-Gram order and counts (the add and delete functions do
 * not do this because it is much too costly).
 */
int ngram_trie_update_counts(ngram_trie_t *t);

#endif /* __NGRAM_TRIE_H__ */
