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
 * @file ms_lattice.h Word lattices for MultiSphinx.
 */

#ifndef __MS_LATTICE_H__
#define __MS_LATTICE_H__

#include <stdio.h>

#include <sphinxbase/prim_type.h>
#include <sphinxbase/garray.h>
#include <sphinxbase/logmath.h>
#include <sphinxbase/ngram_model.h>

#include <multisphinx/dict.h>

#include "nodeid_map.h"

/**
 * Word lattice.
 */
typedef struct ms_lattice_s ms_lattice_t;

/**
 * Lattice node iterator.
 */
typedef struct ms_latnode_iter_s ms_latnode_iter_t;

/**
 * Lattice node structure.
 */
typedef struct ms_latnode_s {
    nodeid_t id;       /**< Node ID (language model state/word ID, sf) */
    int16 fan;         /**< Fan-in count for traversal. */
    garray_t *exits;   /**< Link indices (FIXME: not memory efficient) */
    garray_t *entries; /**< Link indices (FIXME: not memory efficient) */
} ms_latnode_t;

/**
 * Lattice link structure.
 */
typedef struct ms_latlink_t {
    int32 wid;    /**< Word ID. */
    int32 src;    /**< Source node ID. */
    int32 dest;   /**< Destination node ID. */
    int32 ascr;   /**< Acoustic score. */
    int32 lscr;   /**< Language score. */
    int32 alpha;  /**< Forward log-probability. */
    int32 beta;   /**< Backward log-probability. */
} ms_latlink_t;

ms_lattice_t *ms_lattice_init(logmath_t *lmath, dict_t *dict);
ms_lattice_t *ms_lattice_retain(ms_lattice_t *l);
int ms_lattice_free(ms_lattice_t *l);
dict_t *ms_lattice_dict(ms_lattice_t *l);
logmath_t *ms_lattice_lmath(ms_lattice_t *l);

/**
 * Create a language model state.
 */
int32 ms_lattice_lmstate_init(ms_lattice_t *l, int32 w,
                              int32 const *hist, int32 n_hist);

/**
 * Look up a language model state index by word IDs.
 */
int32 ms_lattice_get_lmstate_idx(ms_lattice_t *l, int32 w,
                                 int32 const *hist, int32 n_hist);

/**
 * Rotate a head word into a language model state array.
 *
 * FIXME: Needs actually to take one language model state and generate
 * another, or something suitably encapsulated as that.
 *
 * @param wid Head word of previous language model state.
 * @param hist History component of previous langauge model state.
 * @param n_hist Size of history component.
 * @param max_n_hist Allocation size of @a hist.
 * @return Number of words in newly created language model history
 *         array.
 */
int rotate_lmstate(int32 wid, int32 *hist,
                   int n_hist, int max_n_hist);

/**
 * Look up a language model state's word IDs by index.
 *
 * @return number of history entries.
 */
int ms_lattice_get_lmstate_wids(ms_lattice_t *l, int32 idx,
                                int32 *out_w, int32 *out_hist);

/**
 * Create a node.
 */
ms_latnode_t *ms_lattice_node_init(ms_lattice_t *l, int sf, int32 lmstate);

/**
 * Get a node by index.
 */
ms_latnode_t *ms_lattice_get_node_idx(ms_lattice_t *l, int32 idx);

/**
 * Get a node by ID (start frame + lm state).
 */
ms_latnode_t *ms_lattice_get_node_id(ms_lattice_t *l, int sf, int32 lmstate);

/**
 * Update a node's ID.
 */
ms_latnode_t *ms_lattice_set_node_id(ms_lattice_t *l, ms_latnode_t *node,
                                     int sf, int32 lmstate);

/**
 * Get the index of a node.
 */
int32 ms_lattice_get_idx_node(ms_lattice_t *l, ms_latnode_t *node);

/**
 * Set the start node.
 */
int32 ms_lattice_set_start(ms_lattice_t *l, ms_latnode_t *node);

/**
 * Set the end node.
 */
int32 ms_lattice_set_end(ms_lattice_t *l, ms_latnode_t *node);

/**
 * Get the start node.
 */
ms_latnode_t *ms_lattice_get_start(ms_lattice_t *l);

/**
 * Get the end node.
 */
ms_latnode_t *ms_lattice_get_end(ms_lattice_t *l);

/**
 * Create a link.
 */
ms_latlink_t *ms_lattice_link(ms_lattice_t *l,
			      ms_latnode_t *src, ms_latnode_t *dest,
			      int32 wid, int32 ascr);
/**
 * Get a link by index.
 */
ms_latlink_t *ms_lattice_get_link_idx(ms_lattice_t *l, int32 idx);

/**
 * Get a link's index.
 */
int32 ms_lattice_get_idx_link(ms_lattice_t *l, ms_latlink_t *link);

/**
 * Read a lattice in HTK format.
 */
int ms_lattice_read_htk(ms_lattice_t *l, FILE *fh, int frate);

/**
 * Write a lattice in HTK format.
 */
int ms_lattice_write_htk(ms_lattice_t *l, FILE *fh, int frate);

/**
 * Write a lattice in DOT format.
 */
int ms_lattice_write_dot(ms_lattice_t *l, FILE *fh);

/**
 * Begin a topological traversal of lattice nodes.
 */
ms_latnode_iter_t *ms_lattice_traverse_topo(ms_lattice_t *l,
                                            ms_latnode_t *end);

/**
 * Begin a reverse topological traversal of lattice nodes.
 */
ms_latnode_iter_t *ms_lattice_reverse_topo(ms_lattice_t *l,
                                           ms_latnode_t *start);

/**
 * Traverse all lattice nodes with a given start frame.
 */
ms_latnode_iter_t *ms_lattice_traverse_frame(ms_lattice_t *l,
                                             int frame_idx);

/**
 * Move to the next node in traversal.
 */
ms_latnode_iter_t *ms_latnode_iter_next(ms_latnode_iter_t *itor);

/**
 * Get current node in traversal.
 */
ms_latnode_t *ms_latnode_iter_get(ms_latnode_iter_t *itor);

/**
 * Get index of current node in traversal.
 */
int32 ms_latnode_iter_get_idx(ms_latnode_iter_t *itor);

/**
 * Terminate traversal over lattice nodes.
 */
void ms_latnode_iter_free(ms_latnode_iter_t *itor);

/**
 * Get number of entries.
 */
int ms_latnode_n_entries(ms_latnode_t *node);

/**
 * Get an entry link from a node.
 */
ms_latlink_t *ms_latnode_get_entry(ms_lattice_t *l,
                                   ms_latnode_t *node, int idx);

/**
 * Get an entry index from a node.
 */
int32 ms_latnode_get_entry_idx(ms_lattice_t *l,
                               ms_latnode_t *node, int idx);

/**
 * Get number of exits.
 */
int ms_latnode_n_exits(ms_latnode_t *node);

/**
 * Get an exit link from a node.
 */
ms_latlink_t *ms_latnode_get_exit(ms_lattice_t *l,
                                  ms_latnode_t *node, int idx);

/**
 * Delete a node.
 */
void ms_latnode_unlink(ms_lattice_t *l, ms_latnode_t *node);

/**
 * Delete a link.
 */
void ms_latlink_unlink(ms_lattice_t *l, ms_latlink_t *link);

/**
 * Perform N-Gram expansion on a lattice and assign language model
 * probabilities.
 */
int ms_lattice_expand(ms_lattice_t *l, ngram_model_t *lm);

/**
 * Run the forward algorithm on a lattice.
 *
 * FIXME: This will be done incrementally very soon (hopefully later tonight).
 */
int32 ms_lattice_forward(ms_lattice_t *l, int32 inv_aw);

/**
 * Run the backward algorithm on a lattice.
 *
 * FIXME: Partial version of this coming soon (hopefully later tonight).
 */
int32 ms_lattice_backward(ms_lattice_t *l, int32 inv_aw);

/**
 * Print a description of a lattice node.
 */
int ms_latnode_print(FILE *fh, ms_lattice_t *l, ms_latnode_t *n);

/**
 * Print a description of a lattice link.
 */
int ms_latlink_print(FILE *fh, ms_lattice_t *l, ms_latlink_t *vx);

#endif /* __MS_LATTICE_H__ */
