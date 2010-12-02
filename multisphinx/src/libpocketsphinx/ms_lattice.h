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

#include <multisphinx/dict.h>

#include "nodeid_map.h"

/**
 * Word lattice.
 */
typedef struct ms_lattice_s ms_lattice_t;

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

/**
 * Create a node.
 */
ms_latnode_t *ms_lattice_node_init(ms_lattice_t *l, int sf, int32 lmstate);

/**
 * Get a node by index.
 */
ms_latnode_t *ms_lattice_get_node_idx(ms_lattice_t *l, int32 idx);

/**
 * Get a node by ID.
 */
ms_latnode_t *ms_lattice_get_node_id(ms_lattice_t *l, int sf, int32 lmstate);

/**
 * Get the index of a node.
 */
int32 ms_lattice_get_idx_node(ms_lattice_t *l, ms_latnode_t *node);

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

#endif /* __MS_LATTICE_H__ */
