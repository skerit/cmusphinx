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
 * @file arc_buffer.h Queue passing hypotheses (arcs) between search passes
 */

#ifndef __ARC_BUFFER_H__
#define __ARC_BUFFER_H__

/* SphinxBase headers. */
#include <sphinxbase/garray.h>

/* Local headers. */
#include "bptbl.h"

typedef struct arc_s {
    int32 wid;
    int32 score;
    int32 src;
    int32 dest;
} arc_t;

typedef struct arc_buffer_s {
    int refcount;
    garray_t *arcs;
    garray_t *sf_idx;
    int active_sf; /**< First frame of incoming arcs. */
    int next_sf;   /**< First frame not containing arcs (last frame + 1)
                      (FIXME: same as garray_next_ent(sf_idx)). */
    int active_arc; /**< First incoming arc. */
} arc_buffer_t;

/**
 * Create a new arc buffer.
 */
arc_buffer_t *arc_buffer_init(void);

/**
 * Retain a pointer to an arc buffer.
 */
arc_buffer_t *arc_buffer_retain(arc_buffer_t *fab);

/**
 * Release a pointer to an arc buffer.
 */
int arc_buffer_free(arc_buffer_t *fab);

/**
 * Clear the contents of an arc buffer.
 */
void arc_buffer_reset(arc_buffer_t *fab);

/**
 * Dump contents of arc buffer for debugging.
 */
void arc_buffer_dump(arc_buffer_t *fab, dict_t *dict);

/**
 * Extend the arc buffer up to the given frame index.
 *
 * @param next_sf Index of last start frame to be added to the buffer,
 * plus one.
 */
int arc_buffer_extend(arc_buffer_t *fab, int next_sf);

/**
 * Add arcs to the buffer from a bptbl.
 *
 * Only arcs starting in the newly extended frames will be
 * successfully added to the buffer.
 *
 * @param bptbl Backpointer table to take arcs from.
 * @param start First backpointer index to add to the buffer.
 * @param end One past the last backpointer index to add to the buffer.
 * @return The first backpointer index between start and end which
 *         starts after the next active frame, i.e. the first
 *         backpointer index which must be preserved for the next pass
 *         of arc addition.
 */
bpidx_t arc_buffer_add_bps(arc_buffer_t *fab,
                                   bptbl_t *bptbl, bpidx_t start,
                                   bpidx_t end);
/**
 * Commit extended arcs to the arc buffer.
 *
 * This freezes in place the start frames added since the last call to
 * arc_buffer_extend().  No more arcs with these start frames
 * may be added to the arc buffer.
 */
int arc_buffer_commit(arc_buffer_t *fab);

/**
 * Iterate over arcs in the arc buffer starting at given frame.
 *
 * @param sf Frame to iterate over.
 * @return First arc in frame, or NULL if frame not available.
 */
arc_t *arc_buffer_iter(arc_buffer_t *fab, int sf);

/**
 * Move the arc pointer forward.
 */
arc_t *arc_next(arc_buffer_t *fab, arc_t *ab);

/**
 * Wait until arcs for the given frame are committed.
 *
 * @return First arc in sf, or NULL if the utterance was terminated
 *         before sf was reached.
 */
arc_t *arc_buffer_wait(arc_buffer_t *fab, int sf);

/**
 * Release old arcs from the arc buffer.
 *
 * This releases all arcs starting in frames before first_sf.
 */
int arc_buffer_release(arc_buffer_t *fab, int first_sf);

#endif /* __ARC_BUFFER_H__ */
