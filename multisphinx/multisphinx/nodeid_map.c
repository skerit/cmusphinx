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
 * @file nodeid_map.c Node ID maps.
 */

#include <sphinxbase/prim_type.h>
#include <sphinxbase/garray.h>
#include <sphinxbase/ckd_alloc.h>

#include "nodeid_map.h"

/**
 * Node ID mapping.
 *
 * For each frame we store a set of language model ID to node index
 * mappings, represented as i32p_t.  The number of language model IDs
 * in any given frame is usually fairly small so for the time being
 * this is simply stored as an unsorted list.
 */
struct nodeid_map_s {
    garray_t *frame_maps;
};

struct nodeid_iter_s {
    garray_t *frame_map;
    nodeid_map_t *nmap;
    int32 cf, pos;
};

nodeid_map_t *
nodeid_map_init(void)
{
    nodeid_map_t *nmap = ckd_calloc(1, sizeof(*nmap));
    nmap->frame_maps = garray_init(0, sizeof(garray_t *));
    return nmap;
}

int
nodeid_map_free(nodeid_map_t *nmap)
{
    int sf;
    for (sf = garray_base(nmap->frame_maps);
	 sf < garray_next_idx(nmap->frame_maps); ++sf)
	garray_free(garray_ent(nmap->frame_maps, garray_t *, sf));
    garray_free(nmap->frame_maps);
    ckd_free(nmap);
    return 0;
}

static garray_t *
nodeid_map_get_frame(nodeid_map_t *nmap, int32 sf)
{
    int32 next_sf = garray_next_idx(nmap->frame_maps);
    if (sf >= next_sf || sf < 0)
        return NULL;
    return garray_ent(nmap->frame_maps, garray_t *, sf);
}

static garray_t *
nodeid_map_push_frame(nodeid_map_t *nmap, int32 sf)
{
    int32 next_sf = garray_next_idx(nmap->frame_maps);
    if (sf >= next_sf) {
	garray_expand_to(nmap->frame_maps, sf + 1);
	garray_clear(nmap->frame_maps, next_sf, sf + 1 - next_sf);
    }
    if (garray_ent(nmap->frame_maps, garray_t *, sf) == NULL)
	garray_ent(nmap->frame_maps, garray_t *, sf)
	    = garray_init(0, sizeof(i32p_t));
    return nodeid_map_get_frame(nmap, sf);
}

static int32
frame_map_add(garray_t *frame_map, int32 lmstate, int32 idx)
{
    i32p_t map;
    int32 pos;

    map.a = lmstate;
    map.b = idx;
    for (pos = garray_base(frame_map);
	 pos < garray_next_idx(frame_map); ++pos) {
	i32p_t *m = garray_ptr(frame_map, i32p_t, pos);
	if (m->a == lmstate)
	    return pos;
    }
    garray_append(frame_map, &map);
    return pos;
}

int32
nodeid_map_add(nodeid_map_t *nmap, int sf, int32 lmstate, int32 idx)
{
    garray_t *frame_map;

    frame_map = nodeid_map_push_frame(nmap, sf);
    return frame_map_add(frame_map, lmstate, idx);
}

static int32
frame_map_remap(garray_t *frame_map, int32 lmstate, int32 idx)
{
    int32 pos;
    for (pos = garray_base(frame_map);
	 pos < garray_next_idx(frame_map); ++pos) {
	i32p_t *m = garray_ptr(frame_map, i32p_t, pos);
	if (m->a == lmstate) {
	    m->b = idx;
	    return pos;
	}
    }
    return -1;
}

int32
nodeid_map_remap(nodeid_map_t *nmap, int sf, int32 lmstate, int32 idx)
{
    garray_t *frame_map;

    frame_map = nodeid_map_get_frame(nmap, sf);
    if (frame_map == NULL)
	return -1;
    else
	return frame_map_remap(frame_map, lmstate, idx);
}

int32
nodeid_map_delete_frame(nodeid_map_t *nmap, int sf)
{
    garray_free(garray_ent(nmap->frame_maps, garray_t *, sf));
    garray_ent(nmap->frame_maps, garray_t *, sf) = NULL;
    return sf;
}

int32
frame_map_delete(garray_t *frame_map, int32 lmstate)
{
    int32 pos;
    for (pos = garray_base(frame_map);
	 pos < garray_next_idx(frame_map); ++pos) {
	i32p_t *m = garray_ptr(frame_map, i32p_t, pos);
	if (m->a == lmstate) {
	    garray_delete(frame_map, pos, 1);
	    return pos;
	}
    }
    return -1;
}

int32
nodeid_map_delete(nodeid_map_t *nmap, int sf, int32 lmstate)
{
    garray_t *frame_map;

    frame_map = nodeid_map_push_frame(nmap, sf);
    frame_map_delete(frame_map, lmstate);
    if (garray_size(frame_map) == 0)
	return nodeid_map_delete_frame(nmap, sf);
    else
	return garray_size(frame_map);
}

int32
frame_map_map(garray_t *frame_map, int32 lmstate)
{
    int32 pos;
    for (pos = garray_base(frame_map);
	 pos < garray_next_idx(frame_map); ++pos) {
	i32p_t *m = garray_ptr(frame_map, i32p_t, pos);
	if (m->a == lmstate) {
	    return m->b;
	}
    }
    return -1;
}

int32
nodeid_map_map(nodeid_map_t *nmap, int sf, int32 lmstate)
{
    garray_t *frame_map;

    frame_map = nodeid_map_get_frame(nmap, sf);
    if (frame_map == NULL)
	return -1;
    else
	return frame_map_map(frame_map, lmstate);
}

nodeid_iter_t *
nodeid_map_iter(nodeid_map_t *nmap, int sf)
{
    nodeid_iter_t *itor;
    garray_t *frame_map;

    if (sf == -1) {
	int cf;
	for (cf = garray_base(nmap->frame_maps);
	     cf < garray_next_idx(nmap->frame_maps); ++cf)
	    if ((frame_map = nodeid_map_get_frame(nmap, cf)) != NULL)
		break;
	if (cf == garray_next_idx(nmap->frame_maps))
	    return NULL;
	itor = ckd_calloc(1, sizeof(*itor));
	itor->frame_map = frame_map;
	itor->nmap = nmap;
	itor->cf = cf;
    }
    else {
	/* Need to get a specific frame, or else. */
	frame_map = nodeid_map_get_frame(nmap, sf);
	if (frame_map == NULL)
	    return NULL;
	itor = ckd_calloc(1, sizeof(*itor));
	itor->frame_map = frame_map;
	itor->nmap = NULL;
	itor->cf = sf;
    }
    itor->pos = 0;
    return itor;
}

nodeid_iter_t *
nodeid_iter_next(nodeid_iter_t *itor)
{
    ++itor->pos;
    if (itor->pos >= garray_next_idx(itor->frame_map)) {
	if (itor->nmap == NULL)
	    goto iter_done;
	else {
	    for (; itor->cf < garray_next_idx(itor->nmap->frame_maps);
		 ++itor->cf) {
		if ((itor->frame_map = nodeid_map_get_frame
		     (itor->nmap, itor->cf)) != NULL)
		    break;
	    }
	    if (itor->frame_map == NULL)
		goto iter_done;
	    itor->pos = 0;
	}
    }
    return itor;

iter_done:
    nodeid_iter_free(itor);
    return NULL;
}

int32
nodeid_iter_get(nodeid_iter_t *itor, int32 *out_lmstate)
{
    i32p_t *map = garray_ptr(itor->frame_map, i32p_t, itor->pos);
    if (out_lmstate)
	*out_lmstate = map->a;
    return map->b;
}

void
nodeid_iter_free(nodeid_iter_t *itor)
{
    ckd_free(itor);
}

