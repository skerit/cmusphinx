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
 * @file arc_buffer.c Queue passing hypotheses (arcs) between search passes
 */
#include "arc_buffer.h"

arc_buffer_t *
arc_buffer_init(void)
{
    arc_buffer_t *fab;

    fab = ckd_calloc(1, sizeof(*fab));
    fab->refcount = 1;
    fab->arcs = garray_init(0, sizeof(arc_t));
    fab->sf_idx = garray_init(0, sizeof(int));

    return fab;
}

arc_buffer_t *
arc_buffer_retain(arc_buffer_t *fab)
{
    ++fab->refcount;
    return fab;
}

int
arc_buffer_free(arc_buffer_t *fab)
{
    if (fab == NULL)
        return 0;
    if (--fab->refcount > 0)
        return fab->refcount;

    garray_free(fab->sf_idx);
    garray_free(fab->arcs);
    ckd_free(fab);
    return 0;
}

void
arc_buffer_dump(arc_buffer_t *fab, dict_t *dict)
{
    size_t i, n_arcs;

    n_arcs = garray_next_idx(fab->arcs);
    E_INFO("Arc buffer %p: %d arcs:\n", fab, n_arcs);
    for (i = garray_base(fab->arcs); i < n_arcs; ++i) {
        arc_t *arc = garray_ptr(fab->arcs, arc_t, i);
        E_INFO_NOFN("%s sf %d ef %d\n",
                    dict_wordstr(dict, arc->wid), arc->src, arc->dest);
    }
}

int
arc_buffer_extend(arc_buffer_t *fab, int next_sf)
{
    if (next_sf == fab->next_sf)
        return 0;
    garray_expand_to(fab->sf_idx, next_sf);
    fab->next_sf = next_sf;
    garray_clear(fab->sf_idx, fab->active_sf, fab->next_sf - fab->active_sf);
    return next_sf - fab->active_sf;
}

bpidx_t
arc_buffer_add_bps(arc_buffer_t *fab,
                   bptbl_t *bptbl, bpidx_t start,
                   bpidx_t end)
{
    bpidx_t idx, next_idx;
    int n_arcs;

    n_arcs = 0;
    next_idx = -1;
    for (idx = start; idx < end; ++idx) {
        arc_t arc;
        bp_t ent;

        /* Convert it to an arc. */
        bptbl_get_bp(bptbl, idx, &ent);
        arc.wid = ent.wid;
        arc.src = bptbl_sf(bptbl, idx);
        arc.dest = ent.frame;
        /* If it's inside the appropriate frame span, add it. */
        if (arc.src >= fab->active_sf && arc.src < fab->next_sf) {
            garray_append(fab->arcs, &arc);
            /* Increment the frame counter for its start frame. */
            ++garray_ent(fab->sf_idx, int, arc.src);
            ++n_arcs;
        }
        else {
            if (arc.src >= fab->active_sf && next_idx == -1)
                next_idx = idx;
        }
    }

    E_DEBUG(2,("Added bps from frame %d to %d, index %d to %d\n",
               fab->active_sf, fab->next_sf,
               start, end));
    if (next_idx == -1)
        next_idx = end;
    return next_idx;
}

int
arc_buffer_commit(arc_buffer_t *fab)
{
    size_t n_arcs;
    int n_active_fr, i, prev_count;
    garray_t *active_arc;
    garray_t *active_sf;
    int *sf;

    /* Save frame and arc counts. */
    n_active_fr = fab->next_sf - fab->active_sf;
    n_arcs = garray_next_idx(fab->arcs) - fab->active_arc;

    /* Nothing to do... */
    if (n_active_fr == 0) {
        assert(n_arcs == 0);
        return 0;
    }

    /* Sum forward frame counters to create arc indices. */
    sf = garray_ptr(fab->sf_idx, int, fab->active_sf);
    prev_count = sf[0];
    sf[0] = fab->active_arc;
    for (i = 1; i < n_active_fr; ++i)  {
        int tmp = sf[i];
        sf[i] = sf[i-1] + prev_count;
        prev_count = tmp;
    }

    if (n_arcs > 0) {
        /* Permute incoming arcs to match frame counters */
        active_sf = garray_slice(fab->sf_idx, fab->active_sf, n_active_fr);
        active_arc = garray_slice(fab->arcs, fab->active_arc, n_arcs);

        for (i = 0; i < n_arcs; ++i) {
            arc_t *arc = garray_ptr(active_arc, arc_t, i);
            int *pos = garray_ptr(active_sf, int, arc->src - fab->active_sf);
            /* Copy it into place. */
            garray_ent(fab->arcs, arc_t, *pos) = *arc;
            /* Increment local frame counter. */
            *pos += 1;
        }

        garray_free(active_sf);
        garray_free(active_arc);
    }

    /* Update frame and arc pointers. */
    fab->active_sf += n_active_fr;
    fab->active_arc += n_arcs;
    return n_arcs;
}

arc_t *
arc_buffer_iter(arc_buffer_t *fab, int sf)
{
    int idx;
    if (sf < garray_base(fab->sf_idx) || sf >= fab->active_sf)
        return NULL;
    idx = garray_ent(fab->sf_idx, int, sf);
    if (idx >= fab->active_arc)
        return NULL;
    return garray_ptr(fab->arcs, arc_t, idx);

}

arc_t *
arc_next(arc_buffer_t *fab, arc_t *ab)
{
    ab += 1;
    if (ab >= garray_ptr(fab->arcs, arc_t, fab->active_arc))
        return NULL;
    return ab;
}

arc_t *
arc_buffer_wait(arc_buffer_t *fab, int sf)
{
    /* FIXME: Implement this... */
    return NULL;
}

int
arc_buffer_release(arc_buffer_t *fab, int first_sf)
{
    int next_first_arc;
    if (first_sf == garray_base(fab->sf_idx))
        return 0;

    /* Get the new first arc. */
    next_first_arc = garray_ent(fab->sf_idx, int, first_sf);
    /* Shift back start frames and arcs. */
    garray_shift_from(fab->sf_idx, first_sf);
    garray_set_base(fab->sf_idx, first_sf);
    garray_shift_from(fab->arcs, next_first_arc);
    garray_set_base(fab->arcs, next_first_arc);

    return 0;
}

void
arc_buffer_reset(arc_buffer_t *fab)
{
    fab->active_sf = fab->next_sf = 0;
    fab->active_arc = 0;
    garray_reset(fab->arcs);
    garray_reset(fab->sf_idx);
}
