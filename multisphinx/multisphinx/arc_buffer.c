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

#include "bptbl.h"
#include "arc_buffer.h"

typedef enum arc_buffer_state_s {
    ARC_BUFFER_INITIAL = 0,
    ARC_BUFFER_RUNNING,
    ARC_BUFFER_FINAL,
    ARC_BUFFER_CANCELED
} arc_buffer_state_e;

struct arc_buffer_s {
    int refcount;
    char *name;
    char *uttid;
    /* FIXME: Probably non-orthogonal set of synchronization primitives */
    sbmtx_t *mtx;
    sbsem_t *start, *release;
    sbevent_t *evt;
    garray_t *arcs;
    garray_t *sf_idx;
    garray_t *rc_deltas;
    bptbl_t *input_bptbl;
    ngram_model_t *lm;
    rcdelta_t *tmp_rcdeltas;
    int max_n_rc;
    int state;
    int scores;
    int arc_size;
    int active_sf; /**< First frame of incoming arcs. */
    int next_sf;   /**< First frame not containing arcs (last frame + 1). */
    bpidx_t next_idx;  /**< Next bptbl index to scan from. */
    int active_arc; /**< First incoming arc. */
};

/* Internal functions for producers only. */
static int arc_buffer_extend(arc_buffer_t *fab, int next_sf);
static int arc_buffer_commit(arc_buffer_t *fab);

arc_buffer_t *
arc_buffer_init(char const *name, bptbl_t *input_bptbl,
                ngram_model_t *lm, int keep_scores)
{
    arc_buffer_t *fab;

    fab = ckd_calloc(1, sizeof(*fab));
    fab->refcount = 1;
    fab->name = ckd_salloc(name);
    fab->sf_idx = garray_init(0, sizeof(int));
    fab->start = sbsem_init("arc_buffer:start",0);
    fab->release = sbsem_init("arc_buffer:release",0);
    fab->evt = sbevent_init(FALSE);
    fab->mtx = sbmtx_init();
    fab->input_bptbl = bptbl_retain(input_bptbl);
    if (lm)
        fab->lm = ngram_model_retain(lm);
    fab->scores = keep_scores;
    if (keep_scores) {
        fab->rc_deltas = garray_init(0, sizeof(rcdelta_t));
        fab->max_n_rc = bin_mdef_n_ciphone(input_bptbl->d2p->mdef);
        fab->arc_size = sizeof(sarc_t) + sizeof(bitvec_t) * bitvec_size(fab->max_n_rc);
        fab->arcs = garray_init(0, fab->arc_size);
        fab->tmp_rcdeltas = ckd_calloc(fab->max_n_rc, sizeof(*fab->tmp_rcdeltas));
    }
    else {
        fab->arcs = garray_init(0, sizeof(arc_t));
        fab->arc_size = sizeof(arc_t);
    }
    E_INFO("Initialized arc buffer '%s', each arc occupies %d bytes\n",
           fab->name, fab->arc_size);

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
    garray_free(fab->rc_deltas);
    sbsem_free(fab->start);
    sbsem_free(fab->release);
    sbevent_free(fab->evt);
    sbmtx_free(fab->mtx);
    bptbl_free(fab->input_bptbl);
    ngram_model_free(fab->lm);
    ckd_free(fab->tmp_rcdeltas);
    ckd_free(fab->name);
    ckd_free(fab);
    return 0;
}

void
arc_buffer_lock(arc_buffer_t *fab)
{
    sbmtx_lock(fab->mtx);
}

void
arc_buffer_unlock(arc_buffer_t *fab)
{
    sbmtx_unlock(fab->mtx);
}

bptbl_t *
arc_buffer_input_bptbl(arc_buffer_t *fab)
{
    return fab->input_bptbl;
}

ngram_model_t *
arc_buffer_lm(arc_buffer_t *fab)
{
    return fab->lm;
}


void
arc_buffer_dump(arc_buffer_t *fab, dict_t *dict)
{
    size_t i, n_arcs;

    n_arcs = garray_next_idx(fab->arcs);
    E_INFO("Arc buffer '%s': %d arcs:\n", fab->name, n_arcs);
    for (i = garray_base(fab->arcs); i < n_arcs; ++i) {
        if (fab->scores) {
            sarc_t *arc = (sarc_t *)garray_void(fab->arcs, i);
            rcdelta_t const *d = arc_buffer_get_rcdeltas(fab, arc);
            int i;
            E_INFO_NOFN("%s %d %d %d %d",
                        dict_wordstr(dict, arc->arc.wid),
                        arc->arc.src, arc->arc.dest,
                        arc->score, arc->lscr);
            for (i = 0; i < arc_buffer_max_n_rc(fab); ++i) {
                if (bitvec_is_set(arc->rc_bits, i))
                    E_INFOCONT(" %d:%u", i, *d++);
            }
            E_INFOCONT("\n");
        }
        else {
            arc_t *arc = (arc_t *)garray_void(fab->arcs, i);
            E_INFO_NOFN("%s sf %d ef %d\n",
                        dict_wordstr(dict, arc->wid), arc->src, arc->dest);
        }
    }
}

int
arc_buffer_eou(arc_buffer_t *fab)
{
    return fab->state == ARC_BUFFER_FINAL;
}

int
arc_buffer_consumer_start_utt(arc_buffer_t *fab, int timeout)
{
    int s = (timeout == -1) ? -1 : 0;
    int rc;

    E_INFO("arc_buffer_consumer_start_utt\n");
    if ((rc = sbsem_down(fab->start, s, timeout)) < 0)
        return rc;
    if (fab->state == ARC_BUFFER_CANCELED)
        return -1;
    return 0;
}

int
arc_buffer_consumer_end_utt(arc_buffer_t *fab)
{
    return sbsem_up(fab->release);
}

void
arc_buffer_producer_start_utt(arc_buffer_t *fab, char *uttid)
{
    fab->active_sf = fab->next_sf = 0;
    fab->active_arc = 0;
    fab->next_idx = 0;
    fab->state = ARC_BUFFER_RUNNING;
    fab->uttid = uttid;
    garray_reset(fab->arcs);
    garray_reset(fab->sf_idx);
    /* If we ever have multiple consumers this will be sbsem_set() */
    E_INFO("arc_buffer_producer_start_utt\n");
    sbsem_up(fab->start);
}

static int
arc_buffer_extend(arc_buffer_t *fab, int next_sf)
{
    if (next_sf == fab->next_sf)
        return 0;
    garray_expand_to(fab->sf_idx, next_sf);
    fab->next_sf = next_sf;
    garray_clear(fab->sf_idx, fab->active_sf, fab->next_sf - fab->active_sf);
    return next_sf - fab->active_sf;
}

static bpidx_t
arc_buffer_add_bps(arc_buffer_t *fab,
                   bptbl_t *bptbl, bpidx_t start,
                   bpidx_t end)
{
    bpidx_t idx, next_idx;
    int n_arcs;

    n_arcs = 0;
    next_idx = -1;
    for (idx = start; idx < end; ++idx) {
        sarc_t sarc;
        bp_t ent;

        /* Convert it to an arc. */
        bptbl_get_bp(bptbl, idx, &ent);
        sarc.arc.wid = ent.wid;
        sarc.arc.src = bptbl_sf(bptbl, idx);
        sarc.arc.dest = ent.frame;
        /* If it's inside the appropriate frame span, add it. */
        if (sarc.arc.src >= fab->active_sf && sarc.arc.src < fab->next_sf) {
            /* Have to do this with a pointer because of the variable
             * length array rc_bits. */
            sarc_t *sp = garray_append(fab->arcs, &sarc);

            E_DEBUG(3,("Added arc %s %d -> %d\n",
                       dict_wordstr(bptbl->d2p->dict, sarc.arc.wid),
                       sarc.arc.src, sarc.arc.dest + 1));
            if (fab->scores) {
                /* Compress and add right context deltas. */
                int i, rcsize = bptbl_get_rcdeltas(bptbl, idx, fab->tmp_rcdeltas);
                int n_used;

                sp->score = ent.score;
                if (fab->lm)
                    sp->lscr = bptbl_fake_lmscore(bptbl, fab->lm, idx, &n_used);
                else
                    sp->lscr = 0;
                if (dict_filler_word(bptbl->d2p->dict, sp->arc.wid)
                    || sp->arc.wid == dict_startwid(bptbl->d2p->dict)) {
                    assert(sp->lscr == 0);
                    assert(rcsize == 1);
                }
                if (rcsize == 1)
                    assert(fab->tmp_rcdeltas[0] == 0);
                sp->rc_idx = garray_next_idx(fab->rc_deltas);
                bitvec_clear_all(sp->rc_bits, fab->max_n_rc);
                for (i = 0; i < rcsize; ++i) {
                    if (fab->tmp_rcdeltas[i] != NO_RC) {
                        bitvec_set(sp->rc_bits, i);
                        garray_append(fab->rc_deltas, &fab->tmp_rcdeltas[i]);
                    }
                }
            }
            /* Increment the frame counter for its start frame. */
            ++garray_ent(fab->sf_idx, int, sarc.arc.src);
            ++n_arcs;
        }
        else {
            /* Find the first index of an arc outside the span. */
            if (sarc.arc.src >= fab->active_sf && next_idx == -1)
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

int32
arc_buffer_get_rcscore(arc_buffer_t *fab, sarc_t *ab, int rc)
{
    /* FIXME: THIS IS WRONG!!! */
    return ab->score - garray_ent(fab->rc_deltas, rcdelta_t, ab->rc_idx + rc);
}

rcdelta_t const *
arc_buffer_get_rcdeltas(arc_buffer_t *fab, sarc_t *ab)
{
    return garray_ptr(fab->rc_deltas, rcdelta_t, ab->rc_idx);
}

int
arc_buffer_max_n_rc(arc_buffer_t *fab)
{
    return fab->max_n_rc;
}

bpidx_t
arc_buffer_producer_sweep(arc_buffer_t *fab, int release)
{
    int next_sf;

    arc_buffer_lock(fab);
    next_sf = bptbl_active_sf(fab->input_bptbl);
    if (arc_buffer_extend(fab, next_sf) > 0) {
        E_DEBUG(2,("%s adding arcs to frame %d idx %d:%d\n",
                   fab->name, next_sf, fab->next_idx,
                   bptbl_retired_idx(fab->input_bptbl)));
        fab->next_idx = arc_buffer_add_bps(fab, fab->input_bptbl,
                                           fab->next_idx,
                                           bptbl_retired_idx(fab->input_bptbl));
        if (release && fab->input_bptbl->oldest_bp > 0)
            bptbl_release(fab->input_bptbl, fab->input_bptbl->oldest_bp - 1);
        /* Do this after release since it may wake someone up. */
        arc_buffer_commit(fab);
    }
    arc_buffer_unlock(fab);
    return fab->next_idx;
}

bpidx_t
arc_buffer_producer_end_utt(arc_buffer_t *fab, int release)
{
    int next_sf, i, rc, nth;

    arc_buffer_lock(fab);
    next_sf = bptbl_active_sf(fab->input_bptbl);
    if (arc_buffer_extend(fab, next_sf) > 0) {
        fab->next_idx = arc_buffer_add_bps(fab, fab->input_bptbl,
                                           fab->next_idx,
                                           bptbl_retired_idx(fab->input_bptbl));
        if (release && fab->input_bptbl->oldest_bp > 0)
            bptbl_release(fab->input_bptbl, fab->input_bptbl->oldest_bp - 1);
        E_INFO("%s: marking arc buffer final\n", fab->name);
        fab->state = ARC_BUFFER_FINAL;
        /* Do this after marking the arc buffer final to avoid race conditions. */
        arc_buffer_commit(fab);
    }
    arc_buffer_unlock(fab);
    E_INFO("%s: allocated %d arcs (%d KiB)\n", fab->name,
           garray_alloc_size(fab->arcs),
           garray_alloc_size(fab->arcs) * fab->arc_size / 1024);
    E_INFO("%s: allocated %d start frame entries (%d KiB)\n", fab->name,
           garray_alloc_size(fab->sf_idx),
           garray_alloc_size(fab->sf_idx) * sizeof(int) / 1024);
    if (fab->rc_deltas)
        E_INFO("%s: allocated %d right context deltas (%d KiB)\n", fab->name,
               garray_alloc_size(fab->rc_deltas),
               garray_alloc_size(fab->rc_deltas) * sizeof(rcdelta_t) / 1024);

    nth = fab->refcount - 1;
    E_INFO("Waiting for %d consumers to finish\n", nth);
    for (i = 0; i < nth; ++i)
        if ((rc = sbsem_down(fab->release, -1, -1)) < 0)
            return rc;

    return 0;
}

static int
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
            arc_t *arc = (arc_t *)garray_void(active_arc, i);
            int *pos = garray_ptr(active_sf, int,
                                  arc->src - fab->active_sf);
            /* Copy it into place. */
            memcpy(garray_void(fab->arcs, *pos), arc, fab->arc_size);
            /* Increment local frame counter. */
            *pos += 1;
        }

        garray_free(active_sf);
        garray_free(active_arc);
    }

    /* Update frame and arc pointers. */
    fab->active_sf += n_active_fr;
    fab->active_arc += n_arcs;

    /* Signal consumer thread. */
    sbevent_signal(fab->evt);
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
arc_buffer_iter_next(arc_buffer_t *fab, arc_t *ab)
{
    ab = (arc_t *)((char *)ab + fab->arc_size);
    if (ab >= garray_ptr(fab->arcs, arc_t, fab->active_arc))
        return NULL;
    return ab;
}

int
arc_buffer_producer_shutdown(arc_buffer_t *fab)
{
    fab->state = ARC_BUFFER_CANCELED;
    return sbsem_up(fab->start);
}

int
arc_buffer_consumer_wait(arc_buffer_t *fab, int timeout)
{
    int sec = timeout / 1000000000;
    if (timeout == -1)
        sec = -1;
    if (sbevent_wait(fab->evt, sec, timeout) < 0)
        return -1;
    if (fab->state == ARC_BUFFER_CANCELED)
        return -1;
    return fab->next_sf;
}

int
arc_buffer_consumer_release(arc_buffer_t *fab, int first_sf)
{
    int next_first_arc;
    if (first_sf == garray_base(fab->sf_idx))
        return 0;

    arc_buffer_lock(fab);
    /* Get the new first arc. */
    next_first_arc = garray_ent(fab->sf_idx, int, first_sf);
    /* Shift back start frames and arcs. */
    garray_shift_from(fab->sf_idx, first_sf);
    garray_set_base(fab->sf_idx, first_sf);
    garray_shift_from(fab->arcs, next_first_arc);
    garray_set_base(fab->arcs, next_first_arc);
    /* FIXME: Have to also release rc deltas, but they are not in the
     * same order so this may not be entirely safe. */
    arc_buffer_unlock(fab);

    return 0;
}

char *
arc_buffer_uttid(arc_buffer_t *fab)
{
    return fab->uttid;
}
