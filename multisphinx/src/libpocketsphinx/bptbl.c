/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008-2010 Carnegie Mellon University.  All rights
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
 * @file bptbl.c Forward search lattice for N-Gram search.
 */

/* System headers. */
#include <string.h>
#include <assert.h>

/* SphinxBase headers. */
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/listelem_alloc.h>
#include <sphinxbase/err.h>

/* Local headers. */
#include "bptbl.h"

static int bptbl_rcsize(bptbl_t *bptbl, bp_t *be);

#if 1
#undef E_DEBUG
#define E_DEBUG(level,x) E_INFO x
#undef E_DEBUGCONT
#define E_DEBUGCONT(level,x) E_INFOCONT x
#endif

bptbl_t *
bptbl_init(dict2pid_t *d2p, int n_alloc, int n_frame_alloc)
{
    bptbl_t *bptbl = ckd_calloc(1, sizeof(*bptbl));

    bptbl->d2p = dict2pid_retain(d2p);
    bptbl->n_ent_alloc = n_alloc;
    bptbl->n_retired_alloc = n_alloc;
    bptbl->n_frame_alloc = n_frame_alloc;
    bptbl->n_permute_alloc = n_frame_alloc;

    bptbl->ent = ckd_calloc(bptbl->n_ent_alloc, sizeof(*bptbl->ent));
    bptbl->permute = ckd_calloc(bptbl->n_permute_alloc, sizeof(*bptbl->permute));

    bptbl->retired = ckd_calloc(bptbl->n_retired_alloc, sizeof(*bptbl->retired));
    bptbl->bscore_stack_size = bptbl->n_ent_alloc * 20;
    bptbl->bscore_stack = ckd_calloc(bptbl->bscore_stack_size,
                                     sizeof(*bptbl->bscore_stack));
    bptbl->ef_idx = ckd_calloc(bptbl->n_frame_alloc,
                               sizeof(*bptbl->ef_idx));
    bptbl->valid_fr = bitvec_alloc(bptbl->n_frame_alloc);

    return bptbl;
}

void
bptbl_free(bptbl_t *bptbl)
{
    if (bptbl == NULL)
        return;
    dict2pid_free(bptbl->d2p);
    ckd_free(bptbl->ent);
    ckd_free(bptbl->retired);
    ckd_free(bptbl->permute);
    ckd_free(bptbl->bscore_stack);
    ckd_free(bptbl->ef_idx);
    bitvec_free(bptbl->valid_fr);
    ckd_free(bptbl);
}

void
bptbl_reset(bptbl_t *bptbl)
{
    int i;

    for (i = 0; i < bptbl->n_frame_alloc; ++i) {
        bptbl->ef_idx[i] = -1;
    }
    bitvec_clear_all(bptbl->valid_fr, bptbl->n_frame_alloc);
    bptbl->first_invert_bp = 0;
    bptbl->dest_s_idx = 0;
    bptbl->n_frame = 0;
    bptbl->n_ent = 0;
    bptbl->n_retired = 0;
    bptbl->bss_head = 0;
    bptbl->active_fr = 0;
}

void
dump_bptable(bptbl_t *bptbl)
{
    int i;

    E_INFO("Retired backpointers (%d entries):\n", bptbl->first_invert_bp);
    for (i = 0; i < bptbl->first_invert_bp; ++i) {
        bp_t *ent = bptbl_ent(bptbl, i);
        assert(ent->valid);
        E_INFO_NOFN("%-5d %-10s start %-3d end %-3d score %-8d bp %-3d\n",
                    i, dict_wordstr(bptbl->d2p->dict, ent->wid),
                    bptbl_sf(bptbl, i),
                    ent->frame,
                    ent->score,
                    ent->bp);
    }
    E_INFO("Active backpointers (%d entries starting at %d):\n",
           bptbl->n_ent - bptbl->ef_idx[0], bptbl->ef_idx[0]);
    for (i = bptbl->ef_idx[0]; i < bptbl->n_ent; ++i) {
        bp_t *ent = bptbl_ent(bptbl, i);
        if (!ent->valid)
            E_INFO_NOFN("%-5d INVALID\n", i);
        else
            E_INFO_NOFN("%-5d %-10s start %-3d end %-3d score %-8d bp %-3d\n",
                        i, dict_wordstr(bptbl->d2p->dict, ent->wid),
                        bptbl_sf(bptbl, i),
                        ent->frame,
                        ent->score,
                        ent->bp);
    }
}

/**
 * Mark coaccessible active entries in the backpointer table.
 *
 * @param bptbl Backpointer table
 * @param ef Frame after last frame (as per ef_idx) to mark entries in
 * @param cf Current frame of search
 * @return Number of remaining active entries in marked region
 */
static int
bptbl_mark(bptbl_t *bptbl, int ef, int cf)
{
    int i, j, n_active_fr, last_gc_fr;

    assert(cf >= ef);
    assert(ef > bptbl->active_fr);

    /* Invalidate all active backpointer entries up to ef. */
    E_DEBUG(2,("Invalidating backpointers from %d to %d (%d to %d)\n",
               bptbl->active_fr, ef,
               bptbl_ef_idx(bptbl, bptbl->active_fr), bptbl_ef_idx(bptbl, ef)));
    for (i = bptbl_ef_idx(bptbl, bptbl->active_fr);
         i < bptbl_ef_idx(bptbl, ef); ++i) {
        E_DEBUG(5,("Invalidate bp %d\n", i));
        bptbl->ent[i - bptbl->ef_idx[0]].valid = FALSE;
    }

    /* Now re-activate all ones backwards reachable from the search graph. */
    E_DEBUG(2,("Finding coaccessible frames from backpointers from %d to %d (%d to %d)\n",
               ef, cf,
               bptbl_ef_idx(bptbl, ef), bptbl_ef_idx(bptbl, cf)));
    /* Mark everything immediately reachable from (ef..cf) */
    bitvec_clear_all(bptbl->valid_fr, cf - bptbl->active_fr);
    n_active_fr = 0;
    /* NOTE: This for statement can be sped up at the cost of being
     * less obvious. */
    for (i = bptbl_ef_idx(bptbl, ef);
         i < bptbl_ef_idx(bptbl, cf); ++i) {
        bp_t *ent, *prev;
        ent = bptbl_ent(bptbl, i);
        prev = bptbl_prev(bptbl, ent);
        assert(ent->valid);
        if (prev != NULL) {
            int frame = prev->frame;
            if (frame >= bptbl->active_fr
                && bitvec_is_clear(bptbl->valid_fr, frame - bptbl->active_fr)) {
                E_DEBUG(5,("Validate frame %d\n", frame - bptbl->active_fr));
                bitvec_set(bptbl->valid_fr, frame - bptbl->active_fr);
                ++n_active_fr;
            }
        }
    }
    /* Track the last frame with outgoing backpointers for gc */
    last_gc_fr = ef - 1;
    /* Walk back from every frame marked active up to the last one
     * marked active in the previous round (this means we scan some of
     * them multiple times, if this gets slow this can be
     * optimized). */
    while (n_active_fr > 0) {
        int next_gc_fr = 0;
        n_active_fr = 0;
        for (i = bptbl->active_fr; i <= last_gc_fr; ++i) {
            /* NOTE: We asserted i >= bptbl->active_fr */
            if (bitvec_is_set(bptbl->valid_fr, i - bptbl->active_fr)) {
                bitvec_clear(bptbl->valid_fr, i - bptbl->active_fr);
                /* Add all backpointers in this frame (the bogus
                 * lattice generation algorithm) */
                for (j = bptbl_ef_idx(bptbl, i);
                     j < bptbl_ef_idx(bptbl, i + 1); ++j) {
                    bp_t *ent = bptbl_ent(bptbl, j);
                    bp_t *prev = bptbl_prev(bptbl, ent);
                    ent->valid = TRUE;
                    if (prev != NULL) {
                        int frame = prev->frame;
                        if (frame >= bptbl->active_fr
                            && bitvec_is_clear(bptbl->valid_fr,
                                               frame - bptbl->active_fr)) {
                            bitvec_set(bptbl->valid_fr,
                                       frame - bptbl->active_fr);
                            ++n_active_fr;
                        }
                        if (frame > next_gc_fr)
                            next_gc_fr = frame;
                    }
                }
            }
        }
        E_DEBUG(3,("last_gc_fr %d => %d\n", last_gc_fr, next_gc_fr));
        last_gc_fr = next_gc_fr;
    }
    E_DEBUG(2,("Removed"));
    for (j = 0, i = bptbl_ef_idx(bptbl, bptbl->active_fr);
         i < bptbl_ef_idx(bptbl, ef); ++i) {
        bp_t *ent = bptbl_ent(bptbl, i);
        assert(ent != NULL);
        if (ent->valid)
            ++j;
        else
            E_DEBUGCONT(2,(" %d", i));
    }
    E_DEBUGCONT(2,("\n"));
    return j;
}

/**
 * Retire accessible backpointers
 *
 * @param bptbl Backpointer table.
 * @param n_bp Number of backpointers to be retired.
 * @param eidx Index of first backpointer which cannot be retired.
 * @return Index of last retired backpointer plus one.
 */
static int
bptbl_retire(bptbl_t *bptbl, int n_retired, int eidx)
{
    int src, dest;
    /* bss index for active frames (we don't want to track this in bptbl). */
    int active_dest_s_idx;

    /* First available backpointer index in retired. */
    dest = bptbl->first_invert_bp;
    /* Expand retired if necessary. */
    if (dest + n_retired > bptbl->n_retired_alloc) {
        while (dest + n_retired > bptbl->n_retired_alloc)
            bptbl->n_retired_alloc *= 2;
        assert(dest + n_retired <= bptbl->n_retired_alloc);
        bptbl->retired = ckd_realloc(bptbl->retired,
                                     bptbl->n_retired_alloc
                                     * sizeof(*bptbl->retired));
        E_INFO("Resized retired backpointer table to %d entries\n", bptbl->n_retired_alloc);
    }

    /* Note we use the "raw" backpointer indices here. */
    for (src = 0; src < eidx - bptbl->ef_idx[0]; ++src) {
        if (bptbl->ent[src].valid) {
            int rcsize = bptbl_rcsize(bptbl, bptbl->ent + src);
            E_DEBUG(4,("permute %d => %d\n", src + bptbl->ef_idx[0], dest));
            if (bptbl->ent[src].s_idx != bptbl->dest_s_idx) {
                E_DEBUG(4,("Moving %d rc scores from %d to %d for bptr %d\n",
                           rcsize, bptbl->ent[src].s_idx, bptbl->dest_s_idx, dest));
                assert(bptbl->ent[src].s_idx > bptbl->dest_s_idx);
                if (src < bptbl->n_ent - bptbl->ef_idx[0] - 1)
                    assert(bptbl->dest_s_idx + rcsize <= bptbl->ent[src + 1].s_idx);
                memmove(bptbl->bscore_stack + bptbl->dest_s_idx,
                        bptbl->bscore_stack + bptbl->ent[src].s_idx,
                        rcsize * sizeof(*bptbl->bscore_stack));
            }
            bptbl->retired[dest] = bptbl->ent[src];
            bptbl->retired[dest].s_idx = bptbl->dest_s_idx;

            assert(src < bptbl->n_permute_alloc);
            bptbl->permute[src] = dest;
            bptbl->dest_s_idx += rcsize;
            ++dest;
        }
        else {
            E_DEBUG(4,("permute %d => -1 src %d ef_idx %d\n",
                       src + bptbl->ef_idx[0], src, bptbl->ef_idx[0]));
            assert(src < bptbl->n_permute_alloc);
            bptbl->permute[src] = -1;
        }
    }
    /* We can keep compacting the bscore_stack since it is indirected. */
    if (src < bptbl->n_ent - bptbl->ef_idx[0]
        && bptbl->ent[src].s_idx != bptbl->dest_s_idx) {
        /* Leave dest_s_idx where it is for future compaction. */
        active_dest_s_idx = bptbl->dest_s_idx;
        while (src < bptbl->n_ent - bptbl->ef_idx[0]) {
            int rcsize = bptbl_rcsize(bptbl, bptbl->ent + src);
            E_DEBUG(4,("Moving %d rc scores from %d to %d for bptr %d\n",
                       rcsize, bptbl->ent[src].s_idx, active_dest_s_idx,
                       src + bptbl->ef_idx[0]));
            if (src < bptbl->n_ent - bptbl->ef_idx[0] - 1)
                assert(active_dest_s_idx + rcsize <= bptbl->ent[src + 1].s_idx);
            memmove(bptbl->bscore_stack + active_dest_s_idx,
                    bptbl->bscore_stack + bptbl->ent[src].s_idx,
                    rcsize * sizeof(*bptbl->bscore_stack));
            bptbl->ent[src].s_idx = active_dest_s_idx;
            active_dest_s_idx += rcsize;
            ++src;
        }
        bptbl->bss_head = active_dest_s_idx;
    }
    /* FIXME: This should work but I think it does not (i.e. we'd like
     * to be able to call this to sweep all active backpointers) */
    if (eidx == bptbl->n_ent)
        bptbl->n_ent = dest;
    return dest;
}


/**
 * Remap backpointers in backpointer table.
 *
 * @param bptbl Backpointer table.
 * @param last_retired_bp Index of last retired backpointer plus one
 * @param last_remapped_bp Last backpointer index to be remapped
 */
static void
bptbl_remap(bptbl_t *bptbl, int last_retired_bp, int last_remapped_bp)
{
    int i;

    E_DEBUG(2,("inverting %d:%d from %d to %d and %d to %d\n",
               bptbl->first_invert_bp, last_remapped_bp,
               bptbl->first_invert_bp, last_retired_bp,
               bptbl->ef_idx[0], bptbl->n_ent));
    /* First remap backpointers in newly retired bps. */
    for (i = bptbl->first_invert_bp; i < last_retired_bp; ++i) {
        /* Remember, these are the *source* backpointer indices, so
         * they fall in the range between prev_active_fr (which is the
         * first index of ef_idx) and active_fr. */
        if (bptbl->retired[i].bp >= bptbl->ef_idx[0]
            && bptbl->retired[i].bp < last_remapped_bp) {
            assert(bptbl->retired[i].bp - bptbl->ef_idx[0] < bptbl->n_permute_alloc);
            if (bptbl->retired[i].bp
                != bptbl->permute[bptbl->retired[i].bp - bptbl->ef_idx[0]])
                E_DEBUG(4,("remap retired %d => %d in %d\n",
                           bptbl->retired[i].bp,
                           bptbl->permute[bptbl->retired[i].bp - bptbl->ef_idx[0]], i));
            bptbl->retired[i].bp
                = bptbl->permute[bptbl->retired[i].bp - bptbl->ef_idx[0]];
            E_INFO("sf %d frame %d\n", bptbl_sf(bptbl, i), bptbl->retired[i].frame);
            assert(bptbl_sf(bptbl, i) <= bptbl->retired[i].frame);
        }
    }
    /* Now remap backpointers in still-active bps (which point to the
     * newly retired ones) */
    if (last_retired_bp == bptbl->n_ent)
        return;
    for (i = 0; i < bptbl->n_ent - bptbl->ef_idx[0]; ++i) {
        if (bptbl->ent[i].bp >= bptbl->ef_idx[0]
            && bptbl->ent[i].bp < last_remapped_bp) {
            assert(bptbl->ent[i].bp - bptbl->ef_idx[0] < bptbl->n_permute_alloc);
            if (bptbl->ent[i].bp
                != bptbl->permute[bptbl->ent[i].bp - bptbl->ef_idx[0]])
                E_DEBUG(4,("remap active %d => %d in %d\n",
                           bptbl->ent[i].bp,
                           bptbl->permute[bptbl->ent[i].bp - bptbl->ef_idx[0]],
                           i + bptbl->ef_idx[0]));
            bptbl->ent[i].bp
                = bptbl->permute[bptbl->ent[i].bp - bptbl->ef_idx[0]];
            assert(bptbl_sf(bptbl, i + bptbl->ef_idx[0]) < bptbl->ent[i].frame);
        }
    }
}

/**
 * Update the active frame pointer and backpointer array.
 */
static void
bptbl_update_active(bptbl_t *bptbl, int active_fr, int last_retired_bp)
{
    int frame_delta = active_fr - bptbl->active_fr;
    int bp_delta;
    /* This means nothing happened. */
    if (frame_delta == 0)
        return;
    assert(frame_delta > 0);
    bp_delta = bptbl->ef_idx[frame_delta] - bptbl->ef_idx[0];

    /* Push back active backpointers (eventually this will be circular) */
    E_DEBUG(3,("moving %d ent from %d (%d - %d)\n",
               bptbl->n_ent - bptbl->ef_idx[frame_delta],
               bp_delta, bptbl->ef_idx[frame_delta], bptbl->ef_idx[0]));
    memmove(bptbl->ent, bptbl->ent + bp_delta,
            (bptbl->n_ent - bptbl->ef_idx[frame_delta]) * sizeof(*bptbl->ent));
    
    /* Update ef_idx (implicitly updating ef_idx[0] */
    E_DEBUG(3,("moving %d ef_idx from %d (%d - %d)\n",
               bptbl->n_frame - active_fr,
               frame_delta, active_fr, bptbl->active_fr));
    assert(frame_delta + bptbl->n_frame - active_fr < bptbl->n_frame_alloc);
    memmove(bptbl->ef_idx, bptbl->ef_idx + frame_delta,
            (bptbl->n_frame - active_fr) * sizeof(*bptbl->ef_idx));

    /* And now update stuff. */
    bptbl->active_fr = active_fr;
    bptbl->first_invert_bp = last_retired_bp;
}

static void
bptbl_gc(bptbl_t *bptbl, int oldest_bp, int frame_idx)
{
    int prev_active_fr, active_fr;
    int n_retired, last_retired_bp;

    /* active_fr is the first frame which is still active in search
     * (i.e. for which outgoing word arcs can still be generated).
     * Therefore, any future backpointer table entries will not point
     * backwards to any backpointers before active_fr, and thus any
     * backpointers which are not reachable from those exiting in
     * active_fr will never be reachable. */
    prev_active_fr = bptbl->active_fr;
    if (oldest_bp == NO_BP)
        active_fr = 0;
    else
        active_fr = bptbl_ent(bptbl, oldest_bp)->frame;
    assert(active_fr >= prev_active_fr);
    /* Need at least 2 frames to GC. */
    if (active_fr <= prev_active_fr + 1)
        return;
    /* If there is nothing to GC then finish up. */
    if (bptbl_ef_idx(bptbl, prev_active_fr)
        == bptbl_ef_idx(bptbl, active_fr)) {
        bptbl_update_active(bptbl, active_fr, bptbl->first_invert_bp);
        return;
    }
    E_DEBUG(2,("GC from frame %d to %d\n", prev_active_fr, active_fr));
    /* Expand the permutation table if necessary. */
    if ((bptbl_ef_idx(bptbl, active_fr) - bptbl->ef_idx[0])
        > bptbl->n_permute_alloc) {
        while (bptbl_ef_idx(bptbl, active_fr) - bptbl->ef_idx[0]
               > bptbl->n_permute_alloc)
            bptbl->n_permute_alloc *= 2;
        bptbl->permute = ckd_realloc(bptbl->permute,
                                     bptbl->n_permute_alloc
                                     * sizeof(*bptbl->permute));
        E_INFO("Resized permutation table to %d entries (active ef = %d ef_idx[0] = %d)\n",
               bptbl->n_permute_alloc, bptbl_ef_idx(bptbl, active_fr), bptbl->ef_idx[0]);
    }
    /* Mark, compact, snap pointers. */
    n_retired = bptbl_mark(bptbl, active_fr, frame_idx);
    E_DEBUG(2,("About to retire %d bps\n", n_retired));
    last_retired_bp = bptbl_retire(bptbl, n_retired, bptbl_ef_idx(bptbl, active_fr));
    assert(n_retired == last_retired_bp - bptbl->first_invert_bp);
    bptbl_remap(bptbl, last_retired_bp, bptbl_ef_idx(bptbl, active_fr));
    bptbl_update_active(bptbl, active_fr, last_retired_bp);
    E_INFO("Retired %d bps: now %d retired, %d active\n", n_retired,
           bptbl->first_invert_bp, bptbl->n_ent - bptbl->ef_idx[0]);
}

int
bptbl_push_frame(bptbl_t *bptbl, int oldest_bp)
{
    int frame_idx = bptbl->n_frame;

    E_DEBUG(2,("pushing frame %d, oldest bp %d in frame %d\n",
               frame_idx, oldest_bp,
               oldest_bp == NO_BP
               ? -1 : bptbl_ent(bptbl, oldest_bp)->frame));
    if (frame_idx - bptbl->active_fr >= bptbl->n_frame_alloc) {
        bptbl->n_frame_alloc *= 2;
        E_INFO("Reallocating frame-based bptr arrays to %d\n", bptbl->n_frame_alloc);
        bptbl->ef_idx = ckd_realloc(bptbl->ef_idx,
                                    bptbl->n_frame_alloc * sizeof(*bptbl->ef_idx));
        bptbl->valid_fr = bitvec_realloc(bptbl->valid_fr, bptbl->n_frame_alloc);
    }
    bptbl->ef_idx[frame_idx - bptbl->active_fr] = bptbl->n_ent;
    bptbl->n_frame = frame_idx + 1;
    bptbl_gc(bptbl, oldest_bp, frame_idx);
    return frame_idx;
}

int
bptbl_finalize(bptbl_t *bptbl)
{
    int n_retired, last_retired_bp;

    E_DEBUG(2,("Final GC from frame %d to %d\n",
               bptbl->active_fr, bptbl->n_frame));
    /* If there is nothing to GC then finish up. */
    if (bptbl->n_ent == bptbl->ef_idx[0])
        return 0;
    /* Mark and GC everything from the last frame. */
    n_retired = bptbl_mark(bptbl, bptbl->n_frame - 1, bptbl->n_frame);
    E_DEBUG(2,("About to retire %d bps\n", n_retired));
    last_retired_bp = bptbl_retire(bptbl, n_retired, bptbl->n_ent);
    /* Last retired bp should be the same as bptbl->n_ent */
    bptbl_remap(bptbl, last_retired_bp, bptbl->n_ent);
    /* Just invalidate active entries, no need to move anything. */
    bptbl->first_invert_bp = last_retired_bp;
    bptbl->active_fr = bptbl->n_frame;
    bptbl->ef_idx[0] = bptbl->n_ent;
    E_INFO("Retired %d bps: now %d retired, %d active\n", n_retired,
           bptbl->first_invert_bp, bptbl->n_ent - bptbl->ef_idx[0]);
    return n_retired;
}


int32
bptbl_ef_idx(bptbl_t *bptbl, int frame_idx)
{
    if (frame_idx < bptbl->active_fr)
        return 0;
    else if (frame_idx >= bptbl->n_frame)
        return bptbl->n_ent;
    else {
        return bptbl->ef_idx[frame_idx - bptbl->active_fr];
    }
}

bp_t *
bptbl_ent(bptbl_t *bptbl, bpidx_t bpidx)
{
    if (bpidx == NO_BP)
        return NULL;
    if (bpidx < bptbl->ef_idx[0])
        return bptbl->retired + bpidx;
    else
        return bptbl->ent + (bpidx - bptbl->ef_idx[0]);
}

bpidx_t
bptbl_idx(bptbl_t *bptbl, bp_t *bpe)
{
    if (bpe->frame < bptbl->active_fr)
        return bpe - bptbl->retired;
    else
        return (bpe - bptbl->ent) + bptbl->ef_idx[0];
}

bp_t *
bptbl_prev(bptbl_t *bptbl, bp_t *ent)
{
    return bptbl_ent(bptbl, ent->bp);
}

int
bptbl_sf(bptbl_t *bptbl, bpidx_t bpidx)
{
    bp_t *ent = bptbl_ent(bptbl, bpidx);
    bp_t *prev;

    if (ent == NULL)
        return -1;
    prev = bptbl_prev(bptbl, ent);
    if (prev == NULL)
        return 0;
    else
        return prev->frame + 1;
}

int
bptbl_ef_count(bptbl_t *bptbl, int frame_idx)
{
    return bptbl_ef_idx(bptbl, frame_idx + 1)
        - bptbl_ef_idx(bptbl, frame_idx);
}

static int
bptbl_rcsize(bptbl_t *bptbl, bp_t *be)
{
    int rcsize;

    if (dict_is_single_phone(bptbl->d2p->dict, be->wid)) {
        be->last2_phone = -1;
        rcsize = 1;
    }
    else {
        be->last2_phone = dict_second_last_phone(bptbl->d2p->dict, be->wid);
        rcsize = dict2pid_rssid(bptbl->d2p, be->last_phone, be->last2_phone)->n_ssid;
    }

    return rcsize;
}

bp_t *
bptbl_enter(bptbl_t *bptbl, int32 w, int32 path, int32 score, int rc)
{
    int32 i, rcsize, *bss;
    bp_t *be;

    /* This might happen if recognition fails. */
    if (bptbl->n_ent == NO_BP) {
        E_ERROR("No entries in backpointer table!");
        return NULL;
    }

    /* Expand the backpointer tables if necessary. */
    if (bptbl->n_ent - bptbl->ef_idx[0] >= bptbl->n_ent_alloc) {
        bptbl->n_ent_alloc *= 2;
        assert(bptbl->n_ent - bptbl->ef_idx[0] < bptbl->n_ent_alloc);
        bptbl->ent = ckd_realloc(bptbl->ent,
                                 bptbl->n_ent_alloc
                                 * sizeof(*bptbl->ent));
        E_INFO("Resized backpointer table to %d entries\n", bptbl->n_ent_alloc);
    }
    if (bptbl->bss_head >= bptbl->bscore_stack_size
        - bin_mdef_n_ciphone(bptbl->d2p->mdef)) {
        bptbl->bscore_stack_size *= 2;
        bptbl->bscore_stack = ckd_realloc(bptbl->bscore_stack,
                                          bptbl->bscore_stack_size
                                          * sizeof(*bptbl->bscore_stack));
        E_INFO("Resized score stack to %d entries\n", bptbl->bscore_stack_size);
    }

    be = bptbl_ent(bptbl, bptbl->n_ent);
    be->wid = w;
    be->frame = bptbl->n_frame - 1;
    be->bp = path;
    be->score = score;
    be->s_idx = bptbl->bss_head;
    be->valid = TRUE;
    be->last_phone = dict_last_phone(bptbl->d2p->dict,w);
    bptbl_fake_lmstate(bptbl, bptbl->n_ent);

    /* DICT2PID */
    /* Get diphone ID for final phone and number of ssids corresponding to it. */
    rcsize = bptbl_rcsize(bptbl, be);
    /* Allocate some space on the bptbl->bscore_stack for all of these triphones. */
    for (i = rcsize, bss = bptbl->bscore_stack + bptbl->bss_head; i > 0; --i, bss++)
        *bss = WORST_SCORE;
    bptbl->bscore_stack[bptbl->bss_head + rc] = score;
    E_DEBUG(3,("Entered bp %d sf %d ef %d active_fr %d\n", bptbl->n_ent,
               bptbl_sf(bptbl, bptbl->n_ent), bptbl->n_frame - 1, bptbl->active_fr));
    assert(bptbl_sf(bptbl, bptbl->n_ent) >= bptbl->active_fr);

    bptbl->n_ent++;
    bptbl->bss_head += rcsize;

    return be;
}

void
bptbl_fake_lmstate(bptbl_t *bptbl, int32 bp)
{
    bp_t *ent, *prev;

    assert(bp != NO_BP);
    ent = bptbl_ent(bptbl, bp);
    prev = bptbl_prev(bptbl, ent);
    /* Propagate lm state for fillers, rotate it for words. */
    if (dict_filler_word(bptbl->d2p->dict, ent->wid)) {
        if (prev != NULL) {
            ent->real_wid = prev->real_wid;
            ent->prev_real_wid = prev->prev_real_wid;
        }
    }
    else {
        ent->real_wid = dict_basewid(bptbl->d2p->dict, ent->wid);
        if (prev != NULL)
            ent->prev_real_wid = prev->real_wid;
    }
}
